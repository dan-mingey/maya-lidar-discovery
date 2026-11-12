#!/usr/bin/env python3
"""Rigorous, skeptical structural-vs-karst forensic audit of Candidate #08
(positive control = Nuevo Jalisco) and all 14 surrounding Lacandon complexes
(LC-01 .. LC-14) using 0.33 m DEM, 0.33 m G1, and raw .LAZ ground returns.

Computes objective discriminators separating Classic Maya orthogonal plazas &
closed 360-degree pyramidal/platform mounds from folded Chiapas limestone
bedding planes and karst solution ridges:
  1. True 8-azimuth (360-deg) topographic drop on raw 0.33 m DEM & LRM
  2. Mound rectangularity, flat-top ratio, and aspect ratio
  3. 2D gradient structure-tensor anisotropy (1D parallel strike ridge vs 2D orthogonal walls)
  4. Courtyard floor flatness (area with true DEM slope < 2.0 deg and < 3.0 deg)
  5. 4-quadrant orthogonal enclosure by closed mounds
  6. Raw .LAZ ground-return density strictly on mounds and strictly on the plaza floor
  7. Longitudinal & transverse raw .LAZ + 0.33 m DEM elevation cross-sections
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi

WORK_DIR = Path("/agent/lidar/work/lacandon_deep_dive")
DATA_OUT = Path("/agent/maya-lidar-discovery/data/open_lidar_hunt_v2/lacandon_confidence_audit.json")


def compute_slope_deg(dem: np.ndarray, res_m: float, smooth_m: float = 2.0) -> np.ndarray:
    valid = np.isfinite(dem)
    med = float(np.nanmedian(dem[valid])) if valid.any() else 0.0
    filled = np.where(valid, dem, med)
    sig = max(smooth_m / res_m, 0.8)
    sm = ndi.gaussian_filter(filled, sigma=sig)
    gy, gx = np.gradient(sm, res_m, res_m)
    slp = np.degrees(np.arctan(np.hypot(gx, gy)))
    slp[~valid] = np.nan
    return slp


def compute_gradient_anisotropy(lrm: np.ndarray, res_m: float) -> float:
    """Ratio lambda_1 / lambda_2 of the 2D gradient structure tensor inside central 200x200m.
    1.0 = perfectly balanced 2D orthogonal walls; >2.5 = 1D linear parallel strike ridges.
    """
    H, W = lrm.shape
    r0 = max(0, int(H * 0.13))
    r1 = min(H, int(H * 0.87))
    sub = np.nan_to_num(lrm[r0:r1, r0:r1], nan=0.0)
    sm = ndi.gaussian_filter(sub, sigma=1.2 / res_m)
    gy, gx = np.gradient(sm, res_m, res_m)
    mag = np.hypot(gx, gy)
    thr = np.percentile(mag, 65)
    m = mag >= max(thr, 0.05)
    if m.sum() < 50:
        return 1.5
    Jxx = float(np.mean(gx[m] ** 2))
    Jyy = float(np.mean(gy[m] ** 2))
    Jxy = float(np.mean(gx[m] * gy[m]))
    evals = np.linalg.eigvalsh(np.array([[Jxx, Jxy], [Jxy, Jyy]]))
    l1, l2 = float(max(evals[1], 1e-6)), float(max(evals[0], 1e-6))
    return l1 / l2


def audit_single_complex(
    code: str,
    dem: np.ndarray,
    lrm: np.ndarray,
    g1: np.ndarray,
    px: np.ndarray,
    py: np.ndarray,
    pz: np.ndarray,
    pcls: np.ndarray,
    res_m: float,
    hm: float,
    meta: dict,
) -> dict:
    H, W = dem.shape
    valid = np.isfinite(dem) & np.isfinite(lrm)
    slp_deg = compute_slope_deg(dem, res_m=res_m, smooth_m=2.0)

    # Map raw .laz points to pixel grid
    p_col = np.clip(((px + hm) / (2 * hm) * W).astype(int), 0, W - 1)
    p_row = np.clip(((hm - py) / (2 * hm) * H).astype(int), 0, H - 1)
    pt_grid = np.zeros((H, W), dtype=np.int32)
    np.add.at(pt_grid, (p_row, p_col), 1)
    dist_pt_m = ndi.distance_transform_edt(pt_grid == 0) * res_m

    yy, xx = np.mgrid[0:H, 0:W]
    xm = (xx + 0.5) / W * (2 * hm) - hm
    ym = hm - (yy + 0.5) / H * (2 * hm)

    # Segment candidate mounds inside central 220x220m (|x|,|y| <= 110m)
    in_220 = (np.abs(xm) <= 110.0) & (np.abs(ym) <= 110.0)
    m_mask = valid & in_220 & (lrm >= 1.20)
    op_iter = max(1, int(round(0.7 / res_m)))
    m_mask = ndi.binary_opening(m_mask, iterations=op_iter)
    labeled, nlab = ndi.label(m_mask)

    ray_angles = [math.radians(a) for a in range(0, 360, 45)]
    mounds_detail = []

    for lab in range(1, nlab + 1):
        comp = labeled == lab
        area_m2 = float(comp.sum() * res_m * res_m)
        if area_m2 < 28.0 or area_m2 > 9500.0:
            continue
        ys, xs = np.where(comp)
        pts_xy = np.column_stack([xm[ys, xs], ym[ys, xs]])
        mean_xy = pts_xy.mean(axis=0)
        X = pts_xy - mean_xy
        cov = (X.T @ X) / max(len(X) - 1, 1)
        evals, evecs = np.linalg.eigh(cov)
        order = np.argsort(evals)[::-1]
        evecs = evecs[:, order]
        rot = X @ evecs
        len_maj = float(max(rot[:, 0].max() - rot[:, 0].min() + res_m, res_m))
        len_min = float(max(rot[:, 1].max() - rot[:, 1].min() + res_m, res_m))
        aspect = len_maj / max(len_min, 0.5)
        rect_fill = min(1.0, area_m2 / (len_maj * len_min))
        az_deg = float(math.degrees(math.atan2(evecs[1, 0], evecs[0, 0])) % 180.0)

        # Peak pixel inside mound
        lrm_vals = lrm[ys, xs]
        pk_idx = int(np.nanargmax(lrm_vals))
        pk_r, pk_c = int(ys[pk_idx]), int(xs[pk_idx])
        pk_x, pk_y = float(xm[pk_r, pk_c]), float(ym[pk_r, pk_c])
        pk_dem = float(dem[pk_r, pk_c])
        pk_lrm = float(lrm[pk_r, pk_c])

        # 8-azimuth 360-degree drop test out to max(18m, 0.75*len_maj + 8m), capped at 38m
        max_r_m = min(38.0, max(18.0, 0.65 * len_maj + 8.0))
        radii = np.linspace(4.0, max_r_m, 28)
        dem_drops = []
        lrm_drops = []
        for ang in ray_angles:
            rx = pk_x + radii * math.cos(ang)
            ry = pk_y + radii * math.sin(ang)
            rc = np.clip(((rx + hm) / (2 * hm) * W).astype(int), 0, W - 1)
            rr = np.clip(((hm - ry) / (2 * hm) * H).astype(int), 0, H - 1)
            d_ray = dem[rr, rc]
            l_ray = lrm[rr, rc]
            vd = np.isfinite(d_ray)
            if vd.any():
                dem_drops.append(float(pk_dem - np.nanmin(d_ray[vd])))
                lrm_drops.append(float(pk_lrm - np.nanmin(l_ray[vd])))
            else:
                dem_drops.append(0.0)
                lrm_drops.append(0.0)

        closed_rays_dem = sum(1 for d in dem_drops if d >= 0.75)
        closed_rays_lrm = sum(1 for d in lrm_drops if d >= 0.90)
        # Opposite-pair check along major axis: does a strike ridge continue without dropping?
        min_dem_drop = float(min(dem_drops))
        p25_dem_drop = float(np.percentile(dem_drops, 25))

        # Summit flatness: median slope of top 30% of mound pixels
        top30 = lrm_vals >= np.nanpercentile(lrm_vals, 70)
        summit_slp = float(np.nanmedian(slp_deg[ys[top30], xs[top30]]))

        laz_on = int(pt_grid[comp].sum())
        dil = ndi.binary_dilation(comp, iterations=max(1, int(round(2.0 / res_m))))
        laz_collar = int(pt_grid[dil].sum())

        is_closed_360 = (closed_rays_dem >= 7 and closed_rays_lrm >= 7 and aspect <= 5.0)
        is_ortho_rect = is_closed_360 and (rect_fill >= 0.56) and (aspect <= 4.0)

        mounds_detail.append(
            {
                "cx_m": round(float(mean_xy[0]), 1),
                "cy_m": round(float(mean_xy[1]), 1),
                "pk_x_m": round(pk_x, 1),
                "pk_y_m": round(pk_y, 1),
                "area_m2": round(area_m2, 1),
                "len_maj_m": round(len_maj, 1),
                "len_min_m": round(len_min, 1),
                "aspect": round(aspect, 2),
                "rect_fill": round(rect_fill, 3),
                "az_deg": round(az_deg, 1),
                "max_lrm_m": round(pk_lrm, 2),
                "pk_dem_m": round(pk_dem, 2),
                "closed_rays_dem_of_8": closed_rays_dem,
                "closed_rays_lrm_of_8": closed_rays_lrm,
                "min_dem_drop_m": round(min_dem_drop, 2),
                "p25_dem_drop_m": round(p25_dem_drop, 2),
                "summit_med_slope_deg": round(summit_slp, 2),
                "laz_pts_on": laz_on,
                "laz_pts_collar": laz_collar,
                "laz_pts_m2": round(laz_on / area_m2, 3),
                "is_closed_360": bool(is_closed_360),
                "is_ortho_rect": bool(is_ortho_rect),
            }
        )

    mounds_detail.sort(key=lambda m: m["max_lrm_m"] * math.sqrt(m["area_m2"]), reverse=True)

    # Courtyard / Plaza flatness inside inner 140x140m (|x|,|y| <= 70m, |LRM| <= 0.85m)
    in_140 = valid & (np.abs(xm) <= 70.0) & (np.abs(ym) <= 70.0)
    non_mound_140 = in_140 & (np.abs(lrm) <= 0.85)
    plaza_lt2_mask = non_mound_140 & (slp_deg <= 2.0)
    plaza_lt3_mask = non_mound_140 & (slp_deg <= 3.0)
    plaza_lt4p5_mask = non_mound_140 & (slp_deg <= 4.5)

    plaza_lt2_m2 = float(plaza_lt2_mask.sum() * res_m * res_m)
    plaza_lt3_m2 = float(plaza_lt3_mask.sum() * res_m * res_m)
    plaza_lt4p5_m2 = float(plaza_lt4p5_mask.sum() * res_m * res_m)
    plaza_med_slope = float(np.nanmedian(slp_deg[non_mound_140])) if non_mound_140.any() else 99.0
    bg_med_slope = float(np.nanmedian(slp_deg[valid])) if valid.any() else 99.0

    plaza_laz_pts = int(pt_grid[plaza_lt4p5_mask].sum())
    plaza_pts_m2 = plaza_laz_pts / max(plaza_lt4p5_m2, 1.0)

    # Principal architectural orientation & 4-quadrant enclosure by CLOSED 360° mounds
    closed_m = [m for m in mounds_detail if m["is_closed_360"]]
    ortho_m = [m for m in mounds_detail if m["is_ortho_rect"]]

    use_m = closed_m[:8] if len(closed_m) >= 2 else mounds_detail[:8]
    if use_m:
        angs4 = np.radians([ (m["az_deg"] % 90.0) * 4.0 for m in use_m ])
        dom_az_mod90 = (math.degrees(math.atan2(np.mean(np.sin(angs4)), np.mean(np.cos(angs4)))) / 4.0) % 90.0
        ortho_coh = float(math.hypot(np.mean(np.cos(angs4)), np.mean(np.sin(angs4))))
        cx0 = float(np.mean([m["cx_m"] for m in use_m[:min(5, len(use_m))]]))
        cy0 = float(np.mean([m["cy_m"] for m in use_m[:min(5, len(use_m))]]))
    else:
        dom_az_mod90, ortho_coh, cx0, cy0 = 0.0, 0.0, 0.0, 0.0

    ca, sa = math.cos(math.radians(dom_az_mod90)), math.sin(math.radians(dom_az_mod90))
    quads_closed = set()
    for m in closed_m[:10]:
        dx_q, dy_q = m["cx_m"] - cx0, m["cy_m"] - cy0
        u = dx_q * ca + dy_q * sa
        v = -dx_q * sa + dy_q * ca
        quads_closed.add((1 if u >= 0 else 0, 1 if v >= 0 else 0))

    aniso = compute_gradient_anisotropy(lrm, res_m=res_m)

    return {
        "code": code,
        "lat": meta.get("lat"),
        "lon": meta.get("lon"),
        "strip_id": meta.get("strip_id"),
        "dist_from_card08_km": meta.get("dist_from_card08_km", 0.0),
        "nearest_any_known": meta.get("nearest_any_known"),
        "n_mounds_total_0p33m": len(mounds_detail),
        "n_closed_360_mounds": len(closed_m),
        "n_ortho_rect_mounds": len(ortho_m),
        "closed_quadrants_of_4": len(quads_closed),
        "dom_azimuth_deg": round(dom_az_mod90, 1),
        "ortho_coherence_0p33m": round(ortho_coh, 3),
        "strike_anisotropy_ratio": round(aniso, 2),
        "plaza_area_lt2deg_m2": round(plaza_lt2_m2, 1),
        "plaza_area_lt3deg_m2": round(plaza_lt3_m2, 1),
        "plaza_area_lt4p5deg_m2": round(plaza_lt4p5_m2, 1),
        "plaza_median_slope_deg": round(plaza_med_slope, 2),
        "scene_median_slope_deg": round(bg_med_slope, 2),
        "laz_pts_window": int(len(px)),
        "laz_pts_per_m2": meta.get("laz_pts_per_m2"),
        "mound_pts_per_m2": meta.get("mound_pts_per_m2"),
        "plaza_pts_per_m2": round(plaza_pts_m2, 3),
        "pct_mound_void_gt6m": meta.get("pct_mound_void_gt6m", 0.0),
        "top_mounds": mounds_detail[:8],
    }


def extract_swath_profile(
    dem: np.ndarray,
    px: np.ndarray,
    py: np.ndarray,
    pz: np.ndarray,
    pcls: np.ndarray,
    hm: float,
    p0: tuple[float, float],
    p1: tuple[float, float],
    half_swath_m: float = 7.0,
) -> dict:
    H, W = dem.shape
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    L = math.hypot(dx, dy)
    ux, uy = dx / L, dy / L

    s_pt = (px - p0[0]) * ux + (py - p0[1]) * uy
    d_pt = np.abs(-(px - p0[0]) * uy + (py - p0[1]) * ux)
    m = (s_pt >= 0) & (s_pt <= L) & (d_pt <= half_swath_m)

    s_grid = np.linspace(0, L, 360)
    gx = p0[0] + s_grid * ux
    gy = p0[1] + s_grid * uy
    cols = np.clip(((gx + hm) / (2 * hm) * W).astype(int), 0, W - 1)
    rows = np.clip(((hm - gy) / (2 * hm) * H).astype(int), 0, H - 1)
    dem_prof = dem[rows, cols]

    return {
        "L": L,
        "p0": p0,
        "p1": p1,
        "s_grid": s_grid,
        "dem_prof": dem_prof,
        "pt_s": s_pt[m],
        "pt_z": pz[m],
        "pt_cls": pcls[m],
    }


def main() -> None:
    # 1. Audit #08 Positive Control (270x270m crop from card08_bundle_600m.npz)
    d08 = np.load(WORK_DIR / "card08_bundle_600m.npz")
    cx08, cy08, res08 = float(d08["cx"]), float(d08["cy"]), float(d08["res_m"])
    hm = 135.0
    r0 = int(round((600.0 - hm) / res08))
    r1 = int(round((600.0 + hm) / res08))
    dem08 = d08["dem33"][r0:r1, r0:r1]
    lrm08 = d08["lrm33"][r0:r1, r0:r1]
    g1_08 = d08["g1_33"][r0:r1, r0:r1]
    m08 = (np.abs(d08["ot_x"] - cx08) <= hm) & (np.abs(d08["ot_y"] - cy08) <= hm)
    px08 = d08["ot_x"][m08] - cx08
    py08 = d08["ot_y"][m08] - cy08
    pz08 = d08["ot_z"][m08]
    pc08 = d08["ot_cls"][m08]

    res_all = []
    r_ctrl = audit_single_complex(
        "#08 (Nuevo Jalisco)",
        dem08, lrm08, g1_08, px08, py08, pz08, pc08, res08, hm,
        {
            "lat": 16.719521, "lon": -91.168757, "strip_id": "ChiapsCampeche_NFI_l0s473",
            "dist_from_card08_km": 0.0,
            "nearest_any_known": {"name": "Nuevo Jalisco (Carrasco 1981; Tovalin 2004)", "dist_km": 0.02},
            "laz_pts_per_m2": 0.471, "mound_pts_per_m2": 0.562, "pct_mound_void_gt6m": 0.0,
        },
    )
    res_all.append(r_ctrl)

    # 2. Audit LC-01 .. LC-14
    curated = json.loads((WORK_DIR / "lacandon_top14_curated.json").read_text())
    npz_map = {}
    for item in curated:
        npz = np.load(item["cache_npz"])
        npz_map[item["lc_code"]] = (item, npz)
        rc = audit_single_complex(
            item["lc_code"],
            npz["dem33"], npz["lrm33"], npz["g1_33"],
            npz["px"], npz["py"], npz["pz"], npz["pcls"],
            float(npz["res_m"]), 135.0, item,
        )
        res_all.append(rc)

    DATA_OUT.parent.mkdir(parents=True, exist_ok=True)
    DATA_OUT.write_text(json.dumps(res_all, indent=2))

    print(f"{'ID':18s} | {'TotM':>4s} {'Cls360':>6s} {'OrthoR':>6s} {'Q4':>2s} | {'Aniso':>5s} | {'Plz<2°':>7s} {'Plz<3°':>7s} {'PlzSlp':>6s} {'BgSlp':>5s} | {'LAZ/m2':>6s} {'Mnd/m2':>6s} {'Plz/m2':>6s}")
    print("-" * 118)
    for r in res_all:
        print(
            f"{r['code']:18s} | {r['n_mounds_total_0p33m']:4d} {r['n_closed_360_mounds']:6d} {r['n_ortho_rect_mounds']:6d} {r['closed_quadrants_of_4']:2d} | "
            f"{r['strike_anisotropy_ratio']:5.2f} | {r['plaza_area_lt2deg_m2']:6.0f}m² {r['plaza_area_lt3deg_m2']:6.0f}m² {r['plaza_median_slope_deg']:5.2f}° {r['scene_median_slope_deg']:4.1f}° | "
            f"{r['laz_pts_per_m2']:6.2f} {r['mound_pts_per_m2']:6.2f} {r['plaza_pts_per_m2']:6.2f}"
        )

    # Print detailed mound list for LC-01, LC-02, LC-03, LC-04, LC-08, LC-09, LC-10, LC-14
    for target in ["#08 (Nuevo Jalisco)", "LC-01", "LC-02", "LC-03", "LC-04", "LC-06", "LC-08", "LC-09", "LC-10", "LC-13", "LC-14"]:
        rt = next(x for x in res_all if x["code"] == target)
        print(f"\n=== {target} ({rt['lat']}, {rt['lon']}) Top Mounds ===")
        for i, m in enumerate(rt["top_mounds"][:5], 1):
            print(
                f"  M{i}: pos=({m['cx_m']:+5.1f},{m['cy_m']:+5.1f}) sz={m['len_maj_m']:.1f}x{m['len_min_m']:.1f}m "
                f"(asp={m['aspect']:.2f}, rect={m['rect_fill']:.2f}) LRM=+{m['max_lrm_m']:.2f}m "
                f"DEMrays={m['closed_rays_dem_of_8']}/8 (minDrop={m['min_dem_drop_m']:.2f}m) "
                f"topSlp={m['summit_med_slope_deg']:.1f}° LAZ={m['laz_pts_on']}p ({m['laz_pts_m2']:.2f}/m²) "
                f"closed360={m['is_closed_360']} ortho={m['is_ortho_rect']}"
            )

    # Render diagnostic forensic cross-section sheet for LC-01 + all 14 complexes so we inspect visually!
    fig, axes = plt.subplots(5, 6, figsize=(24, 19), dpi=135, facecolor="#0b0f17")
    ext = [-135, 135, -135, 135]
    for idx, r in enumerate(res_all):
        row = idx // 3
        col0 = (idx % 3) * 2
        ax_map = axes[row, col0]
        ax_prf = axes[row, col0 + 1]

        if r["code"].startswith("#08"):
            dem, lrm, g1, px, py, pz, pcls, res_m = dem08, lrm08, g1_08, px08, py08, pz08, pc08, res08
        else:
            _, npz = npz_map[r["code"]]
            dem, lrm, g1 = npz["dem33"], npz["lrm33"], npz["g1_33"]
            px, py, pz, pcls = npz["px"], npz["py"], npz["pz"], npz["pcls"]
            res_m = float(npz["res_m"])

        ax_map.imshow(g1, extent=ext, origin="upper")
        # Highlight closed 360° mounds in green circles, unclosed ridges in red/orange!
        for m in r["top_mounds"][:8]:
            col_c = "#4ade80" if m["is_closed_360"] else "#fb7185"
            ax_map.scatter([m["pk_x_m"]], [m["pk_y_m"]], s=38, facecolors="none", edgecolors=col_c, linewidths=1.5)

        # Define longitudinal & transverse transects through top 2 mounds or principal axis
        if len(r["top_mounds"]) >= 2:
            m1, m2 = r["top_mounds"][0], r["top_mounds"][1]
            dx_m, dy_m = m2["pk_x_m"] - m1["pk_x_m"], m2["pk_y_m"] - m1["pk_y_m"]
            dist_12 = math.hypot(dx_m, dy_m)
            if dist_12 > 12.0:
                ux, uy = dx_m / dist_12, dy_m / dist_12
            else:
                ang = math.radians(r["dom_azimuth_deg"])
                ux, uy = math.cos(ang), math.sin(ang)
            midx, midy = 0.5 * (m1["pk_x_m"] + m2["pk_x_m"]), 0.5 * (m1["pk_y_m"] + m2["pk_y_m"])
            pL0 = (max(-125.0, min(125.0, midx - 85.0 * ux)), max(-125.0, min(125.0, midy - 85.0 * uy)))
            pL1 = (max(-125.0, min(125.0, midx + 85.0 * ux)), max(-125.0, min(125.0, midy + 85.0 * uy)))
            # Transverse across m1
            pT0 = (max(-125.0, min(125.0, m1["pk_x_m"] - 65.0 * (-uy))), max(-125.0, min(125.0, m1["pk_y_m"] - 65.0 * ux)))
            pT1 = (max(-125.0, min(125.0, m1["pk_x_m"] + 65.0 * (-uy))), max(-125.0, min(125.0, m1["pk_y_m"] + 65.0 * ux)))
        else:
            pL0, pL1 = (-80.0, 0.0), (80.0, 0.0)
            pT0, pT1 = (0.0, -80.0), (0.0, 80.0)

        ax_map.plot([pL0[0], pL1[0]], [pL0[1], pL1[1]], color="#00e5ff", lw=1.4, ls="--")
        ax_map.plot([pT0[0], pT1[0]], [pT0[1], pT1[1]], color="#fde047", lw=1.4, ls=":")
        ax_map.set_xlim(-135, 135)
        ax_map.set_ylim(-135, 135)
        ax_map.set_xticks([])
        ax_map.set_yticks([])
        ax_map.set_title(
            f"{r['code']} | Cls360:{r['n_closed_360_mounds']}/{r['n_mounds_total_0p33m']} (Ortho:{r['n_ortho_rect_mounds']}) {r['closed_quadrants_of_4']}/4Q\n"
            f"Aniso:{r['strike_anisotropy_ratio']:.2f} Plz<2°:{r['plaza_area_lt2deg_m2']:.0f}m² (slp:{r['plaza_median_slope_deg']:.1f}°)",
            color="white", fontsize=7.6, fontweight="bold",
        )

        sL = extract_swath_profile(dem, px, py, pz, pcls, 135.0, pL0, pL1, half_swath_m=6.0)
        sT = extract_swath_profile(dem, px, py, pz, pcls, 135.0, pT0, pT1, half_swath_m=6.0)
        ax_prf.set_facecolor("#111827")
        vpL = np.isfinite(sL["dem_prof"])
        vpT = np.isfinite(sT["dem_prof"])
        if vpL.any():
            z0L = float(np.nanmedian(sL["dem_prof"][vpL]))
            ax_prf.plot(sL["s_grid"][vpL], sL["dem_prof"][vpL] - z0L, color="#00e5ff", lw=1.6, label="Longitudinal DEM")
            ax_prf.scatter(sL["pt_s"], sL["pt_z"] - z0L, c="#00e5ff", s=3.5, alpha=0.55)
        if vpT.any():
            z0T = float(np.nanmedian(sT["dem_prof"][vpT]))
            ax_prf.plot(sT["s_grid"][vpT], sT["dem_prof"][vpT] - z0T, color="#fde047", lw=1.5, ls="--", label="Transverse DEM")
            ax_prf.scatter(sT["pt_s"], sT["pt_z"] - z0T, c="#fde047", s=3.5, alpha=0.55)
        ax_prf.axhline(0.0, color="#64748b", lw=0.7, ls=":")
        ax_prf.tick_params(colors="#cbd5e1", labelsize=7)
        ax_prf.set_title(f"Profiles (ΔZ m) | LAZ:{r['laz_pts_per_m2']:.2f}/m² Mnd:{r['mound_pts_per_m2']:.2f}/m²", color="#cbd5e1", fontsize=7.5)
        ax_prf.grid(True, alpha=0.2)

    fig.tight_layout(pad=0.6)
    diag_out = WORK_DIR / "diag_all14_forensic_karst_test.png"
    fig.savefig(diag_out, facecolor="#0b0f17")
    plt.close(fig)
    print("\nSaved diagnostic forensic sheet:", diag_out)


if __name__ == "__main__":
    main()
