#!/usr/bin/env python3
"""Exhaustive 2-Stage Forensic Scanner of the Surrounding Lacandon Jungle Corridor:
  - Full length of ChiapsCampeche_NFI_l0s473 (containing Candidate #08 / Nuevo Jalisco)
  - Contiguous NW-SE Lacandon strips ChiapsCampeche_NFI_l0s472 .. l0s480
  - Nearby Upper Usumacinta / Montes Azules Chiaps1_GLAS & Chiaps1_NFI strips (<= 38 km)

Enforces strict raw .LAZ ground-point density verification + 0.33m G1/DEM + 1.8km Esri
satellite canopy check so zero TIN/kriging canopy-void triangles can pass!
"""
from __future__ import annotations

import json
import math
import os
import re
import sys
import urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import laspy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.transform import xy as rio_xy
from rasterio.warp import transform as rio_transform
from rasterio.windows import from_bounds
from scipy import ndimage as ndi

sys.path.insert(0, "/agent/maya-lidar-discovery/scripts")
sys.path.insert(0, "/agent/maya-lidar-discovery/scripts/open_lidar_hunt_v2")
import raw_pointcloud_tools as rpt
import v2_scene_tools as vst

OUT_DIR = Path("/agent/lidar/work/lacandon_deep_dive")
OUT_DIR.mkdir(parents=True, exist_ok=True)

LAT0, LON0 = 16.719521, -91.168757


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(min(1.0, a)))


def build_literature_and_registry_catalog() -> list[dict]:
    """Combine EAAMS (5,223), INAH SIC, OSM, AND published INAH Chiapas field-report
    UTM coordinates (Tovalin et al. 2004, 2020; Golden et al. 2016) converted from NAD27/WGS84!
    """
    sites = []
    eaams_gj = json.loads(Path("/agent/lidar/work/ce_deep/eaams/maya_sites.geojson").read_text())
    for feat in eaams_gj["features"]:
        lon, lat = feat["geometry"]["coordinates"][:2]
        if 15.5 <= lat <= 18.0 and -92.5 <= lon <= -90.0:
            p = feat.get("properties") or {}
            sites.append({"name": p.get("name") or "EAAMS", "src": "EAAMS", "lat": float(lat), "lon": float(lon)})

    for p_json, tag in [
        ("/agent/lidar/work/gliht_west_mexico/registries/inah_sic_all_mexico.json", "INAH_SIC"),
        ("/agent/lidar/work/gliht_west_mexico/registries/osm_mexico_archaeology.json", "OSM"),
    ]:
        pth = Path(p_json)
        if pth.exists():
            for s in json.loads(pth.read_text()):
                lat, lon = s.get("lat"), s.get("lon")
                if lat is not None and lon is not None and 15.5 <= float(lat) <= 18.0 and -92.5 <= float(lon) <= -90.0:
                    sites.append({"name": s.get("name") or tag, "src": tag, "lat": float(lat), "lon": float(lon)})

    # Published INAH field-survey UTMs in Tovalin, Ortiz & Velazquez de Leon (2004, Bolom 1:55-90)
    # Note: INEGI 1:50,000 E15D56/E15D66 maps use NAD27 UTM Zone 15N (EPSG:26715)
    tr_nad27 = Transformer.from_crs("EPSG:26715", "EPSG:4326", always_xy=True)
    lit_nad27 = [
        ("Nuevo Jalisco (Carrasco 1981; Tovalin 2004)", 695273.0, 1849232.0),
        ("Plan de Ayutla (Tovalin 2004; Martos 2009)", 683240.0, 1856391.0),
        ("Nuevo Chetumal - La Lucha (Healey 1950; Tovalin 2004)", 678673.0, 1865702.0),
        ("Rancho Ojo de Agua (Blom 1948; Tovalin 2004)", 696602.0, 1864441.0),
    ]
    for nm, ue, un in lit_nad27:
        lon, lat = tr_nad27.transform(ue, un)
        sites.append({"name": nm, "src": "INAH_Lit_NAD27", "lat": float(lat), "lon": float(lon)})

    return sites


def get_chm_path_for_strip(strip_id: str) -> str | None:
    try:
        camp, ls = rpt.nasa_campaign_and_ls(strip_id)
        p = Path(f"/agent/lidar/work/carbon_raw_hunt/chm/extracted/{camp}/{camp}_{ls}_CHM.tif")
        return str(p) if p.exists() else None
    except Exception:
        return None


def norm_lrm_1m(dem: np.ndarray, valid: np.ndarray, sigma_m: float = 22.0) -> np.ndarray:
    h, w = dem.shape
    step = 3
    sv = valid[::step, ::step].astype(np.float32)
    sd = np.where(valid[::step, ::step], dem[::step, ::step], 0.0).astype(np.float32)
    sig = sigma_m / step
    num = ndi.gaussian_filter(sd, sigma=sig, mode="constant", cval=0.0)
    den = ndi.gaussian_filter(sv, sigma=sig, mode="constant", cval=0.0)
    tr_s = np.where(den > 1e-3, num / np.maximum(den, 1e-3), np.nanmedian(dem[valid])).astype(np.float32)
    tr = ndi.zoom(tr_s, (h / tr_s.shape[0], w / tr_s.shape[1]), order=1)[:h, :w]
    lrm = dem - tr
    lrm[~valid] = np.nan
    return lrm


def mound_shape_metrics(comp: np.ndarray, lrm_w: np.ndarray) -> dict:
    ys, xs = np.where(comp)
    area = float(len(xs))
    pts = np.column_stack([xs.astype(np.float64), ys.astype(np.float64)])
    mean = pts.mean(axis=0)
    X = pts - mean
    cov = X.T @ X / max(len(X) - 1, 1)
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]
    evecs = evecs[:, order]
    rot = X @ evecs
    w_len = float(max(rot[:, 0].max() - rot[:, 0].min() + 1.0, 1.0))
    h_len = float(max(rot[:, 1].max() - rot[:, 1].min() + 1.0, 1.0))
    aspect = max(w_len / h_len, h_len / w_len)
    rect_fill = min(1.0, area / (w_len * h_len))
    az_deg = float(math.degrees(math.atan2(evecs[1, 0], evecs[0, 0])) % 90.0)
    vals = lrm_w[ys, xs]
    h_p95 = float(np.nanpercentile(vals, 95))
    h_max = float(np.nanmax(vals))
    return {
        "area_m2": area,
        "w_m": w_len,
        "h_m": h_len,
        "aspect": aspect,
        "rect": rect_fill,
        "az_mod90": az_deg,
        "h_p95": h_p95,
        "h_max": h_max,
        "cx_px": float(xs.mean()),
        "cy_px": float(ys.mean()),
    }


def scan_single_strip_1m(args: tuple[str, str, str | None]) -> list[dict]:
    strip_id, dem_path, chm_path = args
    with rasterio.open(dem_path) as ds:
        dem = ds.read(1).astype(np.float32)
        if ds.nodata is not None:
            dem = np.where(dem == ds.nodata, np.nan, dem)
        dem = np.where((np.abs(dem) > 1e6) | (dem < -100), np.nan, dem)
        valid = np.isfinite(dem)
        if valid.sum() < 5000:
            return []
        transform = ds.transform
        crs = ds.crs

    chm_aligned = None
    if chm_path and Path(chm_path).exists():
        with rasterio.open(chm_path) as cds:
            # Read CHM over DEM bounds
            b = rasterio.transform.array_bounds(dem.shape[0], dem.shape[1], transform)
            cw = from_bounds(b[0], b[1], b[2], b[3], cds.transform)
            chm_aligned = cds.read(
                1, window=cw, out_shape=dem.shape, boundless=True, fill_value=np.nan
            ).astype(np.float32)
            if cds.nodata is not None:
                chm_aligned = np.where(chm_aligned == cds.nodata, np.nan, chm_aligned)
            chm_aligned = np.where((chm_aligned < -2) | (chm_aligned > 65), np.nan, chm_aligned)

    lrm = norm_lrm_1m(dem, valid, sigma_m=22.0)
    fill_dem = np.where(valid, dem, float(np.nanmedian(dem[valid])))
    gy, gx = np.gradient(ndi.gaussian_filter(fill_dem, sigma=1.5))
    slope_deg = np.degrees(np.arctan(np.hypot(gx, gy)))

    H, W = dem.shape
    win_half = 115  # 230x230m window
    step = 55 if strip_id == "ChiapsCampeche_NFI_l0s473" else 70
    cands = []

    for r in range(win_half, H - win_half, step):
        for c in range(win_half, W - win_half, step):
            v_w = valid[r - win_half : r + win_half, c - win_half : c + win_half]
            if v_w.mean() < 0.62:
                continue
            lrm_w = lrm[r - win_half : r + win_half, c - win_half : c + win_half]
            slp_w = slope_deg[r - win_half : r + win_half, c - win_half : c + win_half]

            # Detect discrete mounds / platforms (LRM >= 1.05m)
            m_mask = v_w & (lrm_w >= 1.05)
            m_mask = ndi.binary_opening(m_mask, iterations=1)
            labeled, nlab = ndi.label(m_mask)
            if nlab < 2:
                continue

            mounds = []
            for lab in range(1, nlab + 1):
                comp = labeled == lab
                cnt = int(comp.sum())
                if cnt < 32 or cnt > 9500:
                    continue
                mm = mound_shape_metrics(comp, lrm_w)
                # Reject extremely elongated >6.5:1 knife-edge natural ridges unless compact platform
                if mm["aspect"] > 6.2 or mm["h_max"] < 1.35:
                    continue
                mounds.append(mm)

            if len(mounds) < 2:
                continue

            mounds.sort(key=lambda x: x["h_p95"] * math.sqrt(x["area_m2"]), reverse=True)
            top_m = mounds[:10]
            rect_m = [m for m in top_m if m["rect"] >= 0.56 and m["aspect"] <= 4.2]
            if len(rect_m) < 1:
                continue

            # Check leveled courtyard / plaza area inside inner 150x150m
            in_sl = slice(win_half - 75, win_half + 75)
            plaza_mask = (
                v_w[in_sl, in_sl]
                & (np.abs(lrm_w[in_sl, in_sl]) <= 0.85)
                & (slp_w[in_sl, in_sl] <= 4.5)
            )
            plaza_area = float(plaza_mask.sum())

            # Quadrant surround around centroid of inner plaza/mounds
            mx0 = np.mean([m["cx_px"] for m in top_m[: min(5, len(top_m))]])
            my0 = np.mean([m["cy_px"] for m in top_m[: min(5, len(top_m))]])
            quads = set()
            for m in top_m:
                qx = 1 if m["cx_px"] >= mx0 else 0
                qy = 1 if m["cy_px"] >= my0 else 0
                quads.add((qx, qy))

            # Azimuth coherence (mod 90 deg)
            angs = np.radians([m["az_mod90"] * 4.0 for m in top_m])
            ortho_coh = float(math.hypot(np.mean(np.cos(angs)), np.mean(np.sin(angs))))

            heights = sorted([m["h_max"] for m in top_m], reverse=True)
            h1 = heights[0]
            h2 = heights[1] if len(heights) > 1 else 0.0
            mean_rect = float(np.mean([m["rect"] for m in top_m[:5]]))

            # CHM stats
            chm_med, chm_p90, cleared_frac = 20.0, 28.0, 0.0
            if chm_aligned is not None:
                cw_arr = chm_aligned[r - win_half : r + win_half, c - win_half : c + win_half]
                cv = np.isfinite(cw_arr) & v_w
                if cv.sum() > 200:
                    chm_med = float(np.median(cw_arr[cv]))
                    chm_p90 = float(np.percentile(cw_arr[cv], 90))
                    cleared_frac = float((cw_arr[cv] < 3.0).mean())

            # Require intact forest canopy (median CHM >= 10m and cleared_frac <= 0.18)
            if chm_med < 10.0 or cleared_frac > 0.18:
                continue

            # Architectural score favoring genuine multi-mound rectangular plaza groups!
            score = (
                4.5 * min(len(top_m), 9)
                + 6.5 * min(len(rect_m), 6)
                + 5.5 * len(quads)
                + 14.0 * mean_rect
                + 8.0 * ortho_coh
                + 2.8 * min(h1, 9.0)
                + 2.2 * min(h2, 6.5)
                + 0.0045 * min(plaza_area, 3200.0)
            )
            if score < 46.0:
                continue

            ux, uy = rio_xy(transform, r, c)
            lon_l, lat_l = rio_transform(crs, "EPSG:4326", [ux], [uy])
            cands.append(
                {
                    "strip_id": strip_id,
                    "dem_path": dem_path,
                    "chm_path": chm_path,
                    "lat": round(float(lat_l[0]), 6),
                    "lon": round(float(lon_l[0]), 6),
                    "utm_e": round(float(ux), 1),
                    "utm_n": round(float(uy), 1),
                    "score_1m": round(score, 2),
                    "n_mounds": len(top_m),
                    "n_rect": len(rect_m),
                    "n_quads": len(quads),
                    "mean_rect": round(mean_rect, 3),
                    "ortho_coh": round(ortho_coh, 3),
                    "h1_m": round(h1, 2),
                    "h2_m": round(h2, 2),
                    "plaza_m2": round(plaza_area, 1),
                    "chm_med_m": round(chm_med, 1),
                    "chm_p90_m": round(chm_p90, 1),
                    "cleared_frac": round(cleared_frac, 4),
                    "dist_from_card08_km": round(haversine_km(LAT0, LON0, float(lat_l[0]), float(lon_l[0])), 2),
                }
            )
    return cands


def nms_cands(cands: list[dict], radius_m: float = 280.0) -> list[dict]:
    cands = sorted(cands, key=lambda x: x["score_1m"], reverse=True)
    kept = []
    for c in cands:
        if any(haversine_km(c["lat"], c["lon"], k["lat"], k["lon"]) * 1000.0 < radius_m for k in kept):
            continue
        kept.append(c)
    return kept


def ensure_laz_only(strip_id: str) -> Path:
    laz_key = rpt.LAZ_MAP[strip_id]
    laz_path = rpt.RAW_CACHE / Path(laz_key).name
    if not laz_path.exists():
        url = f"https://opentopography.s3.sdsc.edu/pc-bulk/{laz_key}"
        urllib.request.urlretrieve(url, laz_path)
    return laz_path


def verify_candidate_0p33m_and_laz(cand: dict, laz_pts_by_strip: dict, catalog: list[dict]) -> dict:
    sid = cand["strip_id"]
    lat, lon = cand["lat"], cand["lon"]
    hm = 135.0  # 270x270m standard card window

    u_dem = f"/vsicurl/https://opentopography.s3.sdsc.edu/raster/{rpt.DEM33_MAP[sid]}"
    u_g1 = f"/vsicurl/https://opentopography.s3.sdsc.edu/raster/{rpt.G1_33_MAP[sid]}"

    with rasterio.open(u_dem) as ds:
        xs, ys = rio_transform("EPSG:4326", ds.crs, [lon], [lat])
        cx, cy = float(xs[0]), float(ys[0])
        win = from_bounds(cx - hm, cy - hm, cx + hm, cy + hm, ds.transform)
        dem33 = ds.read(1, window=win, boundless=True, fill_value=np.nan).astype(np.float32)
        if ds.nodata is not None:
            dem33 = np.where(dem33 == ds.nodata, np.nan, dem33)
        dem33 = np.where((np.abs(dem33) > 1e6) | (dem33 < -100), np.nan, dem33)
        res_m = float(abs(ds.res[0]))

    with rasterio.open(u_g1) as gds:
        wg = from_bounds(cx - hm, cy - hm, cx + hm, cy + hm, gds.transform)
        g1_33 = np.moveaxis(
            gds.read([1, 2, 3], window=wg, out_shape=(3, dem33.shape[0], dem33.shape[1]), boundless=True, fill_value=18),
            0,
            -1,
        )

    lrm33 = rpt.compute_lrm(dem33, res_m=res_m, sigma_m=20.0)
    valid = np.isfinite(dem33)
    valid_area_m2 = float(valid.sum() * res_m * res_m)

    # Extract raw .laz ground points inside 270x270m
    lx, ly, lz, lcls = laz_pts_by_strip[sid]
    m_in = (np.abs(lx - cx) <= hm) & (np.abs(ly - cy) <= hm)
    px, py, pz, pcls = lx[m_in], ly[m_in], lz[m_in], lcls[m_in]
    n_laz = int(len(px))
    pts_per_m2 = n_laz / max(valid_area_m2, 1.0)

    # Check raw .laz points STRICTLY on elevated mounds (LRM33 >= 1.2m in central 220x220m)!
    H, W = dem33.shape
    p_col = np.clip(((px - (cx - hm)) / (2 * hm) * W).astype(int), 0, W - 1)
    p_row = np.clip((((cy + hm) - py) / (2 * hm) * H).astype(int), 0, H - 1)
    pt_grid = np.zeros((H, W), dtype=np.int32)
    np.add.at(pt_grid, (p_row, p_col), 1)
    dist_pt_m = ndi.distance_transform_edt(pt_grid == 0) * res_m

    yy, xx = np.mgrid[0:H, 0:W]
    in_core = (np.abs((xx / W - 0.5) * 2 * hm) <= 105.0) & (np.abs((yy / H - 0.5) * 2 * hm) <= 105.0)
    mound_mask = valid & in_core & (lrm33 >= 1.2)
    mound_area_m2 = float(mound_mask.sum() * res_m * res_m)
    mound_pts = int(pt_grid[mound_mask].sum())
    mound_pts_m2 = mound_pts / max(mound_area_m2, 1.0)
    pct_mound_void_gt6m = float((dist_pt_m[mound_mask] > 6.0).mean() * 100.0) if mound_mask.any() else 100.0

    # Fetch 1.8 x 1.8 km Esri World Imagery optical satellite
    dlat = 900.0 / 110574.0
    dlon = 900.0 / (111320.0 * math.cos(math.radians(lat)))
    sat_1800 = vst.fetch_esri_bbox_rgb(lon - dlon, lat - dlat, lon + dlon, lat + dlat, out_w=320, out_h=320, zoom_level=15)
    # Estimate green closed-forest fraction on satellite RGB (reject bright tan roads/pastures)
    r_ch, g_ch, b_ch = sat_1800[..., 0].astype(float), sat_1800[..., 1].astype(float), sat_1800[..., 2].astype(float)
    forest_px = (g_ch >= r_ch * 0.94) & (g_ch >= b_ch * 1.02) & ((r_ch + g_ch + b_ch) < 390)
    sat_forest_pct = round(100.0 * float(forest_px.mean()), 1)

    # Nearest known site across ALL registries + Tovalin/Golden literature!
    best_by_src = {}
    nearest_any = None
    for s in catalog:
        dkm = haversine_km(lat, lon, s["lat"], s["lon"])
        src = s["src"]
        if src not in best_by_src or dkm < best_by_src[src]["dist_km"]:
            best_by_src[src] = {"name": s["name"], "src": src, "dist_km": round(dkm, 2)}
        if nearest_any is None or dkm < nearest_any["dist_km"]:
            nearest_any = {"name": s["name"], "src": src, "dist_km": round(dkm, 2)}

    # Save cached arrays for instant contact-sheet / showcase rendering!
    slug = f"{sid}_{lat:.5f}_{lon:.5f}"
    np.savez_compressed(
        OUT_DIR / f"cand_{slug}.npz",
        dem33=dem33,
        g1_33=g1_33,
        lrm33=lrm33,
        sat_1800=sat_1800,
        px=px - cx,
        py=py - cy,
        pz=pz,
        pcls=pcls,
        cx=cx,
        cy=cy,
        res_m=res_m,
    )

    is_solid_ground = (pts_per_m2 >= 0.25) and (mound_pts >= 80) and (mound_pts_m2 >= 0.22) and (pct_mound_void_gt6m <= 22.0)

    return {
        **cand,
        "valid_area_m2": round(valid_area_m2, 1),
        "laz_pts_270m": n_laz,
        "laz_pts_per_m2": round(pts_per_m2, 3),
        "mound_area_gt1p2m_m2": round(mound_area_m2, 1),
        "mound_laz_pts": mound_pts,
        "mound_pts_per_m2": round(mound_pts_m2, 3),
        "pct_mound_void_gt6m": round(pct_mound_void_gt6m, 1),
        "sat_1800m_forest_pct": sat_forest_pct,
        "passes_laz_ground_gate": bool(is_solid_ground),
        "nearest_any_known": nearest_any,
        "nearest_by_source": best_by_src,
        "cache_npz": str(OUT_DIR / f"cand_{slug}.npz"),
    }


def main() -> None:
    inv = json.loads((OUT_DIR / "lacandon_inventory.json").read_text())
    catalog = build_literature_and_registry_catalog()

    # Select target Lacandon corridor strips:
    # All ChiapsCampeche_NFI_l0s472..l0s481 + all Chiaps1 strips within 38 km of #08!
    target_strips = []
    for s in inv["nearby_strips"]:
        sid = s["strip_id"]
        if (
            re.match(r"ChiapsCampeche_NFI_l0s4(7[2-9]|8[01])$", sid)
            or (s["approx_dist_km"] <= 38.0 and sid.startswith("Chiaps"))
        ):
            target_strips.append((sid, s["path"], get_chm_path_for_strip(sid)))

    print(f"Scanning {len(target_strips)} Lacandon rainforest strips at 1m + CHM...", flush=True)
    all_raw = []
    with ProcessPoolExecutor(max_workers=min(12, os.cpu_count() or 4)) as ex:
        for res in ex.map(scan_single_strip_1m, target_strips):
            all_raw.extend(res)

    print(f"Raw window detections before NMS: {len(all_raw)}", flush=True)

    # Separate l0s473 detections (excluding #08 itself within 220m) vs surrounding strips
    cands_473 = [c for c in all_raw if c["strip_id"] == "ChiapsCampeche_NFI_l0s473" and c["dist_from_card08_km"] >= 0.23]
    cands_473 = nms_cands(cands_473, radius_m=240.0)[:10]

    cands_other = [c for c in all_raw if c["strip_id"] != "ChiapsCampeche_NFI_l0s473"]
    cands_other = nms_cands(cands_other, radius_m=320.0)

    # Diversify slightly across strips while keeping highest architectural scores (max 6 per strip)
    per_strip = {}
    sel_other = []
    for c in cands_other:
        sid = c["strip_id"]
        if per_strip.get(sid, 0) >= 6:
            continue
        per_strip[sid] = per_strip.get(sid, 0) + 1
        sel_other.append(c)
        if len(sel_other) >= 26:
            break

    shortlist = cands_473 + sel_other
    print(f"Shortlisted {len(cands_473)} hits on l0s473 + {len(sel_other)} hits on neighboring Lacandon strips = {len(shortlist)} total for 0.33m + .LAZ audit!", flush=True)

    # Pre-load .laz point clouds once per unique strip
    needed_strips = sorted({c["strip_id"] for c in shortlist})
    print(f"Ensuring .laz ground clouds for {len(needed_strips)} strips: {needed_strips}", flush=True)
    laz_pts_by_strip = {}
    for sid in needed_strips:
        lp = ensure_laz_only(sid)
        las = laspy.read(lp)
        laz_pts_by_strip[sid] = (
            np.array(las.x, dtype=np.float64),
            np.array(las.y, dtype=np.float64),
            np.array(las.z, dtype=np.float32),
            np.array(las.classification, dtype=np.uint8),
        )

    print("Running 0.33m DEM + 0.33m G1 + raw .LAZ ground-density + 1.8km Esri audit...", flush=True)
    verified = []
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(verify_candidate_0p33m_and_laz, c, laz_pts_by_strip, catalog) for c in shortlist]
        for f in futs:
            verified.append(f.result())

    verified.sort(key=lambda x: (x["passes_laz_ground_gate"], x["score_1m"]), reverse=True)
    (OUT_DIR / "lacandon_corridor_verified.json").write_text(json.dumps(verified, indent=2))

    # Render a quick diagnostic visual sheet of ALL verified candidates (0.33m G1 + raw .laz points + LRM + Sat)
    # so we can visually inspect every single candidate and pick only 100% genuine Maya architecture!
    passing = [v for v in verified if v["passes_laz_ground_gate"]]
    failing = [v for v in verified if not v["passes_laz_ground_gate"]]
    print(f"\n=== AUDIT SUMMARY: {len(passing)} PASSED raw .LAZ ground gate | {len(failing)} REJECTED as sparse/void ===", flush=True)
    for i, v in enumerate(verified, 1):
        status = "PASS" if v["passes_laz_ground_gate"] else "VOID-FAIL"
        nk = v["nearest_any_known"]
        print(
            f"{i:2d}. [{status:9s}] {v['strip_id']:26s} ({v['lat']:.5f}, {v['lon']:.5f}) | "
            f"d#08={v['dist_from_card08_km']:5.2f}km | score={v['score_1m']:5.1f} | "
            f"mounds={v['n_mounds']}(rect:{v['n_rect']}) H={v['h1_m']:.1f}m Plz={v['plaza_m2']:.0f}m2 | "
            f"LAZ={v['laz_pts_270m']:5d} ({v['laz_pts_per_m2']:.2f}/m2, mnd={v['mound_laz_pts']}p={v['mound_pts_per_m2']:.2f}/m2, void>6m={v['pct_mound_void_gt6m']:.0f}%) | "
            f"SatForest={v['sat_1800m_forest_pct']:.0f}% CHM={v['chm_med_m']:.0f}m | "
            f"Near={nk['name']} ({nk['dist_km']}km, {nk['src']})"
        )

    # Render visual screening grid of top 24 passing candidates
    top_vis = passing[:24]
    ncols = 3
    nrows = math.ceil(len(top_vis) / ncols)
    fig, axes = plt.subplots(nrows, ncols * 3, figsize=(24, 2.9 * nrows), dpi=130, facecolor="#0b0f17")
    if nrows == 1:
        axes = axes[None, :]
    ext = [-135, 135, -135, 135]
    for idx, v in enumerate(top_vis):
        r_i, c_i = divmod(idx, ncols)
        ax_g1 = axes[r_i, c_i * 3]
        ax_pt = axes[r_i, c_i * 3 + 1]
        ax_sat = axes[r_i, c_i * 3 + 2]
        npz = np.load(v["cache_npz"])
        ax_g1.imshow(npz["g1_33"], extent=ext, origin="upper")
        ax_pt.imshow(npz["lrm33"], extent=ext, origin="upper", cmap="inferno", vmin=-1.5, vmax=5.0)
        ax_pt.scatter(npz["px"], npz["py"], s=0.6, c="#00ffcc", alpha=0.55, linewidths=0)
        ax_sat.imshow(npz["sat_1800"])
        ax_sat.axhline(160, color="#ffd166", lw=0.7, ls="--")
        ax_sat.axvline(160, color="#ffd166", lw=0.7, ls="--")

        nk = v["nearest_any_known"]
        ax_g1.set_title(f"#{idx+1:02d} {v['strip_id'].replace('ChiapsCampeche_NFI_','CC_').replace('Chiaps1_','C1_')}\n({v['lat']:.5f}, {v['lon']:.5f}) d={v['dist_from_card08_km']}km", color="#ffd166", fontsize=7.5, fontweight="bold")
        ax_pt.set_title(f"M:{v['n_mounds']}(R:{v['n_rect']}) H:{v['h1_m']}m Plz:{v['plaza_m2']:.0f}m²\nLAZ:{v['laz_pts_270m']:,} ({v['mound_Pts_str'] if 'mound_Pts_str' in v else str(v['mound_laz_pts'])+' mnd pts'})", color="#4ade80", fontsize=7.5, fontweight="bold")
        ax_sat.set_title(f"Forest:{v['sat_1800m_forest_pct']:.0f}% CHM:{v['chm_med_m']:.0f}m\n{nk['name'][:20]} ({nk['dist_km']}km)", color="#93c5fd", fontsize=7.2, fontweight="bold")
        for ax in (ax_g1, ax_pt, ax_sat):
            ax.set_xticks([])
            ax.set_yticks([])

    for idx in range(len(top_vis), nrows * ncols):
        r_i, c_i = divmod(idx, ncols)
        for k in range(3):
            axes[r_i, c_i * 3 + k].axis("off")

    fig.tight_layout(pad=0.6)
    screen_png = OUT_DIR / "diag_lacandon_corridor_screen24.png"
    fig.savefig(screen_png, facecolor="#0b0f17")
    plt.close(fig)
    print("Saved visual screening sheet:", screen_png, flush=True)


if __name__ == "__main__":
    main()
