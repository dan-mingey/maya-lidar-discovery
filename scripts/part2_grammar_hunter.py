#!/usr/bin/env python3
"""Part 2: Gemini 4 Carbon Multi-Structure Maya Architectural Grammar Hunter.

Scans all 456 G-LiHT 1m bare-earth DEMs paired 1-to-1 with all 456 NASA G-LiHT
Canopy Height Models (CHM), hunting specifically for multi-mound orthogonal
Maya plaza/courtyard/pyramid/acropolis complexes beneath intact forest canopy,
while rejecting modern clearings, ranches, roads, and sascabera borrow pits.
Then verifies contemporary forest cover on Esri World Imagery satellite tiles
and computes exact distances to all 5,223 EAAMS Maya sites + INAH SIC zones.
"""
from __future__ import annotations

import csv
import io
import json
import math
import os
import re
import urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image
from rasterio.enums import Resampling
from rasterio.transform import from_bounds, xy as rio_xy
from rasterio.warp import transform as rio_transform
from rasterio.windows import Window, from_bounds as win_from_bounds, transform as win_transform
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

ROOT = Path("/agent/lidar")
WORK = ROOT / "work/carbon_raw_hunt"
CHM_DIR = WORK / "chm/extracted"
SAT_DIR = WORK / "sat_checks"
SAT_DIR.mkdir(parents=True, exist_ok=True)


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(min(1.0, a)))


def latlon_to_xyz_km(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    r = 6371.0088
    lat_r = np.radians(lats)
    lon_r = np.radians(lons)
    cos_lat = np.cos(lat_r)
    return np.column_stack(
        [r * cos_lat * np.cos(lon_r), r * cos_lat * np.sin(lon_r), r * np.sin(lat_r)]
    )


def build_ot_to_chm_pairs() -> list[tuple[str, str]]:
    ot_dems = sorted(
        list((ROOT / "ot/east/dtm").rglob("*.tif"))
        + list((ROOT / "ot/west/dtm").rglob("*.tif"))
    )
    chms = sorted(CHM_DIR.rglob("*.tif"))
    chm_map: dict[tuple[str, str], Path] = {}
    for c in chms:
        m = re.search(r"AMIGACarb_(.+)_(?:Apr|May)2013_(l\d+s\d+)_CHM", c.name, re.I)
        if not m:
            continue
        camp = m.group(1).lower()
        ls = m.group(2).lower()
        norm = (
            camp.replace("yuc_centro_", "centro_")
            .replace("yuc_norte_", "norte_")
            .replace("yuc_north_2_", "north2_")
            .replace("yuc_south_", "south_")
            .replace("out_of_the_yuc_", "outyuc_")
            .replace("chiap_campeche_", "chiapscampeche_")
            .replace("chiaps_1_", "chiaps1_")
        )
        chm_map[(norm, ls)] = c
        if norm.startswith("outyuc_"):
            chm_map[(norm.replace("outyuc_", "out_yuc_"), ls)] = c

    pairs = []
    for d in ot_dems:
        m = re.match(r"^(.+)_(l\d+s\d+)_DEM_", d.name, re.I)
        k = (m.group(1).lower(), m.group(2).lower()) if m else None
        chm_p = str(chm_map[k]) if k in chm_map else ""
        pairs.append((str(d), chm_p))
    return pairs


def norm_gaussian_lrm(dem: np.ndarray, valid: np.ndarray, sigma_px: float) -> np.ndarray:
    """Fast NaN-safe Local Relief Model via 4x-subsampled normalized Gaussian trend."""
    h, w = dem.shape
    step = 4 if sigma_px >= 12.0 else 1
    if step == 1:
        v_f = valid.astype(np.float32)
        d_zero = np.where(valid, dem, 0.0).astype(np.float32)
        num = ndi.gaussian_filter(d_zero, sigma=sigma_px, mode="constant", cval=0.0)
        den = ndi.gaussian_filter(v_f, sigma=sigma_px, mode="constant", cval=0.0)
        trend = np.where(den > 1e-3, num / np.maximum(den, 1e-3), np.nan)
    else:
        sub_v = valid[::step, ::step].astype(np.float32)
        sub_d = np.where(valid[::step, ::step], dem[::step, ::step], 0.0).astype(np.float32)
        sig_s = sigma_px / step
        num = ndi.gaussian_filter(sub_d, sigma=sig_s, mode="constant", cval=0.0)
        den = ndi.gaussian_filter(sub_v, sigma=sig_s, mode="constant", cval=0.0)
        sub_tr = np.where(den > 1e-3, num / np.maximum(den, 1e-3), 0.0).astype(np.float32)
        trend = ndi.zoom(sub_tr, (h / sub_tr.shape[0], w / sub_tr.shape[1]), order=1)[:h, :w]
    lrm = dem - trend
    lrm[~valid] = np.nan
    return lrm


def split_merged_mounds(lrm: np.ndarray, base_mask: np.ndarray) -> tuple[np.ndarray, int]:
    """Split multi-building basal platforms into individual building footprints
    when distinct peaks (>= 1.3m LRM, >= 16m apart) sit on a shared terrace.
    """
    smooth = ndi.gaussian_filter(np.nan_to_num(lrm, nan=0.0), sigma=2.0)
    local_max = (smooth == ndi.maximum_filter(smooth, size=17)) & base_mask & (smooth >= 1.25)
    markers, n_peaks = ndi.label(local_max)
    if n_peaks <= 1:
        return ndi.label(base_mask)
    # Assign each base_mask pixel to nearest peak weighted by proximity
    base_labels, n_base = ndi.label(base_mask)
    out_labels = np.zeros_like(base_labels)
    cur_id = 0
    for b_id, slc in enumerate(ndi.find_objects(base_labels), start=1):
        if slc is None:
            continue
        comp_sub = base_labels[slc] == b_id
        n_px = int(comp_sub.sum())
        if n_px < 40:
            continue
        ys_r, xs_r = np.where(comp_sub)
        ys = ys_r + slc[0].start
        xs = xs_r + slc[1].start
        lm_sub = local_max[slc] & comp_sub
        n_in = int(lm_sub.sum())
        if n_in <= 1:
            cur_id += 1
            out_labels[ys, xs] = cur_id
        else:
            cpy_r, cpx_r = np.where(lm_sub)
            cpy = cpy_r + slc[0].start
            cpx = cpx_r + slc[1].start
            ch = smooth[cpy, cpx]
            dy = ys[:, None] - cpy[None, :]
            dx = xs[:, None] - cpx[None, :]
            dist2 = (dy * dy + dx * dx) / np.maximum(ch[None, :], 1.0)
            nearest = np.argmin(dist2, axis=1)
            for k in range(len(cpy)):
                sel = nearest == k
                if int(sel.sum()) >= 35:
                    cur_id += 1
                    out_labels[ys[sel], xs[sel]] = cur_id
    return out_labels, cur_id


def extract_mounds_in_window(
    dem: np.ndarray,
    lrm18: np.ndarray,
    slope_deg: np.ndarray,
    interior: np.ndarray,
    res_m: float,
) -> list[dict]:
    base_mask = (lrm18 >= 0.75) & interior
    base_mask = ndi.binary_opening(base_mask, iterations=1)
    base_mask = ndi.binary_closing(base_mask, iterations=1)
    if base_mask.sum() < 120:
        return []

    labeled, nlab = split_merged_mounds(lrm18, base_mask)
    if nlab == 0:
        return []

    objs = ndi.find_objects(labeled)
    mounds = []
    px_area = res_m * res_m
    for lab_idx, slc in enumerate(objs, start=1):
        if slc is None:
            continue
        comp_sub = labeled[slc] == lab_idx
        n_px = int(comp_sub.sum())
        area_m2 = n_px * px_area
        if area_m2 < 45.0 or area_m2 > 7500.0:
            continue

        lrm_sub = lrm18[slc][comp_sub]
        lrm_sub = lrm_sub[np.isfinite(lrm_sub)]
        if len(lrm_sub) < 35:
            continue
        h_max = float(np.max(lrm_sub))
        h_p90 = float(np.percentile(lrm_sub, 90))
        if h_max < 0.95 or h_p90 < 0.78:
            continue

        # Check crisp flank slope (natural broad gentle swells have slope_p80 < 6 deg)
        slp_sub = slope_deg[slc][comp_sub]
        slp_p80 = float(np.nanpercentile(slp_sub, 80)) if np.isfinite(slp_sub).any() else 0.0
        if slp_p80 < 7.2:
            continue

        ys_rel, xs_rel = np.where(comp_sub)
        ys = ys_rel + slc[0].start
        xs = xs_rel + slc[1].start

        pts = np.column_stack([xs.astype(np.float64), ys.astype(np.float64)])
        mean = pts.mean(axis=0)
        X = pts - mean
        cov = X.T @ X / max(len(X) - 1, 1)
        evals, evecs = np.linalg.eigh(cov)
        order = np.argsort(evals)[::-1]
        evecs = evecs[:, order]
        rot = X @ evecs
        w_m = max(float(rot[:, 0].max() - rot[:, 0].min() + 1.0) * res_m, 1.0)
        h_m = max(float(rot[:, 1].max() - rot[:, 1].min() + 1.0) * res_m, 1.0)
        if w_m < 8.0 or w_m > 115.0 or h_m < 5.0:
            continue
        aspect = w_m / h_m
        if aspect > 5.2:
            # Too ribbon-like (road embankment or wall, not a building mound)
            continue

        rect_raw = min(1.0, area_m2 / max(w_m * h_m, 1.0))
        # Major axis azimuth in image coords (0..pi)
        theta = float(math.atan2(evecs[1, 0], evecs[0, 0]) % math.pi)

        # Right-angle / straight-wall boundary orientation score on padded local slice
        pad_y0 = max(0, slc[0].start - 4)
        pad_y1 = min(dem.shape[0], slc[0].stop + 4)
        pad_x0 = max(0, slc[1].start - 4)
        pad_x1 = min(dem.shape[1], slc[1].stop + 4)
        c_pad = labeled[pad_y0:pad_y1, pad_x0:pad_x1] == lab_idx
        eroded = ndi.binary_erosion(c_pad, iterations=1)
        ey, ex = np.where(c_pad & ~eroded)
        ra = 0.0
        if len(ex) >= 18:
            dist_f = ndi.distance_transform_edt(~c_pad)
            gy, gx = np.gradient(dist_f.astype(np.float32))
            ang = np.mod(np.arctan2(gy[ey, ex], gx[ey, ex]), np.pi)
            hist, bins = np.histogram(ang, bins=18, range=(0, np.pi))
            top2 = np.argsort(hist)[-2:]
            if hist[top2].sum() >= 10:
                a0 = (bins[top2[0]] + bins[top2[0] + 1]) / 2
                a1 = (bins[top2[1]] + bins[top2[1] + 1]) / 2
                diff = abs(a0 - a1)
                diff = min(diff, np.pi - diff)
                ra = float(max(0.0, 1.0 - abs(diff - (np.pi / 2)) / (np.pi / 4)))
                conc = float(hist[top2].sum() / max(hist.sum(), 1))
                ra *= 0.5 + 0.5 * conc

        mounds.append(
            {
                "cy": float(mean[1]),
                "cx": float(mean[0]),
                "area_m2": area_m2,
                "length_m": w_m,
                "width_m": h_m,
                "aspect": aspect,
                "h_max": h_max,
                "h_p90": h_p90,
                "slp_p80": slp_p80,
                "rect": rect_raw,
                "ra": ra,
                "theta": theta,
                "ys": ys,
                "xs": xs,
            }
        )
    return mounds


def evaluate_complexes_in_window(
    mounds: list[dict],
    dem: np.ndarray,
    lrm18: np.ndarray,
    lrm45: np.ndarray,
    slope_deg: np.ndarray,
    chm: np.ndarray | None,
    interior: np.ndarray,
    res_m: float,
) -> list[dict]:
    if len(mounds) < 3:
        return []

    centers = np.array([[m["cy"], m["cx"]] for m in mounds])
    tree = cKDTree(centers * res_m)

    # Identify quarry pits (sharp negative cuts < -1.15m with steep >18 deg walls)
    pit_mask = (lrm18 < -1.15) & interior
    quarry_px = pit_mask & (ndi.maximum_filter(slope_deg, size=5) > 18.0)

    # Identify candidate aguadas (smooth flat-bottomed basins < -0.65m with gentle floor)
    aguada_mask = (lrm18 < -0.65) & (slope_deg < 4.0) & interior
    ag_lab, n_ag = ndi.label(ndi.binary_opening(aguada_mask, iterations=2))
    ag_sizes = ndi.sum(aguada_mask, ag_lab, index=np.arange(1, n_ag + 1)) if n_ag > 0 else []
    valid_aguada_px = np.zeros_like(interior)
    for idx_a, sz in enumerate(ag_sizes, start=1):
        if 220 <= sz * res_m * res_m <= 15000:
            valid_aguada_px |= ag_lab == idx_a

    # Identify linear causeway/sacbe hints (elongated 0.35..1.4m ridges)
    low_ridge = (lrm18 >= 0.35) & (lrm18 <= 1.6) & interior
    lr_lab, n_lr = ndi.label(low_ridge)
    sacbe_hint_centers = []
    if n_lr > 0:
        for idx_l, slc in enumerate(ndi.find_objects(lr_lab), start=1):
            if slc is None:
                continue
            h_s = (slc[0].stop - slc[0].start) * res_m
            w_s = (slc[1].stop - slc[1].start) * res_m
            if max(h_s, w_s) < 55.0:
                continue
            c_l = lr_lab[slc] == idx_l
            area_l = c_l.sum() * res_m * res_m
            if area_l < 250 or area_l > 9000:
                continue
            ys_l, xs_l = np.where(c_l)
            pts_l = np.column_stack([xs_l.astype(np.float64), ys_l.astype(np.float64)])
            ev, _ = np.linalg.eigh( np.cov(pts_l.T) )
            elong = math.sqrt(max(ev[1], 1e-3) / max(ev[0], 1e-3))
            if elong >= 6.5:
                sacbe_hint_centers.append(
                    ((slc[0].start + ys_l.mean()), (slc[1].start + xs_l.mean()), max(h_s, w_s))
                )

    # Precompute window-wide mound union & leveled plaza patches ONCE per window!
    all_mound_union = np.zeros_like(interior)
    for m in mounds:
        all_mound_union[m["ys"], m["xs"]] = True
    win_plaza_px = (
        interior
        & (~all_mound_union)
        & (np.abs(lrm18) <= 0.42)
        & (slope_deg <= 4.2)
    )
    win_plaza_px = ndi.binary_opening(win_plaza_px, iterations=2)
    win_pl_lab, win_n_pl = ndi.label(win_plaza_px)
    plaza_patches = []
    if win_n_pl > 0:
        for pid, p_slc in enumerate(ndi.find_objects(win_pl_lab), start=1):
            if p_slc is None:
                continue
            pm = win_pl_lab[p_slc] == pid
            p_area = float(pm.sum() * res_m * res_m)
            if p_area < 180.0:
                continue
            pys, pxs = np.where(pm)
            pcy = p_slc[0].start + float(pys.mean())
            pcx = p_slc[1].start + float(pxs.mean())
            pslp = float(np.nanmedian(slope_deg[p_slc][pm]))
            plaza_patches.append((pcy, pcx, p_area, pslp))

    # Deduplicate anchors within 65m (keep tallest mound in each local neighborhood)
    anchor_order = sorted(range(len(mounds)), key=lambda idx: mounds[idx]["h_max"], reverse=True)
    used_anchor = np.zeros(len(mounds), dtype=bool)

    complexes = []
    for i in anchor_order:
        if used_anchor[i]:
            continue
        anchor = mounds[i]
        if anchor["h_max"] < 1.35:
            break
        close_anchors = tree.query_ball_point(centers[i] * res_m, r=65.0)
        for ca in close_anchors:
            used_anchor[ca] = True
        nbr_idx = tree.query_ball_point(centers[i] * res_m, r=115.0)
        if len(nbr_idx) < 3:
            continue
        grp = [mounds[j] for j in nbr_idx]

        # Require Swath Completeness around group center
        g_cy = float(np.average([m["cy"] for m in grp], weights=[m["area_m2"] for m in grp]))
        g_cx = float(np.average([m["cx"] for m in grp], weights=[m["area_m2"] for m in grp]))
        rad_px = int(90.0 / res_m)
        y0, y1 = max(0, int(g_cy) - rad_px), min(dem.shape[0], int(g_cy) + rad_px)
        x0, x1 = max(0, int(g_cx) - rad_px), min(dem.shape[1], int(g_cx) + rad_px)
        if (y1 - y0) < rad_px or (x1 - x0) < rad_px:
            continue
        swath_fill = float(interior[y0:y1, x0:x1].mean())
        if swath_fill < 0.55:
            continue

        # 1. Canopy Height Model (CHM) check (2013 NASA G-LiHT)
        chm_med = -1.0
        chm_p90 = -1.0
        chm_mound_med = -1.0
        chm_cleared_frac = 1.0
        if chm is not None:
            chm_win = chm[y0:y1, x0:x1]
            c_val = chm_win[np.isfinite(chm_win) & (chm_win >= -1.0) & (chm_win < 60.0)]
            if c_val.size >= 400:
                chm_med = float(np.median(c_val))
                chm_p90 = float(np.percentile(c_val, 90))
                chm_cleared_frac = float((c_val < 3.0).mean())
                m_chm_vals = []
                for m in grp:
                    mv = chm[m["ys"], m["xs"]]
                    mv = mv[np.isfinite(mv) & (mv >= -1.0) & (mv < 60.0)]
                    if mv.size > 0:
                        m_chm_vals.append(mv)
                if m_chm_vals:
                    all_mc = np.concatenate(m_chm_vals)
                    chm_mound_med = float(np.median(all_mc))

        # 2. Orthogonal Grid Coherence (modulo 90 degrees -> 4*theta circular resultant)
        weights = np.array([m["rect"] * math.sqrt(m["area_m2"]) * max(1.0, m["aspect"] - 0.5) for m in grp])
        thetas = np.array([m["theta"] for m in grp])
        c4 = float(np.sum(weights * np.cos(4.0 * thetas)) / max(np.sum(weights), 1e-6))
        s4 = float(np.sum(weights * np.sin(4.0 * thetas)) / max(np.sum(weights), 1e-6))
        ortho_coherence = math.hypot(c4, s4)
        dom_phi = 0.25 * math.atan2(s4, c4)  # dominant orthogonal grid orientation

        # 3. Enclosed Leveled Plaza / Courtyard Test (using precomputed plaza_patches)
        int_patch = interior[y0:y1, x0:x1]
        plaza_area_m2 = 0.0
        plaza_slope_med = 99.0
        pl_cy, pl_cx = g_cy, g_cx
        best_pl_score = -1.0
        for pcy, pcx, p_area, pslp in plaza_patches:
            d_cent = math.hypot(pcy - g_cy, pcx - g_cx) * res_m
            if d_cent > 68.0:
                continue
            ps = min(p_area, 3500.0) * math.exp(-d_cent / 55.0)
            if ps > best_pl_score:
                best_pl_score = ps
                plaza_area_m2 = p_area
                plaza_slope_med = pslp
                pl_cy, pl_cx = pcy, pcx

        # Measure how many of the 4 orthogonal plaza sides (+u, -u, +v, -v) are bounded by mounds
        cos_p, sin_p = math.cos(dom_p := dom_phi), math.sin(dom_p)
        quadrants = set()
        for m in grp:
            dy_m = (m["cy"] - pl_cy) * res_m
            dx_m = (m["cx"] - pl_cx) * res_m
            dist_pl = math.hypot(dx_m, dy_m)
            if 10.0 <= dist_pl <= 90.0:
                u = dx_m * cos_p + dy_m * sin_p
                v = -dx_m * sin_p + dy_m * cos_p
                if abs(u) >= abs(v):
                    quadrants.add(0 if u > 0 else 2)
                else:
                    quadrants.add(1 if v > 0 else 3)
        n_quadrants = len(quadrants)

        # 4. Hierarchy, Rectilinearity, Basal Acropolis, Quarry Penalty
        heights = sorted([m["h_max"] for m in grp], reverse=True)
        max_h = heights[0]
        sec_h = heights[1]
        third_h = heights[2]
        n_rect = sum(1 for m in grp if m["rect"] >= 0.56 and m["ra"] >= 0.35)
        mean_rect = float(np.mean([m["rect"] for m in grp[:5]]))
        mean_ra = float(np.mean([m["ra"] for m in grp[:5]]))

        lrm45_win = lrm45[y0:y1, x0:x1]
        basal_acropolis_m = (
            float(np.nanpercentile(lrm45_win[int_patch], 85)) if int_patch.sum() > 100 else 0.0
        )
        quarry_frac = float(quarry_px[y0:y1, x0:x1].mean())
        has_aguada = bool(valid_aguada_px[y0:y1, x0:x1].any())
        has_sacbe = any(
            math.hypot(sy - g_cy, sx - g_cx) * res_m <= 130.0 for sy, sx, _ in sacbe_hint_centers
        )

        # Require genuine multi-structure architectural qualities
        if n_quadrants < 2 or mean_rect < 0.48:
            continue
        # Reject obvious roadside quarry / spoil heaps
        if quarry_frac > 0.025:
            continue

        # Canopy intactness multiplier: strongly reward closed mature jungle (>= 13m, <8% open)
        if chm_med < 0:
            canopy_mult = 0.7
        elif chm_med < 8.5 or chm_cleared_frac > 0.20:
            # Cleared / pasture / ranch / village -> heavily suppressed (kept separately only for stats)
            canopy_mult = 0.15
        elif chm_med < 12.0 or chm_cleared_frac > 0.10:
            canopy_mult = 0.65
        else:
            canopy_mult = 1.0 + 0.18 * min(1.0, (chm_med - 12.0) / 10.0) * (1.0 - chm_cleared_frac * 4.0)

        grammar_raw = (
            2.2 * min(len(grp), 8)
            + 2.8 * n_quadrants
            + 5.5 * ortho_coherence
            + 4.5 * mean_rect
            + 3.0 * mean_ra
            + 1.8 * min(max_h, 9.0)
            + 1.4 * min(sec_h, 6.0)
            + 1.0 * min(third_h, 4.0)
            + 2.8 * min(1.0, plaza_area_m2 / 1200.0)
            + 1.2 * min(2.5, max(0.0, basal_acropolis_m))
            + (1.5 if has_aguada else 0.0)
            + (1.5 if has_sacbe else 0.0)
            - 40.0 * quarry_frac
        )
        final_score = grammar_raw * canopy_mult

        complexes.append(
            {
                "score": round(float(final_score), 3),
                "grammar_raw": round(float(grammar_raw), 3),
                "canopy_mult": round(float(canopy_mult), 3),
                "win_cy": g_cy,
                "win_cx": g_cx,
                "plaza_cy": pl_cy,
                "plaza_cx": pl_cx,
                "n_mounds": len(grp),
                "n_rect_mounds": n_rect,
                "n_quadrants": n_quadrants,
                "ortho_coherence": round(ortho_coherence, 3),
                "dom_grid_deg": round(math.degrees(dom_phi) % 90.0, 1),
                "mean_rect": round(mean_rect, 3),
                "mean_ra": round(mean_ra, 3),
                "max_height_m": round(max_h, 2),
                "second_height_m": round(sec_h, 2),
                "third_height_m": round(third_h, 2),
                "plaza_area_m2": round(plaza_area_m2, 1),
                "plaza_slope_deg": round(plaza_slope_med, 2),
                "basal_acropolis_m": round(basal_acropolis_m, 2),
                "has_aguada": has_aguada,
                "has_sacbe_hint": has_sacbe,
                "quarry_frac": round(quarry_frac, 4),
                "chm_med_m": round(chm_med, 2),
                "chm_p90_m": round(chm_p90, 2),
                "chm_mound_med_m": round(chm_mound_med, 2),
                "chm_cleared_frac": round(chm_cleared_frac, 4),
            }
        )
    return complexes


def scan_single_transect(args: tuple[str, str]) -> list[dict]:
    dem_path_s, chm_path_s = args
    dem_path = Path(dem_path_s)
    hits_out: list[dict] = []
    try:
        ds = rasterio.open(dem_path)
    except Exception:
        return hits_out

    chm_full = None
    chm_tf = None
    if chm_path_s and Path(chm_path_s).exists():
        try:
            with rasterio.open(chm_path_s) as chm_ds:
                chm_full = chm_ds.read(1).astype(np.float32)
                if chm_ds.nodata is not None:
                    chm_full = np.where(chm_full == chm_ds.nodata, np.nan, chm_full)
                chm_tf = chm_ds.transform
        except Exception:
            chm_full = None

    res_m = float((abs(ds.res[0]) + abs(ds.res[1])) * 0.5)
    # Coarse mask overview (~64m cells) to jump straight to valid diagonal swath blocks!
    ov_step = 64
    ov_h = max(1, math.ceil(ds.height / ov_step))
    ov_w = max(1, math.ceil(ds.width / ov_step))
    ov = ds.read(1, out_shape=(ov_h, ov_w), resampling=Resampling.nearest)
    ov_valid = np.isfinite(ov) & (np.abs(ov) < 1e6) & (ov > -100)
    if ds.nodata is not None:
        ov_valid &= ov != ds.nodata

    win_sz = 1024
    stride = 680  # 344m overlap ensures no 150m complex is ever split across windows
    r_starts = list(range(0, max(1, ds.height - 256), stride))
    if ds.height > win_sz and (ds.height - win_sz) not in r_starts:
        r_starts.append(ds.height - win_sz)
    c_starts = list(range(0, max(1, ds.width - 256), stride))
    if ds.width > win_sz and (ds.width - win_sz) not in c_starts:
        c_starts.append(ds.width - win_sz)

    for r0 in r_starts:
        r1 = min(ds.height, r0 + win_sz)
        or0, or1 = r0 // ov_step, min(ov_h, math.ceil(r1 / ov_step))
        for c0 in c_starts:
            c1 = min(ds.width, c0 + win_sz)
            oc0, oc1 = c0 // ov_step, min(ov_w, math.ceil(c1 / ov_step))
            if int(ov_valid[or0:or1, oc0:oc1].sum()) < 4:
                continue

            win = Window(c0, r0, c1 - c0, r1 - r0)
            dem = ds.read(1, window=win).astype(np.float32)
            if ds.nodata is not None:
                dem = np.where(dem == ds.nodata, np.nan, dem)
            dem = np.where((np.abs(dem) > 1e6) | (dem < -100), np.nan, dem)
            valid = np.isfinite(dem)
            if int(valid.sum()) < 20000:
                continue

            edge_dist = ndi.distance_transform_edt(valid)
            interior = valid & (edge_dist >= max(10, int(12.0 / res_m)))
            if int(interior.sum()) < 15000:
                continue

            lrm18 = norm_gaussian_lrm(dem, valid, sigma_px=18.0 / res_m)
            lrm45 = norm_gaussian_lrm(dem, valid, sigma_px=45.0 / res_m)

            # Smooth DEM slightly (1m) for clean physical slope in degrees
            vf = valid.astype(np.float32)
            dem_sm = ndi.gaussian_filter(np.where(valid, dem, 0.0), 1.0) / np.maximum(
                ndi.gaussian_filter(vf, 1.0), 1e-3
            )
            gy, gx = np.gradient(dem_sm, res_m, res_m)
            slope_deg = np.degrees(np.arctan(np.hypot(gx, gy))).astype(np.float32)
            slope_deg[~valid] = np.nan

            mounds = extract_mounds_in_window(dem, lrm18, slope_deg, interior, res_m)
            if len(mounds) < 3:
                continue

            w_tf = win_transform(win, ds.transform)
            chm_arr = None
            if chm_full is not None and chm_tf is not None:
                # Exact integer pixel offset in shared 1m WGS84 UTM grid
                dc = int(round((w_tf.c - chm_tf.c) / chm_tf.a))
                dr = int(round((w_tf.f - chm_tf.f) / chm_tf.e))
                chm_arr = np.full(dem.shape, np.nan, dtype=np.float32)
                sr0 = max(0, dr)
                sr1 = min(chm_full.shape[0], dr + dem.shape[0])
                sc0 = max(0, dc)
                sc1 = min(chm_full.shape[1], dc + dem.shape[1])
                if sr1 > sr0 and sc1 > sc0:
                    chm_arr[sr0 - dr : sr1 - dr, sc0 - dc : sc1 - dc] = chm_full[sr0:sr1, sc0:sc1]

            cands = evaluate_complexes_in_window(
                mounds, dem, lrm18, lrm45, slope_deg, chm_arr, interior, res_m
            )
            for c in cands:
                # Keep complexes whose center is not right on window boundary (unless strip edge)
                wy, wx = c["win_cy"], c["win_cx"]
                ux, uy = rio_xy(w_tf, wy, wx)
                lons, lats = rio_transform(ds.crs, "EPSG:4326", [ux], [uy])
                c_out = {k: v for k, v in c.items() if not k.startswith("win_") and not k.startswith("plaza_c")}
                c_out.update(
                    {
                        "strip": dem_path.name,
                        "dem_path": str(dem_path),
                        "chm_path": chm_path_s,
                        "crs": str(ds.crs),
                        "utm_e": round(float(ux), 2),
                        "utm_n": round(float(uy), 2),
                        "lat": round(float(lats[0]), 6),
                        "lon": round(float(lons[0]), 6),
                    }
                )
                hits_out.append(c_out)

    ds.close()
    return hits_out


def nms_complexes(cands: list[dict], radius_m: float = 180.0) -> list[dict]:
    cands = sorted(cands, key=lambda x: x["score"], reverse=True)
    kept: list[dict] = []
    for c in cands:
        if not (math.isfinite(c["lat"]) and math.isfinite(c["lon"])):
            continue
        if any(
            haversine_km(c["lat"], c["lon"], k["lat"], k["lon"]) * 1000.0 < radius_m
            for k in kept
        ):
            continue
        kept.append(c)
    return kept


def check_esri_satellite_forest(cand: dict, half_m: float = 220.0) -> dict:
    """Download Esri World Imagery optical tile around candidate (~440m x 440m)
    and verify contemporary forest canopy vs modern agricultural clearing/ranch/road.
    """
    lat, lon = cand["lat"], cand["lon"]
    dlat = half_m / 110574.0
    dlon = half_m / (111320.0 * math.cos(math.radians(lat)))
    bbox = f"{lon - dlon:.6f},{lat - dlat:.6f},{lon + dlon:.6f},{lat + dlat:.6f}"
    url = (
        "https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/export"
        f"?bbox={bbox}&bboxSR=4326&imageSR=4326&size=384,384&format=png&f=image"
    )
    out_png = SAT_DIR / f"sat_{cand['strip'].split('_DEM')[0]}_{lat:.5f}_{lon:.5f}.png"
    try:
        if not out_png.exists():
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            raw = urllib.request.urlopen(req, timeout=20).read()
            out_png.write_bytes(raw)
        img = np.asarray(Image.open(out_png).convert("RGB"), dtype=np.float32) / 255.0
        # Center 220m x 220m patch (middle 50% of 440m window)
        c = img[96:288, 96:288, :]
        r, g, b = c[:, :, 0], c[:, :, 1], c[:, :, 2]
        lum = 0.299 * r + 0.587 * g + 0.114 * b
        ex_green = 2.0 * g - r - b  # Excess Green index
        # Bare limestone soil / white sascab roads / modern roofs have high brightness & R>=G
        bright_bare_frac = float(((lum > 0.48) & (r >= g * 0.96)).mean())
        # Intact tropical forest canopy is dark-to-medium green (lum < 0.42, g > r, g > b)
        forest_px_frac = float(((lum < 0.45) & (g > r * 1.02) & (g >= b)).mean())
        return {
            "sat_png": str(out_png),
            "sat_forest_frac": round(forest_px_frac, 3),
            "sat_bare_soil_road_frac": round(bright_bare_frac, 4),
            "sat_mean_lum": round(float(lum.mean()), 3),
            "sat_ex_green": round(float(ex_green.mean()), 3),
            "sat_intact_canopy": bool(forest_px_frac >= 0.72 and bright_bare_frac <= 0.03),
        }
    except Exception as e:
        return {"sat_error": str(e), "sat_intact_canopy": False}


def main() -> None:
    pairs = build_ot_to_chm_pairs()
    matched_chm = sum(1 for _, c in pairs if c)
    print(f"Prepared {len(pairs)} OT 1m transects ({matched_chm} paired 1-to-1 with NASA CHM)")

    eaams_gj = json.loads((ROOT / "work/ce_deep/eaams/maya_sites.geojson").read_text())
    eaams = []
    for f in eaams_gj["features"]:
        lon, lat = f["geometry"]["coordinates"][:2]
        p = f.get("properties") or {}
        eaams.append(
            {
                "name": p.get("name") or "Unnamed",
                "rank": p.get("rank"),
                "state": p.get("state") or "",
                "lat": float(lat),
                "lon": float(lon),
                "reference": p.get("reference") or "",
            }
        )
    eaams_tree = cKDTree(
        latlon_to_xyz_km(np.array([s["lat"] for s in eaams]), np.array([s["lon"] for s in eaams]))
    )

    sic = []
    for sp in [ROOT / "work/ce_deep/sic_campeche.csv", ROOT / "work/ce_deep/sic_qroo.csv"]:
        raw_b = sp.read_bytes()
        txt = raw_b.decode("utf-8", errors="ignore") if b"\xe1" not in raw_b else raw_b.decode("latin-1")
        for row in csv.DictReader(txt.splitlines()):
            try:
                sic.append(
                    {
                        "name": row["zona_arqueologica_nombre"].strip(),
                        "state": row["nom_ent"].strip(),
                        "lat": float(row["gmaps_latitud"]),
                        "lon": float(row["gmaps_longitud"]),
                    }
                )
            except Exception:
                pass
    sic_tree = cKDTree(
        latlon_to_xyz_km(np.array([s["lat"] for s in sic]), np.array([s["lon"] for s in sic]))
    )

    all_raw: list[dict] = []
    workers = min(os.cpu_count() or 4, 8)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(scan_single_transect, p): p for p in pairs}
        done = 0
        for fut in as_completed(futs):
            done += 1
            res = fut.result()
            all_raw.extend(res)
            if done % 50 == 0 or done == len(pairs):
                print(f"Scanned {done}/{len(pairs)} transects | raw multi-mound hits={len(all_raw)}", flush=True)

    unique_all = nms_complexes(all_raw, radius_m=180.0)
    print(f"Unique multi-structure complexes after 180m NMS: {len(unique_all)}")

    # Enrich with exact EAAMS & INAH SIC distances
    for c in unique_all:
        q = latlon_to_xyz_km(np.array([c["lat"]]), np.array([c["lon"]]))[0]
        _, idxs = eaams_tree.query(q, k=3)
        near_e = []
        for idx in idxs:
            s = eaams[int(idx)]
            dkm = haversine_km(c["lat"], c["lon"], s["lat"], s["lon"])
            near_e.append(
                {
                    "name": s["name"],
                    "rank": s["rank"],
                    "state": s["state"],
                    "dist_km": round(dkm, 2),
                }
            )
        near_e.sort(key=lambda x: x["dist_km"])
        c["nearest_eaams"] = near_e[0]
        c["top3_eaams"] = near_e

        _, s_idx = sic_tree.query(q, k=1)
        ss = sic[int(s_idx)]
        c["nearest_sic"] = {
            "name": ss["name"],
            "state": ss["state"],
            "dist_km": round(haversine_km(c["lat"], c["lon"], ss["lat"], ss["lon"]), 2),
        }

    # Partition into Under-Canopy Intact Forest vs Cleared/Disturbed, and Remote (>5km, >10km) vs Known
    under_canopy = [
        c
        for c in unique_all
        if c["chm_med_m"] >= 11.0 and c["chm_cleared_frac"] <= 0.10 and c["quarry_frac"] <= 0.012
    ]
    remote_5km = [c for c in under_canopy if c["nearest_eaams"]["dist_km"] >= 5.0]
    remote_10km = [c for c in under_canopy if c["nearest_eaams"]["dist_km"] >= 10.0]
    known_prox = [c for c in under_canopy if c["nearest_eaams"]["dist_km"] <= 2.5]

    print(
        f"Strict Under-Canopy (CHM>=11m, cleared<=10%): {len(under_canopy)} | "
        f">=5km from EAAMS: {len(remote_5km)} | >=10km from EAAMS: {len(remote_10km)} | "
        f"<=2.5km of known EAAMS site: {len(known_prox)}"
    )

    # Verify contemporary satellite forest cover on Esri World Imagery for top 45 remote + top 10 known
    to_sat_check = remote_5km[:45] + known_prox[:10]
    print(f"Fetching Esri World Imagery optical canopy checks for {len(to_sat_check)} top complexes...")
    with ThreadPoolExecutor(max_workers=6) as tex:
        sat_futs = {tex.submit(check_esri_satellite_forest, c): c for c in to_sat_check}
        for fut in as_completed(sat_futs):
            c = sat_futs[fut]
            c["esri_sat"] = fut.result()

    remote_verified_both = [
        c for c in remote_5km[:45] if c.get("esri_sat", {}).get("sat_intact_canopy", False)
    ]
    print(f"Top remote complexes passing BOTH 2013 LiDAR CHM AND contemporary Esri Satellite canopy: {len(remote_verified_both)}")

    out_summary = {
        "counts": {
            "total_transects": len(pairs),
            "matched_chm_transects": matched_chm,
            "unique_multi_mound_complexes_all": len(unique_all),
            "strict_under_canopy_2013_chm": len(under_canopy),
            "remote_gt_5km_under_canopy": len(remote_5km),
            "remote_gt_10km_under_canopy": len(remote_10km),
            "known_within_2p5km_under_canopy": len(known_prox),
            "top45_remote_also_intact_in_modern_esri_sat": len(remote_verified_both),
        },
        "top_remote_verified_complexes": remote_verified_both[:25],
        "top_known_site_complexes": known_prox[:15],
        "all_under_canopy_top150": under_canopy[:150],
    }
    (WORK / "part2_grammar_hunt_results.json").write_text(json.dumps(out_summary, indent=2))

    print("\n=== TOP 15 UNMAPPED (>5 KM FROM ALL 5,223 EAAMS) UNDER-CANOPY MAYA COMPLEXES ===")
    for rank, c in enumerate(remote_verified_both[:15], 1):
        ne = c["nearest_eaams"]
        sat = c.get("esri_sat", {})
        print(
            f"#{rank:02d} score={c['score']:5.2f} | ({c['lat']:.5f}, {c['lon']:.5f}) | {c['strip']} | "
            f"mounds={c['n_mounds']} (rect={c['n_rect_mounds']}, quads={c['n_quadrants']}/4, ortho={c['ortho_coherence']:.2f}) | "
            f"Hmax={c['max_height_m']:.1f}m, H2={c['second_height_m']:.1f}m, plaza={c['plaza_area_m2']:.0f}m2 | "
            f"CHM={c['chm_med_m']:.1f}m (open={c['chm_cleared_frac']*100:.1f}%), SatForest={sat.get('sat_forest_frac',0)*100:.0f}% | "
            f"Nearest EAAMS: {ne['name']} ({ne['dist_km']:.1f} km, {ne['state']})"
        )


if __name__ == "__main__":
    main()
