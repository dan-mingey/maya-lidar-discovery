#!/usr/bin/env python3
"""Render high-frame-rate 3D LiDAR canopy-strip & flyover video simulations
for the 5 newly identified unmapped Maya complexes (NC-1..NC-5) and Candidate C.

Fixes & Enhancements:
  1. Right-handed perspective projection:
       depth  = cam_dist + ry * cos(pitch) - rz * sin(pitch)
       sy_cam = +ry * sin(pitch) + rz * cos(pitch)
     so foreground (ry < 0) is at the bottom of the screen, background (ry > 0)
     is at the top of the screen, and pyramids (+rz > 0) rise strictly UPWARD
     with proper front-over-back occlusion.
  2. Real-time camera-locked top-left sun lighting computed from the 3D surface
     normals (nx, ny, nz) at every frame azimuth, preventing the relief-inversion
     (concave crater vs convex mound) optical illusion during 360° rotation.
  3. Solid 3D architectural cutaway block (vertical plinth walls along all 4 edges
     dropping to a 0 m basal datum) + 3D summit elevation pins (+Z UP).
"""
from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi

sys.path.insert(0, "/agent/lidar/work/carbon_raw_hunt")
import raw_pointcloud_tools as rpt

WORK = Path("/agent/lidar/work/carbon_raw_hunt")
VID_WORK = WORK / "video_sims"
VID_WORK.mkdir(parents=True, exist_ok=True)
MEDIA = Path("/cursor/stores/bc-adff6ec3-2baf-40bb-9bab-102a5f567f18/media")

SITES = [
    {
        "code": "NC-1",
        "slug": "sim-nc1-twin-pyramid-acropolis.mp4",
        "strip_id": "South_NFI_l2s450",
        "lat": 18.28691,
        "lon": -90.01315,
        "half_m": 150.0,
        "z_exag": 2.4,
        "peak_label": "Twin Pyramids (+4.9 m & +4.8 m LRM)",
        "title": "NC-1: Twin-Pyramid Rectangular Acropolis & Courtyard",
        "region": "Central-Southern Campeche • 12.1 km from Nearest Named Site (Esperanza)",
        "stats": "17 Mounds (8 Rectangular) • Twin 4.9 m & 4.8 m Pyramids • 979 m² Leveled Plaza • 16.2 m Closed Forest Canopy",
        "lit_status": "UNMAPPED IN NAMED ATLASES (EAAMS / INAH SIC) • Raw 2013 G-LiHT .las + 2026 0.33 m DEM",
    },
    {
        "code": "NC-2",
        "slug": "sim-nc2-monumental-pyramid.mp4",
        "strip_id": "South_NFI_l2s446",
        "lat": 18.42650,
        "lon": -90.19837,
        "half_m": 150.0,
        "z_exag": 1.9,
        "peak_label": "Monumental Pyramid (+8.5 m LRM / +24.7 m Hilltop)",
        "title": "NC-2: 8.5 m Monumental Flat-Topped Pyramid & Terraced Plaza",
        "region": "Balamkú–Tuxpeña Forest Interior, Campeche • 8.0 km from Nearest Named Site (Tuxpeña)",
        "stats": "#1 Grammar Score (68.01) • 9 Mounds • 8.5 m Square-Based Pyramid (+24.7 m Hilltop Prominence) • 14.2 m Canopy",
        "lit_status": "UNMAPPED IN NAMED ATLASES (EAAMS / INAH SIC) • Raw 2013 G-LiHT .las + 2026 0.33 m DEM",
    },
    {
        "code": "NC-3",
        "slug": "sim-nc3-hilltop-temple-plaza.mp4",
        "strip_id": "Centro_GLAS_l1s465",
        "lat": 19.15266,
        "lon": -89.16487,
        "half_m": 145.0,
        "z_exag": 2.3,
        "peak_label": "Summit Temple (+5.6 m LRM)",
        "title": "NC-3: Hilltop Square Temple-Platform & Leveled Courtyard",
        "region": "Campeche–Quintana Roo Borderlands • 12.1 km from Nearest Named Site (Xpaicheil)",
        "stats": "10 Mounds (4/4 Quadrants) • 5.6 m Summit Temple + 3.8 m Mound • 476 m² Leveled Plaza • 12.2 m Canopy",
        "lit_status": "UNMAPPED IN NAMED ATLASES (EAAMS / INAH SIC) • Raw 2013 G-LiHT .las + 2026 0.33 m DEM",
    },
    {
        "code": "NC-4",
        "slug": "sim-nc4-remote-courtyard.mp4",
        "strip_id": "South_GLAS_l2s437",
        "lat": 18.76349,
        "lon": -90.68341,
        "half_m": 145.0,
        "z_exag": 2.5,
        "peak_label": "Orthogonal Courtyard Mound (+3.4 m LRM)",
        "title": "NC-4: Ultra-Remote Orthogonal Platform & Courtyard Complex",
        "region": "West-Central Campeche Interior • 19.6 km from Nearest Named Site (Lake Pan Caliente)",
        "stats": "11 Mounds (5 Crisp Rectangles, 4/4 Quadrants) • 3.4 m & 2.7 m Platforms • 302 m² Courtyard • 15.4 m Canopy",
        "lit_status": "UNMAPPED IN NAMED ATLASES (EAAMS / INAH SIC) • Raw 2013 G-LiHT .las + 2026 0.33 m DEM",
    },
    {
        "code": "NC-5",
        "slug": "sim-nc5-lshaped-courtyard.mp4",
        "strip_id": "South_NFI_l2s445",
        "lat": 18.50660,
        "lon": -90.24067,
        "half_m": 145.0,
        "z_exag": 2.3,
        "peak_label": "Main Platform (+4.3 m LRM)",
        "title": "NC-5: Rectangular Platform & L-Shaped Courtyard Complex",
        "region": "Southern Campeche Forest Interior • 10.7 km from Nearest Named Site (Chen Laguna)",
        "stats": "12 Mounds (4/4 Quadrants) • 4.3 m & 3.7 m Platforms + L-Shaped Range Buildings • 454 m² Plaza • 12.1 m Canopy",
        "lit_status": "UNMAPPED IN NAMED ATLASES (EAAMS / INAH SIC) • Raw 2013 G-LiHT .las + 2026 0.33 m DEM",
    },
    {
        "code": "CAND-C",
        "slug": "sim-candidate-c-borderland-complex.mp4",
        "strip_id": "Centro_GLAS_l1s463",
        "lat": 18.97127,
        "lon": -89.13886,
        "half_m": 155.0,
        "z_exag": 2.8,
        "peak_label": "40×40 m Platform & W. Acropolis (+3.6 m LRM)",
        "title": "Candidate C: 40×40 m Rimmed Platform & Western Acropolis",
        "region": "W. Bacalar Borderlands (~830 m E of Campeche) • 8.6 km from Sacbecan • 20.3 km ENE of Valeriana",
        "stats": "650×300 m Settlement Zone • 40×40 m Square Platform + 140×110 m W. Acropolis • 6–19 m Forest N of Dirt Road",
        "lit_status": "UNMAPPED IN NAMED ATLASES (EAAMS / INAH SIC) • Raw 2013 G-LiHT .las + 2026 0.33 m DEM",
    },
]


def load_fonts():
    try:
        f_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 22)
        f_sub = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 15)
        f_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 13)
        f_mono = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf", 13)
    except Exception:
        f_title = f_sub = f_small = f_mono = ImageFont.load_default()
    return f_title, f_sub, f_small, f_mono


def base_albedo_from_lrm(lrm_norm: np.ndarray, g1_rgb: np.ndarray) -> np.ndarray:
    """Compute unshaded base color (albedo) from LRM height + G1 edge detail so
    that directional 3D lighting can be applied dynamically in camera space.
    """
    t = np.clip(lrm_norm, 0.0, 1.0)
    # Warm limestone/sandstone/gold/terracotta ramp for elevated mounds,
    # cool slate-green for flat plaza floors
    r = np.where(t < 0.35, 0.28 + 0.55 * (t / 0.35), 0.83 + 0.17 * ((t - 0.35) / 0.65))
    g = np.where(t < 0.35, 0.38 + 0.42 * (t / 0.35), 0.80 - 0.32 * ((t - 0.35) / 0.65))
    b = np.where(t < 0.35, 0.36 + 0.18 * (t / 0.35), 0.54 - 0.36 * ((t - 0.35) / 0.65))
    albedo = np.stack([r, g, b], axis=-1).astype(np.float32)

    # Blend subtle architectural edge texture from G1 without baked directional shadows
    g1_lum = (
        0.299 * g1_rgb[:, :, 0] + 0.587 * g1_rgb[:, :, 1] + 0.114 * g1_rgb[:, :, 2]
    ).astype(np.float32) / 255.0
    g1_detail = np.clip(0.78 + 0.38 * (g1_lum - 0.5), 0.55, 1.15)[..., None]
    return np.clip(albedo * g1_detail, 0.05, 1.0)


def prepare_site_scene(site: dict) -> dict:
    b = rpt.fetch_site_bundle(site["strip_id"], site["lat"], site["lon"], half_m=site["half_m"])
    dem = b["dem33"].copy()
    valid = np.isfinite(dem)
    res_m = b["res_m"]

    # Reference basal datum = 5th percentile of valid ground in window so all terrain sits >= 0 m above plinth base
    base_z = float(np.nanpercentile(dem[valid], 5))
    med_z = float(np.nanmedian(dem[valid]))
    dem_fill = np.where(valid, dem, base_z)
    dem_s = ndi.gaussian_filter(dem_fill, sigma=0.75)

    lrm = b["lrm33"]
    lrm_s = ndi.gaussian_filter(np.nan_to_num(lrm, nan=0.0), sigma=0.75)
    lrm_norm = np.clip((lrm_s + 0.4) / 4.8, 0.0, 1.0)

    # Hybrid elevation for 3D display: blend actual topography (above base_z) with amplified LRM
    # so both the natural hill/ridge AND the crisp flat-topped platforms/pyramids stand out boldly
    z_rel_true = np.maximum(0.0, dem_s - base_z)
    z_struct_boost = np.maximum(0.0, lrm_s) * 1.35
    z_surf_grid = 1.6 + 0.65 * z_rel_true + z_struct_boost

    # Compute 3D surface normals (nx, ny, nz) in world coordinates (x=East, y=North, z=Up)
    # Note: row 0 is North (+y), row H-1 is South (-y), so dz/dy = -gradient_row
    dz_drow, dz_dcol = np.gradient(z_surf_grid * site["z_exag"], res_m, res_m)
    dz_dx = dz_dcol
    dz_dy = -dz_drow
    norm_len = np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy + 1.0)
    nx_grid = (-dz_dx / norm_len).astype(np.float32)
    ny_grid = (-dz_dy / norm_len).astype(np.float32)
    nz_grid = (1.0 / norm_len).astype(np.float32)

    albedo_grid = base_albedo_from_lrm(lrm_norm, b["g1_33"])

    # Also build static 2D HUD image
    hs = b["hs33"]
    hs_f = np.clip(np.nan_to_num(hs, nan=0.5), 0.0, 1.0)[..., None]
    hud_bare_arr = np.clip(albedo_grid * (0.35 + 0.80 * hs_f) * 255.0, 0, 255).astype(np.uint8)
    hud_bare_arr[~valid] = np.array([14, 18, 24], dtype=np.uint8)

    sat_img = Image.fromarray(b["sat_rgb"]).resize((dem.shape[1], dem.shape[0]), Image.Resampling.BILINEAR)
    sat_rgb = np.asarray(sat_img, dtype=np.uint8)

    # Subsample dense surface mesh (~370x370)
    H, W = dem.shape
    step = max(1, int(round(min(H, W) / 370.0)))
    yy, xx = np.mgrid[0:H:step, 0:W:step]
    half_m = site["half_m"]
    wx = ((xx.astype(np.float32) / max(W - 1, 1)) * 2.0 - 1.0) * half_m
    wy = (1.0 - (yy.astype(np.float32) / max(H - 1, 1)) * 2.0) * half_m
    wz = z_surf_grid[yy, xx].astype(np.float32)
    v_mask = valid[yy, xx]

    surf_x = wx[v_mask].ravel()
    surf_y = wy[v_mask].ravel()
    surf_z = wz[v_mask].ravel()
    surf_nx = nx_grid[yy, xx][v_mask].ravel()
    surf_ny = ny_grid[yy, xx][v_mask].ravel()
    surf_nz = nz_grid[yy, xx][v_mask].ravel()
    surf_alb = albedo_grid[yy, xx][v_mask].reshape(-1, 3)
    surf_sat_col = sat_rgb[yy, xx][v_mask].reshape(-1, 3)

    # Build 3D Architectural Block Skirt / Plinth Walls around the perimeter of valid pixels
    # so the terrain sits on a solid raised 3D pedestal with zero concave ambiguity
    skirt_x, skirt_y, skirt_z, skirt_nx, skirt_ny, skirt_nz, skirt_alb = [], [], [], [], [], [], []
    n_z_steps = 14
    for edge_name, r_idx, c_idx, enx, eny in [
        ("south", H - 1 - step, slice(0, W, step), 0.0, -1.0),
        ("north", 0, slice(0, W, step), 0.0, 1.0),
        ("west", slice(0, H, step), 0, -1.0, 0.0),
        ("east", slice(0, H, step), W - 1 - step, 1.0, 0.0),
    ]:
        # For each boundary sample, find the outermost valid pixel if edge is nodata
        ex = ((np.arange(W)[c_idx] if isinstance(c_idx, slice) else np.full(len( range(0, H, step)), c_idx)).astype(np.float32) / max(W - 1, 1) * 2.0 - 1.0) * half_m
        ey = (1.0 - (np.arange(H)[r_idx] if isinstance(r_idx, slice) else np.full(len(range(0, W, step)), r_idx)).astype(np.float32) / max(H - 1, 1) * 2.0) * half_m
        ez_top = z_surf_grid[r_idx, c_idx].astype(np.float32)
        ev = valid[r_idx, c_idx]
        if not np.any(ev):
            continue
        ex, ey, ez_top = ex[ev], ey[ev], ez_top[ev]
        for k in range(n_z_steps):
            frac = k / float(n_z_steps - 1)
            zk = ez_top * frac
            skirt_x.append(ex)
            skirt_y.append(ey)
            skirt_z.append(zk)
            skirt_nx.append(np.full_like(ex, enx))
            skirt_ny.append(np.full_like(ex, eny))
            skirt_nz.append(np.zeros_like(ex))
            # Carved stone plinth color with subtle horizontal strata bands
            band_shade = 0.26 + 0.14 * frac + 0.04 * np.sin(zk * 2.5)
            scol = np.column_stack([band_shade * 0.95, band_shade * 0.88, band_shade * 0.76]).astype(np.float32)
            skirt_alb.append(scol)

    if skirt_x:
        sk_x = np.concatenate(skirt_x)
        sk_y = np.concatenate(skirt_y)
        sk_z = np.concatenate(skirt_z)
        sk_nx = np.concatenate(skirt_nx)
        sk_ny = np.concatenate(skirt_ny)
        sk_nz = np.concatenate(skirt_nz)
        sk_alb = np.vstack(skirt_alb)
    else:
        sk_x = sk_y = sk_z = sk_nx = sk_ny = sk_nz = np.zeros(0, dtype=np.float32)
        sk_alb = np.zeros((0, 3), dtype=np.float32)

    # Locate primary monumental peak in central 60% of window for 3D summit pin
    r0, r1 = int(H * 0.20), int(H * 0.80)
    c0, c1 = int(W * 0.20), int(W * 0.80)
    sub_lrm = np.where(valid[r0:r1, c0:c1], lrm_s[r0:r1, c0:c1], -999.0)
    pr, pc = np.unravel_index(int(np.argmax(sub_lrm)), sub_lrm.shape)
    pr += r0
    pc += c0
    peak_x = float(((pc / max(W - 1, 1)) * 2.0 - 1.0) * half_m)
    peak_y = float((1.0 - (pr / max(H - 1, 1)) * 2.0) * half_m)
    peak_z = float(z_surf_grid[pr, pc])

    # Extract raw NASA .las canopy points above ground
    cx, cy = b["cx"], b["cy"]
    raw_nx, raw_ny, raw_nz, _ = b["nasa_pts"]
    na_x = (raw_nx - cx).astype(np.float32)
    na_y = (raw_ny - cy).astype(np.float32)
    na_z_asl = raw_nz.astype(np.float32)

    gx = np.clip(((na_x / half_m + 1.0) * 0.5 * (W - 1)).astype(np.int32), 0, W - 1)
    gy = np.clip(((1.0 - (na_y / half_m + 1.0) * 0.5) * (H - 1)).astype(np.int32), 0, H - 1)
    agl = na_z_asl - dem_s[gy, gx]
    can_mask = (agl >= 2.2) & (agl <= 36.0) & valid[gy, gx]

    can_x = na_x[can_mask]
    can_y = na_y[can_mask]
    can_agl = agl[can_mask]
    # Place canopy points at their true height above the 3D surface mesh
    can_z = z_surf_grid[gy[can_mask], gx[can_mask]].astype(np.float32) + can_agl * 0.85

    if len(can_x) > 90000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(can_x), size=90000, replace=False)
        can_x, can_y, can_z, can_agl = can_x[idx], can_y[idx], can_z[idx], can_agl[idx]

    c_gx = np.clip(((can_x / half_m + 1.0) * 0.5 * (W - 1)).astype(np.int32), 0, W - 1)
    c_gy = np.clip(((1.0 - (can_y / half_m + 1.0) * 0.5) * (H - 1)).astype(np.int32), 0, H - 1)
    sat_c = sat_rgb[c_gy, c_gx].astype(np.float32)
    h_norm = np.clip((can_agl - 2.0) / 22.0, 0.0, 1.0)[:, None]
    emerald = np.column_stack([
        26.0 + 58.0 * h_norm[:, 0],
        92.0 + 118.0 * h_norm[:, 0],
        34.0 + 48.0 * h_norm[:, 0],
    ]).astype(np.float32)
    can_col = np.clip(0.42 * sat_c + 0.58 * emerald, 0, 255).astype(np.uint8)

    hud_bare = Image.fromarray(hud_bare_arr).resize((250, 250), Image.Resampling.BILINEAR)
    hud_sat = Image.fromarray(sat_rgb).resize((250, 250), Image.Resampling.BILINEAR)

    return {
        "site": site,
        "base_z": base_z,
        "med_z": med_z,
        "peak": (peak_x, peak_y, peak_z),
        "surf_x": surf_x,
        "surf_y": surf_y,
        "surf_z": surf_z,
        "surf_nx": surf_nx,
        "surf_ny": surf_ny,
        "surf_nz": surf_nz,
        "surf_alb": surf_alb,
        "surf_sat_col": surf_sat_col,
        "sk_x": sk_x,
        "sk_y": sk_y,
        "sk_z": sk_z,
        "sk_nx": sk_nx,
        "sk_ny": sk_ny,
        "sk_nz": sk_nz,
        "sk_alb": sk_alb,
        "can_x": can_x,
        "can_y": can_y,
        "can_z": can_z,
        "can_col": can_col,
        "hud_bare": np.asarray(hud_bare),
        "hud_sat": np.asarray(hud_sat),
    }


def render_frame(
    scene: dict,
    frame_idx: int,
    total_frames: int,
    out_path: Path,
) -> None:
    site = scene["site"]
    half_m = site["half_m"]
    z_exag = site["z_exag"]

    W_img, H_img = 1280, 720
    vp_x0, vp_y0, vp_w, vp_h = 20, 95, 950, 600

    u = frame_idx / max(total_frames - 1, 1)

    # Smooth 360° orbit with natural oblique aerial pitch (26° to 34° above horizon)
    az_deg = -35.0 + 330.0 * u
    pitch_deg = 29.0 + 4.5 * math.sin(u * math.pi * 2.0)
    az = math.radians(az_deg)
    pitch = math.radians(pitch_deg)
    ca, sa = math.cos(az), math.sin(az)
    cp, sp = math.cos(pitch), math.sin(pitch)

    # Canopy stripping sweep schedule (diag_c in [-1.42, +1.42]):
    if u < 0.14:
        sweep_pos = -1.50
        phase_label = "PHASE 1: OPTICAL SATELLITE & RAW .LAS FOREST CANOPY (12–20 m TREES)"
        phase_col = (72, 219, 132)
    elif u < 0.44:
        sweep_pos = -1.45 + 2.95 * ((u - 0.14) / 0.30)
        phase_label = "PHASE 2: AIRBORNE LiDAR LASER SWEEP — STRIPPING JUNGLE CANOPY"
        phase_col = (0, 230, 255)
    else:
        sweep_pos = 2.00
        phase_label = "PHASE 3: BARE-EARTH 0.33 m MAYA ARCHITECTURE & 3D ORBIT (+Z UP)"
        phase_col = (255, 202, 75)

    # Camera-locked top-left sun vector in world space so upper-left slopes of every
    # pyramid are always sunlit (eliminating the concave relief-inversion illusion!)
    # In camera space: sun comes from left (-x_cam), slightly front (-y_cam), and high (+z_cam)
    sun_rx, sun_ry, sun_rz = -0.52, -0.34, 0.78
    # Rotate camera-space sun vector back to world space (inverse of rx, ry rotation)
    sun_wx = sun_rx * ca + sun_ry * sa
    sun_wy = -sun_rx * sa + sun_ry * ca
    sun_wz = sun_rz

    # Shade terrain surface with camera-locked sun
    dot_sun = np.clip(
        scene["surf_nx"] * sun_wx + scene["surf_ny"] * sun_wy + scene["surf_nz"] * sun_wz,
        -0.2,
        1.0,
    )
    shade = np.clip(0.38 + 0.74 * np.maximum(0.0, dot_sun), 0.25, 1.12)[:, None]
    surf_lit_col = np.clip(scene["surf_alb"] * shade * 255.0, 0, 255).astype(np.uint8)

    # Shade 3D block skirt walls
    sk_dot = np.clip(
        scene["sk_nx"] * sun_wx + scene["sk_ny"] * sun_wy + scene["sk_nz"] * sun_wz,
        -0.2,
        1.0,
    )
    sk_shade = np.clip(0.42 + 0.58 * np.maximum(0.0, sk_dot), 0.28, 1.0)[:, None]
    sk_lit_col = np.clip(scene["sk_alb"] * sk_shade * 255.0, 0, 255).astype(np.uint8)

    sx, sy, sz = scene["surf_x"], scene["surf_y"], scene["surf_z"] * z_exag
    diag_s = (sx - sy) / (1.4142 * half_m)

    revealed = diag_s <= sweep_pos
    band = np.abs(diag_s - sweep_pos) < 0.038

    s_col = np.where(
        revealed[:, None],
        surf_lit_col,
        (scene["surf_sat_col"].astype(np.float32) * 0.58).astype(np.uint8),
    )
    if np.any(band):
        s_col = s_col.copy()
        s_col[band] = np.array([0, 245, 255], dtype=np.uint8)

    # Combine skirt + surface + unstripped canopy points
    cx, cy, cz = scene["can_x"], scene["can_y"], scene["can_z"] * z_exag
    diag_c = (cx - cy) / (1.4142 * half_m)
    keep_c = diag_c > (sweep_pos - 0.025)

    sk_x, sk_y, sk_z = scene["sk_x"], scene["sk_y"], scene["sk_z"] * z_exag

    if np.any(keep_c):
        cx_k = cx[keep_c]
        cy_k = cy[keep_c]
        cz_k = cz[keep_c]
        cc_k = scene["can_col"][keep_c].copy()
        c_band = np.abs(diag_c[keep_c] - sweep_pos) < 0.035
        if np.any(c_band):
            cc_k[c_band] = np.array([120, 255, 255], dtype=np.uint8)
        all_x = np.concatenate([sk_x, sx, cx_k])
        all_y = np.concatenate([sk_y, sy, cy_k])
        all_z = np.concatenate([sk_z, sz, cz_k])
        all_c = np.vstack([sk_lit_col, s_col, cc_k])
    else:
        all_x = np.concatenate([sk_x, sx])
        all_y = np.concatenate([sk_y, sy])
        all_z = np.concatenate([sk_z, sz])
        all_c = np.vstack([sk_lit_col, s_col])

    # Center vertically around mid-height of terrain block
    z_pivot = float(np.median(sz)) * 0.65

    # Right-handed 3D rotation & perspective projection:
    #   ry < 0 is foreground (closer to camera) -> sy_cam < 0 -> bottom of screen (py larger)
    #   ry > 0 is background (farther from camera) -> sy_cam > 0 -> top of screen (py smaller)
    #   rz > 0 is elevation above ground -> sy_cam > 0 -> rises UPWARD on screen (py smaller)
    rx = all_x * ca - all_y * sa
    ry = all_x * sa + all_y * ca
    rz = all_z - z_pivot

    cam_dist = half_m * 3.25
    depth = cam_dist + ry * cp - rz * sp
    sy_cam = +ry * sp + rz * cp
    fov = vp_w * 1.15
    cy_screen = vp_h * 0.56

    valid_p = depth > 10.0
    px = (vp_w * 0.5 + (rx[valid_p] * fov) / depth[valid_p]).astype(np.int32)
    py = (cy_screen - (sy_cam[valid_p] * fov) / depth[valid_p]).astype(np.int32)
    p_depth = depth[valid_p]
    p_col = all_c[valid_p]

    in_vp = (px >= 1) & (px < vp_w - 2) & (py >= 1) & (py < vp_h - 2)
    px, py, p_depth, p_col = px[in_vp], py[in_vp], p_depth[in_vp], p_col[in_vp]

    # Sort back-to-front (largest depth first, smallest depth last)
    order = np.argsort(p_depth)[::-1]
    px = px[order]
    py = py[order]
    p_col = p_col[order]

    vp_buf = np.zeros((vp_h, vp_w, 3), dtype=np.uint8)
    bg_grad = np.linspace(12, 24, vp_h, dtype=np.uint8)[:, None]
    vp_buf[:, :, 0] = bg_grad
    vp_buf[:, :, 1] = (bg_grad * 1.15).astype(np.uint8)
    vp_buf[:, :, 2] = (bg_grad * 1.45).astype(np.uint8)

    # Splat 2x2 pixels
    vp_buf[py, px] = p_col
    vp_buf[py, px + 1] = p_col
    vp_buf[py + 1, px] = p_col
    vp_buf[py + 1, px + 1] = p_col

    canvas = np.full((H_img, W_img, 3), (11, 15, 22), dtype=np.uint8)
    canvas[vp_y0 : vp_y0 + vp_h, vp_x0 : vp_x0 + vp_w] = vp_buf

    hud_x0 = 995
    hud_y0_top = 115
    hud_y0_bot = 420
    canvas[hud_y0_top : hud_y0_top + 250, hud_x0 : hud_x0 + 250] = scene["hud_bare"]
    canvas[hud_y0_bot : hud_y0_bot + 250, hud_x0 : hud_x0 + 250] = scene["hud_sat"]

    img = Image.fromarray(canvas)
    draw = ImageDraw.Draw(img)
    f_title, f_sub, f_small, f_mono = load_fonts()

    # Project 3D monument peak pin (+Z UP) when peak is revealed
    pk_x, pk_y, pk_z_raw = scene["peak"]
    pk_diag = (pk_x - pk_y) / (1.4142 * half_m)
    if pk_diag <= sweep_pos:
        pk_z = pk_z_raw * z_exag
        prx = pk_x * ca - pk_y * sa
        pry = pk_x * sa + pk_y * ca
        prz_base = pk_z - z_pivot
        prz_top = (pk_z + 18.0 * z_exag) - z_pivot

        d_base = cam_dist + pry * cp - prz_base * sp
        sy_base = +pry * sp + prz_base * cp
        d_top = cam_dist + pry * cp - prz_top * sp
        sy_top = +pry * sp + prz_top * cp

        sx_b = int(vp_x0 + vp_w * 0.5 + (prx * fov) / d_base)
        sy_b = int(vp_y0 + cy_screen - (sy_base * fov) / d_base)
        sx_t = int(vp_x0 + vp_w * 0.5 + (prx * fov) / d_top)
        sy_t = int(vp_y0 + cy_screen - (sy_top * fov) / d_top)

        if (vp_x0 + 40 <= sx_b <= vp_x0 + vp_w - 40) and (vp_y0 + 50 <= sy_t <= vp_y0 + vp_h - 60):
            draw.line([sx_b, sy_b, sx_t, sy_t], fill=(255, 230, 80), width=2)
            draw.ellipse([sx_b - 4, sy_b - 4, sx_b + 4, sy_b + 4], fill=(255, 80, 60), outline=(255, 255, 255))
            lbl = f"▲ {site['peak_label']}"
            tw = len(lbl) * 7 + 14
            bx0 = max(vp_x0 + 14, min(sx_t - tw // 2, vp_x0 + vp_w - tw - 14))
            by0 = max(vp_y0 + 48, sy_t - 24)
            draw.rectangle([bx0, by0, bx0 + tw, by0 + 22], fill=(15, 23, 36), outline=(255, 215, 60), width=1)
            draw.text((bx0 + 7, by0 + 4), lbl, fill=(255, 235, 110), font=f_small)

    # Draw +Z UP vertical reference indicator in top-right of 3D viewport
    ax_x, ax_y = vp_x0 + vp_w - 65, vp_y0 + 75
    draw.line([ax_x, ax_y + 22, ax_x, ax_y - 18], fill=(255, 215, 70), width=2)
    draw.polygon([(ax_x, ax_y - 24), (ax_x - 5, ax_y - 14), (ax_x + 5, ax_y - 14)], fill=(255, 215, 70))
    draw.text((ax_x - 22, ax_y + 26), "+Z UP (ASL)", fill=(255, 215, 70), font=f_small)

    # Header bar
    draw.rectangle([0, 0, W_img, 84], fill=(16, 22, 32))
    draw.line([0, 84, W_img, 84], fill=(55, 75, 105), width=2)
    draw.text((22, 10), site["title"], fill=(255, 248, 225), font=f_title)
    draw.text((22, 39), f"{site['region']}   |   {site['stats']}", fill=(185, 205, 230), font=f_sub)
    draw.text((22, 61), site["lit_status"], fill=(110, 230, 175), font=f_small)

    # Viewport border & phase banner
    draw.rectangle([vp_x0, vp_y0, vp_x0 + vp_w, vp_y0 + vp_h], outline=(65, 88, 120), width=2)
    draw.rectangle([vp_x0 + 12, vp_y0 + 12, vp_x0 + 650, vp_y0 + 42], fill=(12, 18, 26))
    draw.text((vp_x0 + 22, vp_y0 + 19), phase_label, fill=phase_col, font=f_small)

    # Telemetry overlay at bottom of 3D viewport
    draw.rectangle([vp_x0 + 12, vp_y0 + vp_h - 42, vp_x0 + vp_w - 12, vp_y0 + vp_h - 10], fill=(12, 18, 26))
    telem = (
        f"Strip: {site['strip_id']}  |  Azimuth: {az_deg % 360:5.1f}°  |  Pitch: {pitch_deg:4.1f}°  |  "
        f"Base Datum: {scene['base_z']:.0f} m ASL  |  Coords: WITHHELD (INAH)"
    )
    draw.text((vp_x0 + 22, vp_y0 + vp_h - 33), telem, fill=(210, 225, 245), font=f_mono)

    # Right HUD labels & camera viewing vector indicator
    draw.text((hud_x0, hud_y0_top - 20), "0.33 m BARE-EARTH G1 + LRM", fill=(255, 210, 90), font=f_small)
    draw.rectangle([hud_x0, hud_y0_top, hud_x0 + 250, hud_y0_top + 250], outline=(80, 100, 130), width=2)

    draw.text((hud_x0, hud_y0_bot - 20), "ESRI OPTICAL SATELLITE CANOPY", fill=(100, 230, 150), font=f_small)
    draw.rectangle([hud_x0, hud_y0_bot, hud_x0 + 250, hud_y0_bot + 250], outline=(80, 100, 130), width=2)

    for hy in (hud_y0_top, hud_y0_bot):
        mcx, mcy = hud_x0 + 125, hy + 125
        # Camera sits at ry < 0 in camera coords => world (-sa, -ca) looking toward center (mcx, mcy)
        cam_wx = -math.sin(az)
        cam_wy = -math.cos(az)
        # On 2D map: +x is East (right), +y is North (up, so -y in image pixels)
        x_cam = int(mcx + cam_wx * 95)
        y_cam = int(mcy - cam_wy * 95)
        draw.line([x_cam, y_cam, mcx, mcy], fill=(0, 240, 255), width=2)
        draw.ellipse([mcx - 4, mcy - 4, mcx + 4, mcy + 4], outline=(255, 220, 0), width=2)

    img.save(out_path, format="PNG", compress_level=1)


def _worker_chunk(args):
    scene, indices, frames_dir, total_frames = args
    for idx in indices:
        render_frame(scene, idx, total_frames, frames_dir / f"frame_{idx:04d}.png")
    return len(indices)


def build_site_video(site: dict, total_frames: int = 180, fps: int = 24) -> Path:
    print(f"Preparing 3D scene for {site['code']} ({site['strip_id']})...", flush=True)
    scene = prepare_site_scene(site)
    frames_dir = VID_WORK / f"frames_{site['code']}"
    if frames_dir.exists():
        shutil.rmtree(frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    n_workers = min(4, os.cpu_count() or 4)
    chunks = [list(range(w, total_frames, n_workers)) for w in range(n_workers)]
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        list(ex.map(_worker_chunk, [(scene, ch, frames_dir, total_frames) for ch in chunks]))

    out_mp4 = MEDIA / site["slug"]
    cmd = [
        "ffmpeg",
        "-y",
        "-framerate",
        str(fps),
        "-i",
        str(frames_dir / "frame_%04d.png"),
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        "19",
        "-preset",
        "fast",
        "-movflags",
        "+faststart",
        str(out_mp4),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    shutil.rmtree(frames_dir, ignore_errors=True)
    print(f"Wrote {out_mp4} ({out_mp4.stat().st_size / 1e6:.2f} MB)", flush=True)
    return out_mp4


def main():
    outputs = []
    for s in SITES:
        mp4 = build_site_video(s, total_frames=180, fps=24)
        outputs.append(str(mp4))

    concat_list = VID_WORK / "concat_list.txt"
    concat_list.write_text("".join(f"file '{p}'\n" for p in outputs))
    master_mp4 = MEDIA / "sim-all-new-maya-complexes-showcase.mp4"
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-f",
            "concat",
            "-safe",
            "0",
            "-i",
            str(concat_list),
            "-c",
            "copy",
            "-movflags",
            "+faststart",
            str(master_mp4),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    print(f"Wrote master showcase: {master_mp4} ({master_mp4.stat().st_size / 1e6:.2f} MB)", flush=True)


if __name__ == "__main__":
    main()
