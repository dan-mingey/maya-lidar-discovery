#!/usr/bin/env python3
"""Render Hunt-v1-matching 3D LiDAR Canopy-Strip & 360° Orbit GIFs + MP4s
for Open LiDAR Hunt v2 (Motagua UAV + Amazon EBA Acre).

Matches the exact 3D visual grammar of NC-1..NC-6 (`render_3d_video_simulations.py`):
  - Phase 1: Optical Satellite + Real 3D Laser Forest Canopy Points (`.laz` / EPT)
  - Phase 2: Diagonal Cyan Airborne LiDAR Laser Sweep Stripping Jungle Canopy
  - Phase 3: 360° Oblique Orbit of Solid Cutaway Bare-Earth Block (+Z UP, camera-locked
             top-left sun lighting eliminating crater/mound relief inversion, 3D summit/ditch
             telemetry pin, and dual 2D Bare-Earth + Satellite HUD with live camera vector)
  - Explicit WGS84 coordinates in header & bottom telemetry bar.

Outputs (GIF + MP4 written to `media/open-lidar-hunt-v2/` in store + repo):
  1. `sim-v2-qc1-quirigua-acropolis-control.{gif,mp4}`
  2. `sim-v2-ac1-acre-ring-geoglyph-control.{gif,mp4}`
  3. `sim-v2-mg-gualan-a-terrace-lead.{gif,mp4}`
  4. `sim-v2-am05-acre-canopy-audit.{gif,mp4}`
  5. `sim-v2-mg-riotambor-hills-audit.{gif,mp4}`
  6. `sim-v2-all-showcase.mp4` (concatenated master showcase)
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2_scene_tools as vst

WORK = vst.WORK
VID_WORK = WORK / "video_sims_v2"
VID_WORK.mkdir(parents=True, exist_ok=True)
STORE_OUT = vst.STORE_MEDIA
REPO_OUT = vst.REPO_MEDIA


SITES = [
    {
        "code": "QC-1",
        "slug_base": "sim-v2-qc1-quirigua-acropolis-control",
        "kind": "motagua",
        "src": "quirigua_museo_dtm_cog.tif",
        "ept_idx": 7,
        "lat": 15.26885,
        "lon": -89.04025,
        "half_m": 190.0,
        "z_exag": 2.1,
        "mode": "mound",
        "pin_xy": (-42.0, -48.0),
        "peak_label": "UNESCO Southern Acropolis (+11.2 m Masonry)",
        "title": "QC-1 POSITIVE CONTROL: Quiriguá UNESCO Maya Acropolis & Great Plaza",
        "region": "Izabal, Guatemala (15.26885° N, 89.04025° W) • Missouri S&T Fault UAV LiDAR",
        "stats": "0.25 m DTM • 166 pts/m² • +11.2 m Masonry Acropolis & Ballcourt under 35 m Rainforest",
        "lit_status": "VERIFIED POSITIVE CONTROL (UNESCO World Heritage Core — Confirms Sub-Canopy Pipeline)",
        "badge_col": (74, 222, 128),
    },
    {
        "code": "AC-1",
        "slug_base": "sim-v2-ac1-acre-ring-geoglyph-control",
        "kind": "eba",
        "src": WORK / "amazon_eba/laz/NP_T-0414.laz",
        "ept_idx": 0,
        "lat": -9.85446,
        "lon": -67.27142,
        "half_m": 220.0,
        "z_exag": 3.4,
        "mode": "ditch",
        "pin_xy": (68.0, -28.0),
        "peak_label": "267 m Circular Geoglyph Ditch (-2.2 m) & Radial Roads",
        "title": "AC-1 POSITIVE CONTROL: 267 m Circular Ditch Geoglyph & Mounded Village",
        "region": "Acre, SW Amazonia, Brazil (-9.85446° S, -67.27142° W) • INPE EBA Strip NP_T-0414",
        "stats": "1.0 m Circular-Disk PMF Bare-Earth DTM • -2.2 m Ring Ditch + Outer Bank + Radial Causeways",
        "lit_status": "VERIFIED AMAZON CONTROL (Jacobs acrc33 + mouv65 / Peripato 2023 — Confirms Ditch Detector)",
        "badge_col": (74, 222, 128),
    },
    {
        "code": "MG-GUALAN-A",
        "slug_base": "sim-v2-mg-gualan-a-terrace-lead",
        "kind": "motagua",
        "src": "gualan_dtm_cog.tif",
        "ept_idx": 8,
        "lat": 15.13427,
        "lon": -89.31641,
        "half_m": 165.0,
        "z_exag": 2.2,
        "mode": "mound",
        "pin_xy": (-22.0, 12.0),
        "peak_label": "Rectilinear Terrace Platform (+3.5–6.0 m above Channel)",
        "title": "MG-Gualán A: Elevated River-Terrace Platform & Mound Complex",
        "region": "Gualán, Zacapa, Guatemala (15.13427° N, 89.31641° W) • 33 km WSW of Quiriguá",
        "stats": "0.50 m UAV DTM (119 pts/m²) • Rectilinear Platform Overlooking Stream under 21 m Canopy",
        "lit_status": "STRONGEST UNCATALOGUED HUNT v2 LEAD • Survives 3D Sub-Canopy Tree Removal • Needs Field Check",
        "badge_col": (250, 204, 21),
    },
    {
        "code": "AM-05-AUDIT",
        "slug_base": "sim-v2-am05-acre-canopy-audit",
        "kind": "eba",
        "src": WORK / "amazon_eba/laz/NP_T-0414.laz",
        "ept_idx": 0,
        "lat": -9.89700,
        "lon": -67.30296,
        "half_m": 200.0,
        "z_exag": 2.4,
        "mode": "audit",
        "pin_xy": (15.0, 0.0),
        "peak_label": "True Ground Is Smooth! (22 m Forest Edge Mislabeled Class 2)",
        "title": "AM-05 FORENSIC AUDIT: Why Pass-1 Looked Rectangular (Canopy-Edge Artifact)",
        "region": "Acre, Brazil (-9.89700° S, -67.30296° W) • INPE EBA Strip NP_T-0414",
        "stats": "West = 0 m Pasture • East = 22 m Rainforest Mislabeled 'Class 2' in Raw LAZ • <0.5 m True Relief",
        "lit_status": "AUDIT VERDICT: REJECTED FALSE POSITIVE — Stripping 22 m Trees Reveals Natural Smooth Slope",
        "badge_col": (251, 146, 60),
    },
    {
        "code": "MG-RIOTAMBOR",
        "slug_base": "sim-v2-mg-riotambor-hills-audit",
        "kind": "motagua",
        "src": "rio_tambor_dtm_cog.tif",
        "ept_idx": 9,
        "lat": 14.93932,
        "lon": -89.81672,
        "half_m": 165.0,
        "z_exag": 2.0,
        "mode": "mound",
        "pin_xy": (0.0, 0.0),
        "peak_label": "Dissected Ridge Knoll (~21 km E of Taube Jade Zone)",
        "title": "MG-Río Tambor Audit: Dissected Forested Spurs in Motagua Jade Hinterland",
        "region": "Zacapa, Guatemala (14.93932° N, 89.81672° W) • ~21 km E of Carrizal Grande Jade Sites",
        "stats": "0.50 m UAV DTM (165 pts/m²) • Eroded Ridge Spurs under 17 m Dry-Forest Canopy",
        "lit_status": "AUDIT VERDICT: AMBIGUOUS / MOSTLY NATURAL TECTONIC SPURS (Lacks Level Maya Courtyard)",
        "badge_col": (248, 113, 113),
    },
]


def load_fonts():
    try:
        f_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 21)
        f_sub = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14)
        f_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 13)
        f_mono = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf", 12)
    except Exception:
        f_title = f_sub = f_small = f_mono = ImageFont.load_default()
    return f_title, f_sub, f_small, f_mono


def compute_albedo(lrm_s: np.ndarray, g1_rgb: np.ndarray, mode: str) -> np.ndarray:
    if mode == "ditch":
        # Highlight negative circular geoglyph ditches in crisp sapphire-cyan,
        # neutral plateau in warm olive-sand, and positive embankments/roads in terracotta-gold!
        d_neg = np.clip(-lrm_s / 1.8, 0.0, 1.0)
        d_pos = np.clip(lrm_s / 1.4, 0.0, 1.0)
        r = 0.56 - 0.42 * d_neg + 0.38 * d_pos
        g = 0.58 - 0.14 * d_neg + 0.14 * d_pos
        b = 0.46 + 0.46 * d_neg - 0.22 * d_pos
        albedo = np.stack([r, g, b], axis=-1).astype(np.float32)
    else:
        t = np.clip((lrm_s + 0.6) / 5.2, 0.0, 1.0)
        r = np.where(t < 0.35, 0.28 + 0.55 * (t / 0.35), 0.83 + 0.17 * ((t - 0.35) / 0.65))
        g = np.where(t < 0.35, 0.38 + 0.42 * (t / 0.35), 0.80 - 0.32 * ((t - 0.35) / 0.65))
        b = np.where(t < 0.35, 0.36 + 0.18 * (t / 0.35), 0.54 - 0.36 * ((t - 0.35) / 0.65))
        albedo = np.stack([r, g, b], axis=-1).astype(np.float32)

    g1_lum = (
        0.299 * g1_rgb[:, :, 0] + 0.587 * g1_rgb[:, :, 1] + 0.114 * g1_rgb[:, :, 2]
    ).astype(np.float32) / 255.0
    g1_detail = np.clip(0.76 + 0.42 * (g1_lum - 0.5), 0.52, 1.16)[..., None]
    return np.clip(albedo * g1_detail, 0.05, 1.0)


def prepare_site_scene(site: dict) -> dict:
    if site["kind"] == "motagua":
        b = vst.fetch_motagua_bundle(
            site["code"], site["src"], site["ept_idx"], site["lat"], site["lon"], half_m=site["half_m"]
        )
    else:
        b = vst.fetch_eba_acre_bundle(
            site["code"], site["src"], site["lat"], site["lon"], half_m=site["half_m"]
        )

    dem = b["dem"].copy()
    valid = np.isfinite(dem)
    res_m = b["res_m"]
    half_m = site["half_m"]

    base_z = float(np.nanpercentile(dem[valid], 4))
    dem_fill = np.where(valid, dem, base_z)
    dem_s = ndi.gaussian_filter(dem_fill, sigma=0.75)

    lrm_s = ndi.gaussian_filter(np.nan_to_num(b["lrm"], nan=0.0), sigma=0.75)

    if site["mode"] == "ditch":
        # Amplify both negative geoglyph ditch and positive outer bank so 2.2m ditch is dramatic in 3D
        z_rel_true = np.maximum(0.0, dem_s - base_z)
        z_surf_grid = np.maximum(0.6, 2.5 + 0.65 * z_rel_true + 1.55 * lrm_s)
    elif site["mode"] == "audit":
        # Honest 1:1 true topography above base datum (no artificial boost) so flat ground looks flat!
        z_surf_grid = 1.8 + 0.80 * np.maximum(0.0, dem_s - base_z)
    else:
        z_rel_true = np.maximum(0.0, dem_s - base_z)
        z_struct_boost = np.maximum(0.0, lrm_s) * 1.15
        z_surf_grid = 1.8 + 0.65 * z_rel_true + z_struct_boost

    dz_drow, dz_dcol = np.gradient(z_surf_grid * site["z_exag"], res_m, res_m)
    dz_dx = dz_dcol
    dz_dy = -dz_drow
    norm_len = np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy + 1.0)
    nx_grid = (-dz_dx / norm_len).astype(np.float32)
    ny_grid = (-dz_dy / norm_len).astype(np.float32)
    nz_grid = (1.0 / norm_len).astype(np.float32)

    albedo_grid = compute_albedo(lrm_s, b["g1_rgb"], site["mode"])

    hs_f = np.clip(np.nan_to_num(b["hs"], nan=0.5), 0.0, 1.0)[..., None]
    hud_bare_arr = np.clip(albedo_grid * (0.36 + 0.80 * hs_f) * 255.0, 0, 255).astype(np.uint8)
    hud_bare_arr[~valid] = np.array([14, 18, 24], dtype=np.uint8)

    sat_img = Image.fromarray(b["sat_rgb"]).resize((dem.shape[1], dem.shape[0]), Image.Resampling.BILINEAR)
    sat_rgb = np.asarray(sat_img, dtype=np.uint8)

    H, W = dem.shape
    step = max(1, int(round(min(H, W) / 360.0)))
    yy, xx = np.mgrid[0:H:step, 0:W:step]
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

    # Build 3D cutaway plinth walls along outermost valid pixels on all 4 sides
    skirt_x, skirt_y, skirt_z, skirt_nx, skirt_ny, skirt_nz, skirt_alb = [], [], [], [], [], [], []
    n_z_steps = 14
    for r_idx, c_idx, enx, eny in [
        (H - 1 - step, slice(0, W, step), 0.0, -1.0),
        (0, slice(0, W, step), 0.0, 1.0),
        (slice(0, H, step), 0, -1.0, 0.0),
        (slice(0, H, step), W - 1 - step, 1.0, 0.0),
    ]:
        ex = (
            (np.arange(W)[c_idx] if isinstance(c_idx, slice) else np.full(len(range(0, H, step)), c_idx)).astype(
                np.float32
            )
            / max(W - 1, 1)
            * 2.0
            - 1.0
        ) * half_m
        ey = (
            1.0
            - (np.arange(H)[r_idx] if isinstance(r_idx, slice) else np.full(len(range(0, W, step)), r_idx)).astype(
                np.float32
            )
            / max(H - 1, 1)
            * 2.0
        ) * half_m
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
            band_shade = 0.26 + 0.14 * frac + 0.04 * np.sin(zk * 2.2)
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

    # Pin world coordinate
    pin_wx, pin_wy = site["pin_xy"]
    pc = int(np.clip(((pin_wx / half_m + 1.0) * 0.5) * (W - 1), 0, W - 1))
    pr = int(np.clip(((1.0 - (pin_wy / half_m + 1.0) * 0.5)) * (H - 1), 0, H - 1))
    peak_z = float(z_surf_grid[pr, pc])

    # Extract real 3D canopy points above ground
    na_x = (b["pts_x"] - b["cx"]).astype(np.float32)
    na_y = (b["pts_y"] - b["cy"]).astype(np.float32)
    agl = b["pts_agl"].astype(np.float32)

    gx = np.clip(((na_x / half_m + 1.0) * 0.5 * (W - 1)).astype(np.int32), 0, W - 1)
    gy = np.clip(((1.0 - (na_y / half_m + 1.0) * 0.5) * (H - 1)).astype(np.int32), 0, H - 1)
    can_mask = (agl >= 1.8) & (agl <= 42.0) & valid[gy, gx]

    can_x = na_x[can_mask]
    can_y = na_y[can_mask]
    can_agl = agl[can_mask]
    can_z = z_surf_grid[gy[can_mask], gx[can_mask]].astype(np.float32) + can_agl * 0.72

    if len(can_x) > 85000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(can_x), size=85000, replace=False)
        can_x, can_y, can_z, can_agl = can_x[idx], can_y[idx], can_z[idx], can_agl[idx]

    c_gx = np.clip(((can_x / half_m + 1.0) * 0.5 * (W - 1)).astype(np.int32), 0, W - 1)
    c_gy = np.clip(((1.0 - (can_y / half_m + 1.0) * 0.5) * (H - 1)).astype(np.int32), 0, H - 1)
    sat_c = sat_rgb[c_gy, c_gx].astype(np.float32)
    h_norm = np.clip((can_agl - 1.8) / 24.0, 0.0, 1.0)[:, None]
    emerald = np.column_stack(
        [26.0 + 58.0 * h_norm[:, 0], 92.0 + 118.0 * h_norm[:, 0], 34.0 + 48.0 * h_norm[:, 0]]
    ).astype(np.float32)
    can_col = np.clip(0.42 * sat_c + 0.58 * emerald, 0, 255).astype(np.uint8)

    hud_bare = Image.fromarray(hud_bare_arr).resize((250, 250), Image.Resampling.BILINEAR)
    hud_sat = Image.fromarray(sat_rgb).resize((250, 250), Image.Resampling.BILINEAR)

    return {
        "site": site,
        "base_z": base_z,
        "peak": (float(pin_wx), float(pin_wy), peak_z),
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


def render_frame(scene: dict, frame_idx: int, total_frames: int, out_path: Path) -> None:
    site = scene["site"]
    half_m = site["half_m"]
    z_exag = site["z_exag"]

    W_img, H_img = 1280, 720
    vp_x0, vp_y0, vp_w, vp_h = 20, 95, 950, 600

    u = frame_idx / max(total_frames - 1, 1)
    az_deg = -35.0 + 330.0 * u
    pitch_deg = 30.0 + 4.5 * math.sin(u * math.pi * 2.0)
    az = math.radians(az_deg)
    pitch = math.radians(pitch_deg)
    ca, sa = math.cos(az), math.sin(az)
    cp, sp = math.cos(pitch), math.sin(pitch)

    if u < 0.15:
        sweep_pos = -1.50
        phase_label = "PHASE 1: OPTICAL SATELLITE & RAW 3D LASER CANOPY (.LAZ / EPT)"
        phase_col = (72, 219, 132)
    elif u < 0.45:
        sweep_pos = -1.45 + 2.95 * ((u - 0.15) / 0.30)
        phase_label = "PHASE 2: AIRBORNE LiDAR LASER SWEEP — STRIPPING VEGETATION"
        phase_col = (0, 230, 255)
    else:
        sweep_pos = 2.00
        phase_label = "PHASE 3: CLEANED BARE-EARTH 3D TOPOGRAPHY & ORBIT (+Z UP)"
        phase_col = (255, 202, 75)

    sun_rx, sun_ry, sun_rz = -0.52, -0.34, 0.78
    sun_wx = sun_rx * ca + sun_ry * sa
    sun_wy = -sun_rx * sa + sun_ry * ca
    sun_wz = sun_rz

    dot_sun = np.clip(
        scene["surf_nx"] * sun_wx + scene["surf_ny"] * sun_wy + scene["surf_nz"] * sun_wz, -0.2, 1.0
    )
    shade = np.clip(0.38 + 0.74 * np.maximum(0.0, dot_sun), 0.25, 1.12)[:, None]
    surf_lit_col = np.clip(scene["surf_alb"] * shade * 255.0, 0, 255).astype(np.uint8)

    sk_dot = np.clip(scene["sk_nx"] * sun_wx + scene["sk_ny"] * sun_wy + scene["sk_nz"] * sun_wz, -0.2, 1.0)
    sk_shade = np.clip(0.42 + 0.58 * np.maximum(0.0, sk_dot), 0.28, 1.0)[:, None]
    sk_lit_col = np.clip(scene["sk_alb"] * sk_shade * 255.0, 0, 255).astype(np.uint8)

    sx, sy, sz = scene["surf_x"], scene["surf_y"], scene["surf_z"] * z_exag
    diag_s = (sx - sy) / (1.4142 * half_m)

    revealed = diag_s <= sweep_pos
    band = np.abs(diag_s - sweep_pos) < 0.038

    s_col = np.where(
        revealed[:, None],
        surf_lit_col,
        (scene["surf_sat_col"].astype(np.float32) * 0.68).astype(np.uint8),
    )
    if np.any(band):
        s_col = s_col.copy()
        s_col[band] = np.array([0, 245, 255], dtype=np.uint8)

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

    z_pivot = float(np.median(sz)) * 0.65

    rx = all_x * ca - all_y * sa
    ry = all_x * sa + all_y * ca
    rz = all_z - z_pivot

    cam_dist = half_m * 3.30
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

    order = np.argsort(p_depth)[::-1]
    px = px[order]
    py = py[order]
    p_col = p_col[order]

    vp_buf = np.zeros((vp_h, vp_w, 3), dtype=np.uint8)
    bg_grad = np.linspace(12, 24, vp_h, dtype=np.uint8)[:, None]
    vp_buf[:, :, 0] = bg_grad
    vp_buf[:, :, 1] = (bg_grad * 1.15).astype(np.uint8)
    vp_buf[:, :, 2] = (bg_grad * 1.45).astype(np.uint8)

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

    # Project 3D feature pin when revealed
    pk_x, pk_y, pk_z_raw = scene["peak"]
    pk_diag = (pk_x - pk_y) / (1.4142 * half_m)
    if pk_diag <= sweep_pos:
        pk_z = pk_z_raw * z_exag
        prx = pk_x * ca - pk_y * sa
        pry = pk_x * sa + pk_y * ca
        prz_base = pk_z - z_pivot
        prz_top = (pk_z + 22.0 * z_exag) - z_pivot

        d_base = cam_dist + pry * cp - prz_base * sp
        sy_base = +pry * sp + prz_base * cp
        d_top = cam_dist + pry * cp - prz_top * sp
        sy_top = +pry * sp + prz_top * cp

        sx_b = int(vp_x0 + vp_w * 0.5 + (prx * fov) / d_base)
        sy_b = int(vp_y0 + cy_screen - (sy_base * fov) / d_base)
        sx_t = int(vp_x0 + vp_w * 0.5 + (prx * fov) / d_top)
        sy_t = int(vp_y0 + cy_screen - (sy_top * fov) / d_top)

        if (vp_x0 + 35 <= sx_b <= vp_x0 + vp_w - 35) and (vp_y0 + 48 <= sy_t <= vp_y0 + vp_h - 55):
            draw.line([sx_b, sy_b, sx_t, sy_t], fill=(255, 230, 80), width=2)
            draw.ellipse([sx_b - 4, sy_b - 4, sx_b + 4, sy_b + 4], fill=(255, 80, 60), outline=(255, 255, 255))
            lbl = f"▼ {site['peak_label']}" if site["mode"] == "ditch" else f"▲ {site['peak_label']}"
            tw = len(lbl) * 7 + 16
            bx0 = max(vp_x0 + 14, min(sx_t - tw // 2, vp_x0 + vp_w - tw - 14))
            by0 = max(vp_y0 + 48, sy_t - 24)
            draw.rectangle([bx0, by0, bx0 + tw, by0 + 22], fill=(15, 23, 36), outline=(255, 215, 60), width=1)
            draw.text((bx0 + 8, by0 + 4), lbl, fill=(255, 235, 110), font=f_small)

    ax_x, ax_y = vp_x0 + vp_w - 65, vp_y0 + 75
    draw.line([ax_x, ax_y + 22, ax_x, ax_y - 18], fill=(255, 215, 70), width=2)
    draw.polygon([(ax_x, ax_y - 24), (ax_x - 5, ax_y - 14), (ax_x + 5, ax_y - 14)], fill=(255, 215, 70))
    draw.text((ax_x - 22, ax_y + 26), "+Z UP (m)", fill=(255, 215, 70), font=f_small)

    draw.rectangle([0, 0, W_img, 84], fill=(16, 22, 32))
    draw.line([0, 84, W_img, 84], fill=(55, 75, 105), width=2)
    draw.text((22, 9), site["title"], fill=(255, 248, 225), font=f_title)
    draw.text((22, 38), f"{site['region']}   |   {site['stats']}", fill=(185, 205, 230), font=f_sub)
    draw.text((22, 60), site["lit_status"], fill=site["badge_col"], font=f_small)

    draw.rectangle([vp_x0, vp_y0, vp_x0 + vp_w, vp_y0 + vp_h], outline=(65, 88, 120), width=2)
    draw.rectangle([vp_x0 + 12, vp_y0 + 12, vp_x0 + 660, vp_y0 + 42], fill=(12, 18, 26))
    draw.text((vp_x0 + 22, vp_y0 + 19), phase_label, fill=phase_col, font=f_small)

    draw.rectangle([vp_x0 + 12, vp_y0 + vp_h - 42, vp_x0 + vp_w - 12, vp_y0 + vp_h - 10], fill=(12, 18, 26))
    telem = (
        f"ID: {site['code']}  |  WGS84: ({site['lat']:.5f}, {site['lon']:.5f})  |  "
        f"Az: {az_deg % 360:5.1f}°  |  Pitch: {pitch_deg:4.1f}°  |  Datum: {scene['base_z']:.0f} m"
    )
    draw.text((vp_x0 + 20, vp_y0 + vp_h - 33), telem, fill=(210, 225, 245), font=f_mono)

    draw.text((hud_x0, hud_y0_top - 20), "CLEANED BARE-EARTH LRM + SLOPE", fill=(255, 210, 90), font=f_small)
    draw.rectangle([hud_x0, hud_y0_top, hud_x0 + 250, hud_y0_top + 250], outline=(80, 100, 130), width=2)

    draw.text((hud_x0, hud_y0_bot - 20), "ESRI OPTICAL SATELLITE CANOPY", fill=(100, 230, 150), font=f_small)
    draw.rectangle([hud_x0, hud_y0_bot, hud_x0 + 250, hud_y0_bot + 250], outline=(80, 100, 130), width=2)

    for hy in (hud_y0_top, hud_y0_bot):
        mcx, mcy = hud_x0 + 125, hy + 125
        cam_wx = -math.sin(az)
        cam_wy = -math.cos(az)
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


def build_site_media(site: dict, total_frames: int = 120, fps_mp4: int = 20) -> tuple[Path, Path]:
    print(f"Preparing 3D scene for {site['code']} ({site['lat']:.5f}, {site['lon']:.5f})...", flush=True)
    scene = prepare_site_scene(site)
    frames_dir = VID_WORK / f"frames_{site['code']}"
    if frames_dir.exists():
        shutil.rmtree(frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    n_workers = min(6, os.cpu_count() or 4)
    chunks = [list(range(w, total_frames, n_workers)) for w in range(n_workers)]
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        list(ex.map(_worker_chunk, [(scene, ch, frames_dir, total_frames) for ch in chunks]))

    out_mp4 = STORE_OUT / f"{site['slug_base']}.mp4"
    out_gif = STORE_OUT / f"{site['slug_base']}.gif"

    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps_mp4),
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
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # High-quality palettegen/paletteuse GIF matching Hunt v1 (840px wide, 12 fps)
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-i",
            str(out_mp4),
            "-vf",
            "fps=12,scale=840:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=192:stats_mode=diff[p];[s1][p]paletteuse=dither=bayer:bayer_scale=5",
            "-loop",
            "0",
            str(out_gif),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    shutil.copy2(out_mp4, REPO_OUT / out_mp4.name)
    shutil.copy2(out_gif, REPO_OUT / out_gif.name)
    shutil.rmtree(frames_dir, ignore_errors=True)
    print(
        f"Wrote {out_gif.name} ({out_gif.stat().st_size / 1e6:.2f} MB) & {out_mp4.name} ({out_mp4.stat().st_size / 1e6:.2f} MB)",
        flush=True,
    )
    return out_gif, out_mp4


def main() -> None:
    mp4s = []
    for s in SITES:
        _, mp4 = build_site_media(s, total_frames=120, fps_mp4=20)
        mp4s.append(str(mp4))

    concat_list = VID_WORK / "concat_v2.txt"
    concat_list.write_text("".join(f"file '{p}'\n" for p in mp4s))
    master_mp4 = STORE_OUT / "sim-v2-all-showcase.mp4"
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
    shutil.copy2(master_mp4, REPO_OUT / master_mp4.name)
    print(f"Wrote master v2 showcase: {master_mp4.name} ({master_mp4.stat().st_size / 1e6:.2f} MB)", flush=True)


if __name__ == "__main__":
    main()
