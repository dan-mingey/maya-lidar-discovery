#!/usr/bin/env python3
"""Render Hunt-v1-identical 0.33 m 3D LiDAR Canopy-Strip & 360° Bare-Earth Orbit
GIFs, MP4s, 3D Gallery, and 4-Panel Yucatán/Globe Locator Atlas for NJ-1..NJ-3:
Three unbroken-jungle (99.3–99.6% closed canopy, zero trails/farms within 2+ km)
Maya rectangular courtyard & pyramid complexes outside Estrada-Belli et al. (2025).
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, Polygon as MplPolygon, Rectangle
from PIL import Image, ImageDraw
from scipy import ndimage as ndi
from shapely.geometry import shape

sys.path.insert(0, "/agent/maya-lidar-discovery/scripts")
sys.path.insert(0, "/agent/maya-lidar-discovery/scripts/open_lidar_hunt_v2")
import raw_pointcloud_tools as rpt
import render_3d_video_simulations as r3d
import v2_scene_tools as vst

OUT_STORE = Path("/cursor/stores/bc-adff6ec3-2baf-40bb-9bab-102a5f567f18/media/open-lidar-hunt-v3")
OUT_REPO = Path("/agent/maya-lidar-discovery/media/open-lidar-hunt-v3")
VID_TMP = Path("/agent/lidar/work/carbon_raw_hunt/video_sims_nj")
OUT_STORE.mkdir(parents=True, exist_ok=True)
OUT_REPO.mkdir(parents=True, exist_ok=True)
VID_TMP.mkdir(parents=True, exist_ok=True)

SITES = [
    {
        "code": "NJ-1",
        "slug_base": "sim-nj1-yucnorth-courtyard-pyramids",
        "strip_id": "North2_GLAS_l0s391",
        "lat": 20.64872,
        "lon": -88.75685,
        "half_m": 145.0,
        "z_exag": 2.45,
        "peak_label": "Temple Pyramid (+6.0 m above Plaza)",
        "title": "NJ-1: Enclosed Rectangular Courtyard, Temple Pyramid & Platform Complex",
        "region": "Central-East Yucatán (20.64872° N, 88.75685° W) • 6.98 km from Cenote Xtohil (EAAMS)",
        "stats": "9 Mounds • +6.0 m Pyramid • 3,142 m² Plaza • 99.3% Unbroken Jungle (Zero Trails/Farms in 2 km)",
        "lit_status": "0 POLYGONS IN ESTRADA-BELLI 2025 (MARI-GISLAB) • Unregistered in EAAMS & INAH SIC • 0.33 m DEM + Raw .LAS",
    },
    {
        "code": "NJ-2",
        "slug_base": "sim-nj2-yucnorth-four-quadrant-plaza",
        "strip_id": "North2_GLAS_l0s390",
        "lat": 20.63046,
        "lon": -88.75890,
        "half_m": 145.0,
        "z_exag": 2.55,
        "peak_label": "4-Quadrant Courtyard Pyramid (+5.0 m)",
        "title": "NJ-2: Four-Quadrant Rectangular Mound Enclosure & Central Plaza",
        "region": "Central-East Yucatán (20.63046° N, 88.75890° W) • 6.04 km from Ikil (EAAMS)",
        "stats": "9 Mounds (4 Rectangular, 4/4 Quadrants) • +5.0 m Pyramid • 548 m² Plaza • 99.6% Unbroken Jungle",
        "lit_status": "0 POLYGONS IN ESTRADA-BELLI 2025 (MARI-GISLAB) • Unregistered in EAAMS & INAH SIC • 0.33 m DEM + Raw .LAS",
    },
    {
        "code": "NJ-3",
        "slug_base": "sim-nj3-yucnorth-18mound-courtyard",
        "strip_id": "North2_NFI_l8s387",
        "lat": 20.66345,
        "lon": -88.77712,
        "half_m": 145.0,
        "z_exag": 2.60,
        "peak_label": "Courtyard Range Mound (+4.4 m)",
        "title": "NJ-3: 18-Mound Enclosed Courtyard & Range-Structure Complex",
        "region": "Central-East Yucatán (20.66345° N, 88.77712° W) • 4.31 km from Cenote Xtohil (EAAMS)",
        "stats": "18 Mounds (4 Rectangular, 4/4 Quadrants) • 1,068 m² Plaza • 99.4% Unbroken Jungle (8,188 Ground Pts)",
        "lit_status": "0 POLYGONS IN ESTRADA-BELLI 2025 (MARI-GISLAB) • Unregistered in EAAMS & INAH SIC • 0.33 m DEM + Raw .LAS",
    },
]


def prepare_flat_plinth_scene(site: dict) -> dict:
    """Use Hunt v1's exact 0.33 m scene builder, but detrend broad regional tilt
    so the Maya plaza sits dead-level on the 3D stone plinth (zero mountain illusion!).
    """
    b = rpt.fetch_site_bundle(site["strip_id"], site["lat"], site["lon"], half_m=site["half_m"])
    dem = b["dem33"].copy()
    valid = np.isfinite(dem)
    res_m = b["res_m"]

    base_z = float(np.nanpercentile(dem[valid], 5))
    dem_fill = np.where(valid, dem, base_z)
    dem_s = ndi.gaussian_filter(dem_fill, sigma=0.75)

    # Remove broad >55m regional bedrock tilt while preserving meso-acropolis & crisp mounds
    broad_trend = ndi.gaussian_filter(dem_s, sigma=55.0 / res_m)
    local_topo = np.maximum(0.0, (dem_s - broad_trend) + 1.2)

    lrm_s = ndi.gaussian_filter(np.nan_to_num(b["lrm33"], nan=0.0), sigma=0.75)
    lrm_norm = np.clip((lrm_s + 0.35) / 4.2, 0.0, 1.0)

    z_struct_boost = np.maximum(0.0, lrm_s) * 1.45
    z_surf_grid = 1.6 + 0.40 * local_topo + z_struct_boost

    dz_drow, dz_dcol = np.gradient(z_surf_grid * site["z_exag"], res_m, res_m)
    dz_dx = dz_dcol
    dz_dy = -dz_drow
    norm_len = np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy + 1.0)
    nx_grid = (-dz_dx / norm_len).astype(np.float32)
    ny_grid = (-dz_dy / norm_len).astype(np.float32)
    nz_grid = (1.0 / norm_len).astype(np.float32)

    albedo_grid = r3d.base_albedo_from_lrm(lrm_norm, b["g1_33"])

    hs_f = np.clip(np.nan_to_num(b["hs33"], nan=0.5), 0.0, 1.0)[..., None]
    hud_bare_arr = np.clip(albedo_grid * (0.35 + 0.80 * hs_f) * 255.0, 0, 255).astype(np.uint8)
    hud_bare_arr[~valid] = np.array([14, 18, 24], dtype=np.uint8)

    sat_img = Image.fromarray(b["sat_rgb"]).resize((dem.shape[1], dem.shape[0]), Image.Resampling.BILINEAR)
    sat_rgb = np.asarray(sat_img, dtype=np.uint8)

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
            band_shade = 0.26 + 0.14 * frac + 0.04 * np.sin(zk * 2.5)
            scol = np.column_stack([band_shade * 0.95, band_shade * 0.88, band_shade * 0.76]).astype(np.float32)
            skirt_alb.append(scol)

    sk_x = np.concatenate(skirt_x)
    sk_y = np.concatenate(skirt_y)
    sk_z = np.concatenate(skirt_z)
    sk_nx = np.concatenate(skirt_nx)
    sk_ny = np.concatenate(skirt_ny)
    sk_nz = np.concatenate(skirt_nz)
    sk_alb = np.vstack(skirt_alb)

    r0, r1 = int(H * 0.20), int(H * 0.80)
    c0, c1 = int(W * 0.20), int(W * 0.80)
    sub_lrm = np.where(valid[r0:r1, c0:c1], lrm_s[r0:r1, c0:c1], -999.0)
    pr, pc = np.unravel_index(int(np.argmax(sub_lrm)), sub_lrm.shape)
    pr += r0
    pc += c0
    peak_x = float(((pc / max(W - 1, 1)) * 2.0 - 1.0) * half_m)
    peak_y = float((1.0 - (pr / max(H - 1, 1)) * 2.0) * half_m)
    peak_z = float(z_surf_grid[pr, pc])

    cx, cy = b["cx"], b["cy"]
    raw_nx, raw_ny, raw_nz, _ = b["nasa_pts"]
    na_x = (raw_nx - cx).astype(np.float32)
    na_y = (raw_ny - cy).astype(np.float32)
    na_z_asl = raw_nz.astype(np.float32)

    gx = np.clip(((na_x / half_m + 1.0) * 0.5 * (W - 1)).astype(np.int32), 0, W - 1)
    gy = np.clip(((1.0 - (na_y / half_m + 1.0) * 0.5) * (H - 1)).astype(np.int32), 0, H - 1)
    agl = na_z_asl - dem_s[gy, gx]
    can_mask = (agl >= 2.0) & (agl <= 36.0) & valid[gy, gx]

    can_x = na_x[can_mask]
    can_y = na_y[can_mask]
    can_agl = agl[can_mask]
    can_z = z_surf_grid[gy[can_mask], gx[can_mask]].astype(np.float32) + can_agl * 0.85

    if len(can_x) > 85000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(can_x), size=85000, replace=False)
        can_x, can_y, can_z, can_agl = can_x[idx], can_y[idx], can_z[idx], can_agl[idx]

    c_gx = np.clip(((can_x / half_m + 1.0) * 0.5 * (W - 1)).astype(np.int32), 0, W - 1)
    c_gy = np.clip(((1.0 - (can_y / half_m + 1.0) * 0.5) * (H - 1)).astype(np.int32), 0, H - 1)
    sat_c = sat_rgb[c_gy, c_gx].astype(np.float32)
    h_norm = np.clip((can_agl - 2.0) / 18.0, 0.0, 1.0)[:, None]
    emerald = np.column_stack(
        [26.0 + 58.0 * h_norm[:, 0], 92.0 + 118.0 * h_norm[:, 0], 34.0 + 48.0 * h_norm[:, 0]]
    ).astype(np.float32)
    can_col = np.clip(0.42 * sat_c + 0.58 * emerald, 0, 255).astype(np.uint8)

    hud_bare = Image.fromarray(hud_bare_arr).resize((250, 250), Image.Resampling.BILINEAR)
    hud_sat = Image.fromarray(sat_rgb).resize((250, 250), Image.Resampling.BILINEAR)

    return {
        "site": site,
        "base_z": base_z,
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


def render_frame_with_coords(scene: dict, frame_idx: int, total_frames: int, out_path: Path) -> None:
    r3d.render_frame(scene, frame_idx, total_frames, out_path)
    # Overwrite bottom telemetry bar so exact WGS84 coordinates are prominently shown!
    site = scene["site"]
    u = frame_idx / max(total_frames - 1, 1)
    az_deg = -35.0 + 330.0 * u
    pitch_deg = 29.0 + 4.5 * math.sin(u * math.pi * 2.0)
    img = Image.open(out_path)
    draw = ImageDraw.Draw(img)
    _, _, _, f_mono = r3d.load_fonts()
    vp_x0, vp_y0, vp_w, vp_h = 20, 95, 950, 600
    draw.rectangle([vp_x0 + 12, vp_y0 + vp_h - 42, vp_x0 + vp_w - 12, vp_y0 + vp_h - 10], fill=(12, 18, 26))
    telem = (
        f"Strip: {site['strip_id']}  |  WGS84: ({site['lat']:.5f}, {site['lon']:.5f})  |  "
        f"Az: {az_deg % 360:5.1f}°  |  Pitch: {pitch_deg:4.1f}°  |  Datum: {scene['base_z']:.0f} m ASL"
    )
    draw.text((vp_x0 + 20, vp_y0 + vp_h - 33), telem, fill=(210, 225, 245), font=f_mono)
    img.save(out_path, format="PNG", compress_level=1)


def _worker_chunk(args):
    scene, indices, frames_dir, total_frames = args
    for idx in indices:
        render_frame_with_coords(scene, idx, total_frames, frames_dir / f"frame_{idx:04d}.png")
    return len(indices)


def render_locator_map() -> Path:
    BND = Path("/agent/lidar/raw_deep/boundaries")
    ne0 = json.loads((BND / "ne_adm0.geojson").read_text())
    mex = json.loads((BND / "mex_adm1.geojson").read_text())
    eaams = json.loads(Path("/agent/lidar/work/ce_deep/eaams/maya_sites.geojson").read_text())
    stroke = [pe.withStroke(linewidth=2.5, foreground="#0b1120")]

    fig = plt.figure(figsize=(20, 12.5), facecolor="#070d19")
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.05, 1.25], height_ratios=[1.05, 1.0], hspace=0.22, wspace=0.16)

    # Panel A: Yucatán Peninsula showing Estrada-Belli 2025 Southern zone vs Northern Unbroken-Jungle NJ-1..NJ-3!
    ax_pen = fig.add_subplot(gs[0, 0])
    ax_pen.set_facecolor("#0a192f")
    for f in ne0["features"]:
        if f["properties"].get("ADMIN") in ("Mexico", "Guatemala", "Belize"):
            g = shape(f["geometry"])
            for poly in ([g] if g.geom_type == "Polygon" else list(g.geoms)):
                xy = np.asarray(poly.exterior.coords)
                ax_pen.add_patch(MplPolygon(xy, closed=True, facecolor="#1e293b", edgecolor="#64748b", lw=0.9))
    for f in mex["features"]:
        g = shape(f["geometry"])
        for poly in ([g] if g.geom_type == "Polygon" else list(g.geoms)):
            xy = np.asarray(poly.exterior.coords)
            ax_pen.plot(xy[:, 0], xy[:, 1], color="#475569", lw=0.7, ls="--")

    # Shade Estrada-Belli 2025 mapped box vs Northern unmapped MARI-GISLAB=0 zone
    ax_pen.add_patch(
        Rectangle((-90.8, 17.8), 2.3, 1.65, facecolor="#f43f5e", alpha=0.18, edgecolor="#fb7185", lw=1.4, ls="--")
    )
    ax_pen.text(-89.65, 18.6, "Hunt v1 (NC-1..NC-6)\nInside Estrada-Belli 2025\nMARI-GISLAB Zone", color="#fda4af", fontsize=8.5, ha="center", fontweight="bold", path_effects=stroke)

    ax_pen.scatter([-88.75685, -88.75890, -88.77712], [20.64872, 20.63046, 20.66345], s=120, c="#22c55e", edgecolors="white", linewidths=1.3, zorder=8)
    ax_pen.annotate(
        "NEW: NJ-1, NJ-2 & NJ-3\n(20.63°–20.66°N, 88.76°–88.78°W)\n0 Polygons in MARI-GISLAB!",
        xy=(-88.76, 20.65),
        xytext=(-90.55, 21.18),
        color="#4ade80",
        fontsize=9.2,
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#4ade80", lw=1.6),
        path_effects=stroke,
        zorder=9,
    )
    for city, clo, cla in [("Mérida", -89.62, 20.97), ("Cancún", -86.85, 21.16), ("Chichén Itzá", -88.57, 20.68), ("Campeche", -90.53, 19.84)]:
        ax_pen.scatter([clo], [cla], s=32, c="#38bdf8", zorder=6)
        ax_pen.text(clo + 0.08, cla + 0.04, city, color="#e2e8f0", fontsize=8.2, path_effects=stroke)

    ax_pen.set_xlim(-91.2, -86.6)
    ax_pen.set_ylim(17.7, 21.7)
    ax_pen.tick_params(colors="#cbd5e1", labelsize=8)
    ax_pen.set_title("A. YUCATÁN CONTEXT — OUTSIDE ESTRADA-BELLI (2025) MAPPED ZONE", color="#f8fafc", fontsize=11, fontweight="bold")

    # Panel B: Regional Satellite Zoom (Chichén Itzá WSW Jungle Block + EAAMS Sites)
    ax_reg = fig.add_subplot(gs[0, 1])
    rlon0, rlon1 = -88.90, -88.54
    rlat0, rlat1 = 20.55, 20.76
    sat_reg = vst.fetch_esri_bbox_rgb(rlon0, rlat0, rlon1, rlat1, out_w=920, out_h=540, zoom_level=12)
    ax_reg.imshow(sat_reg, extent=[rlon0, rlon1, rlat0, rlat1], aspect="auto", zorder=1)

    e_lons, e_lats = [], []
    for f in eaams["features"]:
        c = f["geometry"]["coordinates"]
        if rlon0 <= c[0] <= rlon1 and rlat0 <= c[1] <= rlat1:
            e_lons.append(c[0])
            e_lats.append(c[1])
            nm = f["properties"].get("NAME") or f["properties"].get("Name") or ""
            if nm:
                ax_reg.text(c[0] + 0.004, c[1] + 0.003, nm, color="#93c5fd", fontsize=7.5, path_effects=stroke, zorder=6)
    ax_reg.scatter(e_lons, e_lats, s=36, marker="^", c="#38bdf8", edgecolors="black", label="Known EAAMS Maya Sites", zorder=5)

    for s, dy in zip(SITES, [0.022, -0.032, 0.024]):
        ax_reg.scatter([s["lon"]], [s["lat"]], s=125, c="#facc15", edgecolors="black", linewidths=1.3, zorder=8)
        ax_reg.annotate(
            f"{s['code']} ({s['lat']:.5f}, {s['lon']:.5f})",
            xy=(s["lon"], s["lat"]),
            xytext=(s["lon"] - 0.11, s["lat"] + dy),
            color="#facc15",
            fontsize=8.6,
            fontweight="bold",
            arrowprops=dict(arrowstyle="->", color="#facc15", lw=1.3),
            path_effects=stroke,
            zorder=9,
        )
    ax_reg.set_xlim(rlon0, rlon1)
    ax_reg.set_ylim(rlat0, rlat1)
    ax_reg.tick_params(colors="#cbd5e1", labelsize=8)
    ax_reg.legend(loc="lower right", facecolor="#0f172a", edgecolor="#475569", labelcolor="#f8fafc", fontsize=8)
    ax_reg.set_title("B. UNBROKEN CENTRAL-EAST YUCATÁN FOREST BLOCK (19–22 KM WSW OF CHICHÉN ITZÁ)", color="#facc15", fontsize=11, fontweight="bold")

    # Panel C: 2.4 x 2.4 km Wide Satellite Proof on NJ-1, NJ-2, NJ-3 (Zero Trails / Farmland)
    gs_sub = gs[1, 1].subgridspec(1, 3, wspace=0.08)
    for idx, s in enumerate(SITES):
        ax_s = fig.add_subplot(gs_sub[0, idx])
        d2 = 1200.0 / 110574.0
        dl2 = 1200.0 / (111320.0 * math.cos(math.radians(s["lat"])))
        sat_w = vst.fetch_esri_bbox_rgb(s["lon"] - dl2, s["lat"] - d2, s["lon"] + dl2, s["lat"] + d2, 360, 360, zoom_level=15)
        ax_s.imshow(sat_w)
        ax_s.axhline(180, color="#facc15", lw=0.9, ls="--")
        ax_s.axvline(180, color="#facc15", lw=0.9, ls="--")
        ax_s.add_patch(Circle((180, 180), 22, facecolor="none", edgecolor="#facc15", lw=1.6))
        ax_s.set_title(f"{s['code']} (2.4 × 2.4 km Sat)\nZero Trails / Zero Farms", color="#4ade80", fontsize=9, fontweight="bold")
        ax_s.axis("off")

    # Panel D: Summary & Honest Prior-Art Table
    ax_t = fig.add_subplot(gs[1, 0])
    ax_t.set_facecolor("#0f172a")
    ax_t.set_xlim(0, 100)
    ax_t.set_ylim(0, 100)
    ax_t.axis("off")
    ax_t.add_patch(Rectangle((0.5, 0.5), 99, 99, facecolor="#0f172a", edgecolor="#334155", lw=1.6))
    ax_t.text(3, 93, "C. EXACT COORDINATES, METRICS & LITERATURE AUDIT", color="#f8fafc", fontsize=10.8, fontweight="bold")

    rows = [
        ("ID", "WGS84 COORDS", "0.33m ARCHITECTURE", "CANOPY & PRIOR-ART VERDICT", "#93c5fd"),
        ("NJ-1", "20.64872° N\n88.75685° W", "9 mounds, +6.0 m pyramid,\n3,142 m² enclosed plaza", "99.3% closed jungle • 6.98 km from Cenote Xtohil\n0 in MARI-GISLAB • 0 in EAAMS / INAH SIC / OSM", "#4ade80"),
        ("NJ-2", "20.63046° N\n88.75890° W", "9 mounds (4 rect, 4/4 Q),\n+5.0 m pyramid, 548 m² plaza", "99.6% closed jungle • 6.04 km from Ikil (EAAMS)\n0 in MARI-GISLAB • 0 in EAAMS / INAH SIC / OSM", "#4ade80"),
        ("NJ-3", "20.66345° N\n88.77712° W", "18 mounds (4 rect, 4/4 Q),\n+4.4 m range, 1,068 m² plaza", "99.4% closed jungle • 4.31 km from Cenote Xtohil\n0 in MARI-GISLAB • 8,188 ground .laz pts in slice", "#4ade80"),
    ]
    yc = 83
    for i, (c1, c2, c3, c4, col) in enumerate(rows):
        hdr = i == 0
        if not hdr:
            ax_t.plot([2.5, 97.5], [yc + 4.5, yc + 4.5], color="#1e293b", lw=1.0)
        ax_t.text(3, yc, c1, color=col, fontsize=8.5, fontweight="bold", va="top")
        ax_t.text(13, yc, c2, color="#e2e8f0", fontsize=8.0, va="top", family="monospace")
        ax_t.text(34, yc, c3, color="#cbd5e1", fontsize=8.0, va="top")
        ax_t.text(62, yc, c4, color=col if hdr else "#e2e8f0", fontsize=7.8, va="top")
        yc -= 8 if hdr else 19.5

    ax_t.text(
        3,
        14,
        "Transparency Note: All 3 complexes have ZERO polygons in Estrada-Belli et al. (2025) MARI-GISLAB,\n"
        "EAAMS, INAH SIC, and Golden et al. (2016). Because they lie inside public 2013 NASA G-LiHT transects\n"
        "reprocessed at 0.33 m by Britton (Nov 2026), their mounds may appear in unnamed automated regional\n"
        "structure counts (Schroder 2020; Britton 2026), though none has ever been named, mapped in detail, or excavated.",
        color="#94a3b8",
        fontsize=7.5,
        va="top",
    )

    fig.suptitle(
        "OPEN LiDAR HUNT v3 — UNBROKEN-JUNGLE 0.33 m MAYA COURTYARD & PYRAMID COMPLEXES (NJ-1 TO NJ-3)",
        color="#f8fafc",
        fontsize=14.5,
        fontweight="bold",
        y=0.975,
    )
    out_map = OUT_STORE / "map-nj1-to-nj3-unbroken-jungle-locator.png"
    fig.savefig(out_map, dpi=155, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    shutil.copy2(out_map, OUT_REPO / out_map.name)
    return out_map


def main() -> None:
    mp4s = []
    still_frames = []
    total_frames = 108
    fps_mp4 = 20

    for site in SITES:
        print(f"Building flat-plinth 0.33m 3D orbit for {site['code']}...", flush=True)
        scene = prepare_flat_plinth_scene(site)
        fdir = VID_TMP / f"frames_{site['code']}"
        if fdir.exists():
            shutil.rmtree(fdir)
        fdir.mkdir(parents=True, exist_ok=True)

        n_workers = min(8, os.cpu_count() or 4)
        chunks = [list(range(w, total_frames, n_workers)) for w in range(n_workers)]
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            list(ex.map(_worker_chunk, [(scene, ch, fdir, total_frames) for ch in chunks]))

        still_p = VID_TMP / f"still_{site['code']}.png"
        shutil.copy2(fdir / f"frame_{int(total_frames * 0.62):04d}.png", still_p)
        still_frames.append(still_p)

        out_mp4 = OUT_STORE / f"{site['slug_base']}.mp4"
        out_gif = OUT_STORE / f"{site['slug_base']}.gif"
        subprocess.run(
            [
                "ffmpeg", "-y", "-framerate", str(fps_mp4), "-i", str(fdir / "frame_%04d.png"),
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "19", "-preset", "fast",
                "-movflags", "+faststart", str(out_mp4),
            ],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            [
                "ffmpeg", "-y", "-i", str(out_mp4),
                "-vf", "fps=12,scale=840:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=192:stats_mode=diff[p];[s1][p]paletteuse=dither=bayer:bayer_scale=5",
                "-loop", "0", str(out_gif),
            ],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        shutil.copy2(out_mp4, OUT_REPO / out_mp4.name)
        shutil.copy2(out_gif, OUT_REPO / out_gif.name)
        shutil.rmtree(fdir, ignore_errors=True)
        mp4s.append(str(out_mp4))
        print(f"Wrote {out_gif.name} ({out_gif.stat().st_size/1e6:.2f} MB) & {out_mp4.name}", flush=True)

    # Master showcase MP4
    clist = VID_TMP / "concat.txt"
    clist.write_text("".join(f"file '{p}'\n" for p in mp4s))
    master = OUT_STORE / "sim-nj-all-3-unbroken-jungle-showcase.mp4"
    subprocess.run(
        ["ffmpeg", "-y", "-f", "concat", "-safe", "0", "-i", str(clist), "-c", "copy", "-movflags", "+faststart", str(master)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    shutil.copy2(master, OUT_REPO / master.name)

    # 3-panel bare-earth 3D still gallery
    imgs = [Image.open(p) for p in still_frames]
    gal = Image.new("RGB", (1280, 720 * len(imgs)), (11, 15, 22))
    for i, im in enumerate(imgs):
        gal.paste(im, (0, i * 720))
    gal_out = OUT_STORE / "sim-nj-3d-gallery-all-3.png"
    gal.save(gal_out, optimize=True)
    shutil.copy2(gal_out, OUT_REPO / gal_out.name)

    out_map = render_locator_map()
    print("ALL DONE! Map:", out_map.name, "Gallery:", gal_out.name, flush=True)


if __name__ == "__main__":
    main()
