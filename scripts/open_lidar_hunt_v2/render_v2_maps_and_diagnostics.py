#!/usr/bin/env python3
"""Render Hunt-v1-grade cartographic locator atlas AND 5-panel raw 3D point-cloud
cross-section diagnostic sheets for Open LiDAR Hunt v2 (Motagua UAV + Acre EBA).

Outputs (written to both agent store `media/open-lidar-hunt-v2/` and repo):
  1. `v2-motagua-and-acre-globe-locator-map.png`
     4-panel cartographic atlas with 3D Orthographic Western Hemisphere Globe,
     Motagua Valley (Guatemala) regional zoom + 8 UAV footprints + UNESCO Quirigua +
     Taube/Hruby jade zone, Acre (Brazil) regional zoom + EBA LAZ swaths + 6,298 Jacobs
     geoglyphs + Peripato 2023 points, and full coordinate/verdict reference table.
  2. Six 5-panel Hunt-v1-style point-cloud cross-section sheets:
     - `sheet-qc-quirigua-acropolis-control.png`
     - `sheet-ac-acre-ring-geoglyph-control.png`
     - `sheet-am05-acre-canopy-audit.png`
     - `sheet-mg-gualan-a-terrace-complex.png`
     - `sheet-mg-gualan-b-roadside-audit.png`
     - `sheet-mg-rio-tambor-jade-hinterland.png`
"""
from __future__ import annotations

import csv
import json
import math
import shutil
import sys
import zipfile
from pathlib import Path

import laspy
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import rdata
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, Polygon as MplPolygon, Rectangle
from pyproj import Transformer
from rasterio.warp import transform_bounds
from shapely.geometry import shape

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2_scene_tools as vst

WORK = vst.WORK
STORE_OUT = vst.STORE_MEDIA
REPO_OUT = vst.REPO_MEDIA
BND = Path("/agent/lidar/raw_deep/boundaries")

stroke_dark = [pe.withStroke(linewidth=2.6, foreground="#0b1120")]
stroke_med = [pe.withStroke(linewidth=1.9, foreground="#0b1120")]


def sync_out(src_path: Path) -> None:
    dst = REPO_OUT / src_path.name
    if src_path.resolve() != dst.resolve():
        shutil.copy2(src_path, dst)


def hav_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2.0 * r * math.asin(math.sqrt(a))


def plot_ortho_globe(ax_globe, center_lon: float = -78.0, center_lat: float = 4.0) -> None:
    """Render 3D orthographic hemisphere showing Yucatán (v1), Motagua (v2), and Acre (v2)."""
    ne0 = json.loads((BND / "ne_adm0.geojson").read_text())
    ax_globe.set_facecolor("#070d19")
    ax_globe.set_xlim(-1.08, 1.08)
    ax_globe.set_ylim(-1.08, 1.08)
    ax_globe.set_aspect("equal")
    ax_globe.axis("off")

    ocean = Circle((0, 0), 1.0, facecolor="#0a192f", edgecolor="#38bdf8", lw=1.8, zorder=1)
    ax_globe.add_patch(ocean)

    lam0 = math.radians(center_lon)
    phi0 = math.radians(center_lat)

    def proj_pts(lons_deg, lats_deg):
        lam = np.radians(lons_deg)
        phi = np.radians(lats_deg)
        cos_c = np.sin(phi0) * np.sin(phi) + np.cos(phi0) * np.cos(phi) * np.cos(lam - lam0)
        x = np.cos(phi) * np.sin(lam - lam0)
        y = np.cos(phi0) * np.sin(phi) - np.sin(phi0) * np.cos(phi) * np.cos(lam - lam0)
        return x, y, cos_c

    for lat_g in range(-60, 80, 20):
        lons = np.linspace(-180, 180, 360)
        lats = np.full_like(lons, lat_g)
        x, y, c = proj_pts(lons, lats)
        x[c < 0] = np.nan
        y[c < 0] = np.nan
        ax_globe.plot(x, y, color="#1e3a5f", lw=0.5, zorder=2)

    for lon_g in range(-160, 20, 20):
        lats = np.linspace(-80, 80, 240)
        lons = np.full_like(lats, lon_g)
        x, y, c = proj_pts(lons, lats)
        x[c < 0] = np.nan
        y[c < 0] = np.nan
        ax_globe.plot(x, y, color="#1e3a5f", lw=0.5, zorder=2)

    highlight = {"Guatemala": "#f59e0b", "Brazil": "#10b981", "Mexico": "#0ea5e9"}
    for f in ne0["features"]:
        adm = f["properties"].get("ADMIN", "")
        fc = highlight.get(adm, "#293548")
        ec = "#e2e8f0" if adm in highlight else "#64748b"
        lw = 0.85 if adm in highlight else 0.4
        alpha = 0.85 if adm in highlight else 0.95
        g = shape(f["geometry"])
        polys = [g] if g.geom_type == "Polygon" else list(g.geoms)
        for poly in polys:
            coords = np.asarray(poly.exterior.coords)
            if len(coords) > 300:
                coords = coords[:: max(1, len(coords) // 300)]
            x, y, c = proj_pts(coords[:, 0], coords[:, 1])
            if np.mean(c > 0) < 0.55:
                continue
            vis = c > 0
            if np.sum(vis) >= 3:
                xy = np.column_stack([x[vis], y[vis]])
                ax_globe.add_patch(
                    MplPolygon(xy, closed=True, facecolor=fc, edgecolor=ec, lw=lw, alpha=alpha, zorder=3)
                )

    for label, lon_p, lat_p, col, dy_lbl in [
        ("Hunt v1: Yucatán\n(Campeche NC-1..6)", -89.8, 18.7, "#38bdf8", 0.11),
        ("Hunt v2A: Motagua Fault UAV\n(Guatemala: Quiriguá / Gualán)", -89.3, 15.2, "#facc15", -0.15),
        ("Hunt v2B: Acre EBA LiDAR\n(SW Amazonia, Brazil)", -67.6, -9.9, "#34d399", -0.14),
    ]:
        px, py, _ = proj_pts(np.array([lon_p]), np.array([lat_p]))
        ax_globe.scatter(px, py, s=95, c=col, edgecolors="white", linewidths=1.3, zorder=8)
        ax_globe.plot([px[0], px[0] + 0.14], [py[0], py[0] + dy_lbl], color=col, lw=1.2, zorder=7)
        ax_globe.text(
            px[0] + 0.16,
            py[0] + dy_lbl,
            label,
            color=col,
            fontsize=8.5,
            fontweight="bold",
            va="center",
            path_effects=stroke_dark,
            zorder=9,
        )

    ax_globe.set_title(
        "A. WESTERN HEMISPHERE CONTEXT — HUNT v1 vs. HUNT v2 REGIONS",
        color="#f8fafc",
        fontsize=11.5,
        fontweight="bold",
        pad=8,
    )


def render_locator_atlas() -> Path:
    print("Rendering 4-panel v2 globe & regional locator atlas...", flush=True)
    fig = plt.figure(figsize=(21, 13.5), facecolor="#070d19")
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.0, 1.32], height_ratios=[1.05, 1.0], hspace=0.22, wspace=0.16)

    # Panel A: Globe
    ax_globe = fig.add_subplot(gs[0, 0])
    plot_ortho_globe(ax_globe)

    # Panel B: Motagua Fault Corridor (Guatemala)
    ax_mot = fig.add_subplot(gs[0, 1])
    mlon0, mlon1 = -90.02, -88.70
    mlat0, mlat1 = 14.82, 15.63
    sat_mot = vst.fetch_esri_bbox_rgb(mlon0, mlat0, mlon1, mlat1, out_w=960, out_h=590, zoom_level=9)
    ax_mot.imshow(sat_mot, extent=[mlon0, mlon1, mlat0, mlat1], aspect="auto", alpha=0.86, zorder=1)

    ne0 = json.loads((BND / "ne_adm0.geojson").read_text())
    for f in ne0["features"]:
        adm = f["properties"].get("ADMIN", "")
        if adm in ("Guatemala", "Honduras", "Belize", "El Salvador"):
            g = shape(f["geometry"])
            polys = [g] if g.geom_type == "Polygon" else list(g.geoms)
            for poly in polys:
                xy = np.asarray(poly.exterior.coords)
                ax_mot.plot(xy[:, 0], xy[:, 1], color="#cbd5e1", lw=1.0, ls="--", alpha=0.7, zorder=2)

    # Plot all 8 Motagua UAV DTM footprints
    for tif in sorted((WORK / "motagua/dtm").glob("*_dtm_cog.tif")):
        with rasterio.open(tif) as ds:
            b = transform_bounds(ds.crs, "EPSG:4326", *ds.bounds)
        rect = Rectangle(
            (b[0], b[1]),
            max(b[2] - b[0], 0.020),
            max(b[3] - b[1], 0.020),
            facecolor="#38bdf8",
             alpha=0.28,
            edgecolor="#38bdf8",
            lw=1.5,
            zorder=4,
        )
        ax_mot.add_patch(rect)

    mot_pts = [
        ("QC-1: Quiriguá UNESCO Acropolis\n(15.2689°N, 89.0403°W) [CONTROL]", 15.26885, -89.04025, "#22c55e", (0.04, 0.05)),
        ("MG-Gualan A: River-Terrace Complex\n(15.1343°N, 89.3164°W) [STRONGEST LEAD]", 15.13427, -89.31641, "#facc15", (-0.18, -0.14)),
        ("MG-Gualan B: S. of Roadway\n(15.1257°N, 89.3288°W) [DOWNGRADED]", 15.12571, -89.32880, "#fb923c", (0.05, -0.11)),
        ("MG-RioTambor: Dissected Hills\n(14.9393°N, 89.8167°W) [JADE ZONE ~21km]", 14.93932, -89.81672, "#f87171", (0.05, 0.03)),
    ]
    for lbl, lat, lon, col, (dx, dy) in mot_pts:
        ax_mot.scatter([lon], [lat], s=95, c=col, edgecolors="white", linewidths=1.3, zorder=8)
        ax_mot.annotate(
            lbl,
            xy=(lon, lat),
            xytext=(lon + dx, lat + dy),
            color=col,
            fontsize=8.2,
            fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=col, lw=1.2),
            path_effects=stroke_dark,
            zorder=9,
        )

    # Known archaeological reference: Carrizal Grande / Río El Tambor Jade Zone (Taube et al. 2004)
    ax_mot.scatter([-89.96], [14.96], s=65, marker="D", c="#c084fc", edgecolors="white", zorder=7)
    ax_mot.text(
        -89.95,
        15.01,
        "Taube/Hruby 2004\nJade Workshops",
        color="#e9d5ff",
        fontsize=7.8,
        fontweight="bold",
        path_effects=stroke_med,
        zorder=8,
    )
    # Lake Izabal label
    ax_mot.text(-89.16, 15.51, "Lago de Izabal", color="#7dd3fc", fontsize=9, style="italic", path_effects=stroke_med)
    ax_mot.text(-89.55, 15.07, "Río Motagua Fault Trace →", color="#cbd5e1", fontsize=8.5, rotation=15, path_effects=stroke_med)

    ax_mot.set_xlim(mlon0, mlon1)
    ax_mot.set_ylim(mlat0, mlat1)
    ax_mot.set_xlabel("Longitude (°W)", color="#cbd5e1", fontsize=9)
    ax_mot.set_ylabel("Latitude (°N)", color="#cbd5e1", fontsize=9)
    ax_mot.tick_params(colors="#cbd5e1", labelsize=8)
    ax_mot.set_title(
        "B. MOTAGUA–POLOCHIC UAV LiDAR CORRIDOR (GUATEMALA — JAN 2024 / PUB. JUL 2026)",
        color="#facc15",
        fontsize=11.5,
        fontweight="bold",
        pad=8,
    )

    # Panel C: Acre EBA Transects & Known Geoglyph Province (Brazil)
    ax_acre = fig.add_subplot(gs[1, 1])
    alon0, alon1 = -68.35, -66.75
    alat0, alat1 = -10.50, -9.55
    sat_acre = vst.fetch_esri_bbox_rgb(alon0, alat0, alon1, alat1, out_w=960, out_h=570, zoom_level=9)
    ax_acre.imshow(sat_acre, extent=[alon0, alon1, alat0, alat1], aspect="auto", alpha=0.84, zorder=1)

    # Plot Jacobs geoglyphs + Peripato points in bbox
    j_lons, j_lats = [], []
    with open(WORK / "known_sites/amazon_geoglyphs_points.csv") as f:
        for r in csv.DictReader(f):
            la, lo = float(r["lat"]), float(r["lon"])
            if alat0 <= la <= alat1 and alon0 <= lo <= alon1:
                j_lats.append(la)
                j_lons.append(lo)
    ax_acre.scatter(j_lons, j_lats, s=10, c="#38bdf8", alpha=0.55, label="Jacobs Geoglyph Catalog (n=6,298 total)", zorder=3)

    df_p = rdata.read_rds(WORK / "known_sites/peripato/Vperipato-ade2541-78f685a/Database/Earthworks.rds")
    trees_df = df_p[df_p["Database"] == "TREES/INPE"]
    ax_acre.scatter(
        trees_df["Longitude"],
        trees_df["Latitude"],
        s=42,
        marker="^",
        c="#c084fc",
        edgecolors="black",
        linewidths=0.5,
        label="Peripato et al. 2023 Science (24 EBA Hits)",
        zorder=5,
    )

    # Plot EBA Acre LAZ transect footprints in bbox
    zf = zipfile.ZipFile(WORK / "amazon_eba/transectos_acre.zip")
    for info in zf.infolist():
        if not info.filename.endswith(".laz"):
            continue
        with zf.open(info) as fp:
            with laspy.open(fp) as rdr:
                h = rdr.header
                epsg = h.parse_crs().to_epsg()
                tr_ll = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
                lo0, la0 = tr_ll.transform(h.mins[0], h.mins[1])
                lo1, la1 = tr_ll.transform(h.maxs[0], h.maxs[1])
        if lo1 >= alon0 and lo0 <= alon1 and la1 >= alat0 and la0 <= alat1:
            ax_acre.add_patch(
                Rectangle(
                    (lo0, la0),
                    max(lo1 - lo0, 0.02),
                    max(la1 - la0, 0.02),
                    facecolor="#f43f5e",
                    alpha=0.26,
                    edgecolor="#fb7185",
                    lw=1.5,
                    zorder=4,
                )
            )
            ax_acre.text(lo0, la1 + 0.012, info.filename.replace(".laz", ""), color="#fda4af", fontsize=7.5, path_effects=stroke_med, zorder=6)

    acre_pts = [
        ("AC-1 Control: Circular Geoglyph + Mound Village\n(-9.8545°S, -67.2714°W | acrc33 + mouv65)", -9.85446, -67.27142, "#22c55e", (-0.58, 0.14)),
        ("AM-05 Audit: Forest/Pasture Canopy Boundary\n(-9.8970°S, -67.3030°W | Class-2 LAZ Artifact!)", -9.89700, -67.30296, "#f97316", (-0.62, -0.18)),
        ("Rio Branco (State Capital)", -9.975, -67.810, "#f8fafc", (-0.28, 0.09)),
    ]
    for lbl, lat, lon, col, (dx, dy) in acre_pts:
        ax_acre.scatter([lon], [lat], s=95, c=col, edgecolors="white", linewidths=1.3, zorder=8)
        ax_acre.annotate(
            lbl,
            xy=(lon, lat),
            xytext=(lon + dx, lat + dy),
            color=col,
            fontsize=8.2,
            fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=col, lw=1.2),
            path_effects=stroke_dark,
            zorder=9,
        )

    ax_acre.set_xlim(alon0, alon1)
    ax_acre.set_ylim(alat0, alat1)
    ax_acre.set_xlabel("Longitude (°W)", color="#cbd5e1", fontsize=9)
    ax_acre.set_ylabel("Latitude (°S)", color="#cbd5e1", fontsize=9)
    ax_acre.tick_params(colors="#cbd5e1", labelsize=8)
    leg = ax_acre.legend(loc="lower right", facecolor="#0f172a", edgecolor="#475569", fontsize=7.8, labelcolor="#e2e8f0")
    leg.set_zorder(10)
    ax_acre.set_title(
        "C. ACRE GEOGLYPH PROVINCE & INPE EBA LiDAR TRANSECTS (SW AMAZONIA, BRAZIL)",
        color="#34d399",
        fontsize=11.5,
        fontweight="bold",
        pad=8,
    )

    # Panel D: Forensic Summary & Exact Coordinates Table
    ax_tbl = fig.add_subplot(gs[1, 0])
    ax_tbl.set_facecolor("#0f172a")
    ax_tbl.set_xlim(0, 100)
    ax_tbl.set_ylim(0, 100)
    ax_tbl.axis("off")
    ax_tbl.add_patch(Rectangle((0.5, 0.5), 99, 99, facecolor="#0f172a", edgecolor="#334155", lw=1.6))

    ax_tbl.text(
        3.0,
        94.0,
        "D. EXACT COORDINATES & FORENSIC AUDIT (WHY PASS 1 FELT WEAK)",
        color="#f8fafc",
        fontsize=10.8,
        fontweight="bold",
    )

    rows = [
        ("SITE / CODE", "WGS84 COORDS", "RELIEF / SCALE", "FORENSIC AUDIT VERDICT", "#93c5fd"),
        (
            "QC-1 Quiriguá\n(Museo UAV)",
            "15.26885° N\n89.04025° W",
            "0.25 m DTM\n+11 m Acropolis",
            "KNOWN UNESCO CONTROL: Pass 1 cropped modern\nbanana drainage canals S of park! Fixed here.",
            "#4ade80",
        ),
        (
            "AC-1 Geoglyph\n(NP_T-0414)",
            "9.85446° S\n67.27142° W",
            "1.0 m PMF DTM\n267 m Ring Ditch",
            "KNOWN GEOGLYPH CONTROL (acrc33 + mouv65):\n2.2 m ditch + radial roads cleanly resolved.",
            "#4ade80",
        ),
        (
            "MG-Gualan A\n(Gualán UAV)",
            "15.13427° N\n89.31641° W",
            "0.50 m DTM\n+2.5–6 m Terrace",
            "STRONGEST UNCATALOGUED LEAD: Rectilinear\nriver-terrace platform surviving 3D tree strip.",
            "#facc15",
        ),
        (
            "AM-05 Acre\n(NP_T-0414)",
            "9.89700° S\n67.30296° W",
            "1.0 m PMF DTM\n<0.5 m Bare Earth",
            "ARTIFACT EXPOSED: EBA LAZ Class 2 mixed 22 m\nforest w/ 0 m pasture! Flat after PMF filter.",
            "#fb923c",
        ),
        (
            "MG-Gualan B\n(Gualán UAV)",
            "15.12571° N\n89.32880° W",
            "0.50 m DTM\n+1.8 m Platform",
            "DOWNGRADED: +16 m 'red blob' was leaked tree\ncrown; low platform sits by modern road.",
            "#fb923c",
        ),
        (
            "MG-RioTambor\n(Río Tambor)",
            "14.93932° N\n89.81672° W",
            "0.50 m DTM\nDissected Spurs",
            "AMBIGUOUS / NATURAL: ~21 km from Taube jade\nsites, but 3D stripping shows natural spurs.",
            "#f87171",
        ),
    ]

    y_cursor = 86.5
    for idx, (c1, c2, c3, c4, col) in enumerate(rows):
        is_hdr = idx == 0
        fs = 8.0 if is_hdr else 7.6
        fw = "bold" if is_hdr else "normal"
        if not is_hdr:
            ax_tbl.plot([2.5, 97.5], [y_cursor + 5.2, y_cursor + 5.2], color="#1e293b", lw=1.0)
        ax_tbl.text(3.0, y_cursor, c1, color=col, fontsize=fs, fontweight="bold", va="top")
        ax_tbl.text(23.5, y_cursor, c2, color="#e2e8f0", fontsize=fs, fontweight=fw, va="top", family="monospace")
        ax_tbl.text(42.5, y_cursor, c3, color="#cbd5e1", fontsize=fs, fontweight=fw, va="top")
        ax_tbl.text(61.5, y_cursor, c4, color=col if is_hdr else "#e2e8f0", fontsize=fs, fontweight=fw, va="top")
        y_cursor -= 7.2 if is_hdr else 12.8

    fig.suptitle(
        "OPEN LiDAR HUNT v2 — GLOBAL & REGIONAL LOCATOR ATLAS, FOOTPRINTS & FORENSIC AUDIT",
        color="#f8fafc",
        fontsize=15,
        fontweight="bold",
        y=0.975,
    )

    out_path = STORE_OUT / "v2-motagua-and-acre-globe-locator-map.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    sync_out(out_path)
    print(f"Wrote {out_path} ({out_path.stat().st_size / 1e6:.2f} MB)", flush=True)
    return out_path


def render_diagnostic_sheet(
    b: dict,
    out_name: str,
    title: str,
    subtitle: str,
    verdict_badge: str,
    badge_color: str,
    transect_p0: tuple[float, float],
    transect_p1: tuple[float, float],
    lrm_lim: float = 3.0,
    feature_note: str = "",
) -> Path:
    """Render Hunt-v1-style 5-panel technical proof sheet (4 top maps + raw 3D laser profile)."""
    dem = b["dem"]
    lrm = b["lrm"]
    g1 = b["g1_rgb"]
    chm = b["chm"]
    sat = b["sat_rgb"]
    half_m = b["half_m"]
    lat, lon = b["lat"], b["lon"]

    fig = plt.figure(figsize=(18.5, 10.8), facecolor="#0b1120")
    gs = GridSpec(2, 4, figure=fig, height_ratios=[1.12, 0.95], hspace=0.26, wspace=0.15)

    ext = [-half_m, half_m, -half_m, half_m]
    x0, y0 = transect_p0
    x1, y1 = transect_p1

    # Panel 1: Bare-Earth Hillshade + LRM blend
    ax1 = fig.add_subplot(gs[0, 0])
    hs_f = np.nan_to_num(b["hs"], nan=0.5)
    cmap_lrm = plt.get_cmap("RdYlBu_r")
    lrm_n = np.clip((np.nan_to_num(lrm, nan=0.0) + lrm_lim) / (2.0 * lrm_lim), 0.0, 1.0)
    lrm_rgb = cmap_lrm(lrm_n)[..., :3]
    blend = np.clip(lrm_rgb * (0.42 + 0.75 * hs_f[..., None]), 0.0, 1.0)
    blend[~np.isfinite(dem)] = [0.06, 0.08, 0.12]
    ax1.imshow(blend, extent=ext)
    ax1.plot([x0, x1], [y0, y1], color="#00f5ff", ls="--", lw=2.2, zorder=6)
    ax1.scatter([x0, x1], [y0, y1], c="#00f5ff", s=38, edgecolors="black", zorder=7)
    ax1.text(x0, y0 + 10, "A", color="#00f5ff", fontsize=10, fontweight="bold", path_effects=stroke_dark)
    ax1.text(x1, y1 + 10, "A'", color="#00f5ff", fontsize=10, fontweight="bold", path_effects=stroke_dark)
    ax1.set_title(f"1. Bare-Earth Hillshade + LRM (±{lrm_lim:.1f} m)", color="#facc15", fontsize=10, fontweight="bold")

    # Panel 2: G1 Architectural Slope + Multi-Hillshade
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.imshow(g1, extent=ext)
    ax2.plot([x0, x1], [y0, y1], color="#00f5ff", ls="--", lw=1.8, alpha=0.85)
    ax2.set_title(f"2. {b['res_m']:.2f} m Architectural Slope Relief", color="#e2e8f0", fontsize=10, fontweight="bold")

    # Panel 3: Canopy Height Model (CHM)
    ax3 = fig.add_subplot(gs[0, 2])
    im3 = ax3.imshow(chm, extent=ext, cmap="YlGn", vmin=0, vmax=30)
    ax3.plot([x0, x1], [y0, y1], color="#00f5ff", ls="--", lw=1.8, alpha=0.85)
    ax3.set_title(f"3. 3D Laser Canopy Height (p90={np.percentile(chm, 90):.1f} m)", color="#4ade80", fontsize=10, fontweight="bold")

    # Panel 4: Esri World Imagery Satellite
    ax4 = fig.add_subplot(gs[0, 3])
    ax4.imshow(sat, extent=ext)
    ax4.plot([x0, x1], [y0, y1], color="#00f5ff", ls="--", lw=1.8, alpha=0.9)
    ax4.set_title("4. Esri Optical Satellite (XYZ Tiles)", color="#38bdf8", fontsize=10, fontweight="bold")

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim(-half_m, half_m)
        ax.set_ylim(-half_m, half_m)
        ax.tick_params(colors="#94a3b8", labelsize=7.5)
        ax.set_xlabel("Meters E–W", color="#94a3b8", fontsize=8)
    ax1.set_ylabel("Meters N–S", color="#94a3b8", fontsize=8)

    # Panel 5: Raw 3D LiDAR Point-Cloud Cross-Section along A -> A' (10 m swath width)
    ax5 = fig.add_subplot(gs[1, :])
    ax5.set_facecolor("#0f172a")
    px = b["pts_x"] - b["cx"]
    py = b["pts_y"] - b["cy"]
    pz = b["pts_z"]
    is_gnd = b["pts_is_gnd"]

    dx_t, dy_t = x1 - x0, y1 - y0
    seg_len = math.hypot(dx_t, dy_t)
    ux, uy = dx_t / seg_len, dy_t / seg_len
    s_along = (px - x0) * ux + (py - y0) * uy
    d_perp = np.abs(-(px - x0) * uy + (py - y0) * ux)

    swath_w = 7.0
    m_sw = (s_along >= 0.0) & (s_along <= seg_len) & (d_perp <= swath_w)
    s_pts = s_along[m_sw]
    z_pts = pz[m_sw]
    g_pts = is_gnd[m_sw]

    # Sample cleaned bare-earth DTM curve along transect
    n_prof = 320
    s_grid = np.linspace(0.0, seg_len, n_prof)
    wx_g = x0 + ux * s_grid
    wy_g = y0 + uy * s_grid
    H, W = dem.shape
    c_idx = np.clip(((wx_g + half_m) / (2.0 * half_m) * (W - 1)).astype(np.int32), 0, W - 1)
    r_idx = np.clip(((half_m - wy_g) / (2.0 * half_m) * (H - 1)).astype(np.int32), 0, H - 1)
    z_dtm_prof = dem[r_idx, c_idx]

    can_m = ~g_pts
    if np.sum(can_m) > 35000:
        rng = np.random.default_rng(7)
        sel = rng.choice(np.where(can_m)[0], size=35000, replace=False)
        ax5.scatter(s_pts[sel], z_pts[sel], s=2.2, c="#22c55e", alpha=0.35, label="Raw 3D Canopy Returns (.laz / EPT)")
    else:
        ax5.scatter(s_pts[can_m], z_pts[can_m], s=2.8, c="#22c55e", alpha=0.40, label="Raw 3D Canopy Returns (.laz / EPT)")

    ax5.scatter(s_pts[g_pts], z_pts[g_pts], s=5.0, c="#fb923c", alpha=0.65, label="Sub-Canopy Ground Returns")
    ax5.plot(s_grid, z_dtm_prof, color="#facc15", lw=2.6, label="Cleaned Bare-Earth DTM Surface", zorder=6)

    # If Motagua had leaked trees clamped, plot dashed line of uncleaned DTM for complete transparency
    if "dem_raw_delivered" in b:
        z_raw_prof = b["dem_raw_delivered"][r_idx, c_idx]
        if np.nanmax(np.abs(z_raw_prof - z_dtm_prof)) > 1.2:
            ax5.plot(
                s_grid,
                z_raw_prof,
                color="#f43f5e",
                ls=":",
                lw=1.7,
                alpha=0.85,
                label="Unfiltered Delivered DTM (Leaked Tree Crowns Removed)",
                zorder=5,
            )

    z_lo = float(np.nanmin(z_dtm_prof)) - 2.5
    z_hi = max(float(np.percentile(z_pts, 99.2)) if len(z_pts) else z_lo + 15.0, float(np.nanmax(z_dtm_prof)) + 8.0)
    ax5.set_xlim(0, seg_len)
    ax5.set_ylim(z_lo, z_hi)
    ax5.grid(True, color="#1e293b", ls="--", lw=0.8)
    ax5.set_xlabel(f"Distance Along Transect A → A' (meters, ±{swath_w:.0f} m laser swath)", color="#cbd5e1", fontsize=9.5)
    ax5.set_ylabel("Elevation (m, WGS84 / UTM)", color="#cbd5e1", fontsize=9.5)
    ax5.tick_params(colors="#cbd5e1", labelsize=8.5)
    ax5.set_title(
        f"5. RAW 3D LiDAR POINT-CLOUD CROSS-SECTION (A → A', {len(s_pts):,} LASER RETURNS IN SWATH)",
        color="#facc15",
        fontsize=10.8,
        fontweight="bold",
        pad=6,
    )
    if feature_note:
        ax5.text(
            0.015,
            0.91,
            feature_note,
            transform=ax5.transAxes,
            color="#fef08a",
            fontsize=9.2,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#1e293b", edgecolor="#facc15", lw=1.1),
            zorder=9,
        )
    ax5.legend(loc="upper right", facecolor="#0b1120", edgecolor="#475569", fontsize=8.5, labelcolor="#f8fafc")

    # Header + Verdict Badge
    fig.suptitle(title, color="#f8fafc", fontsize=14, fontweight="bold", x=0.04, y=0.975, ha="left")
    fig.text(
        0.04,
        0.942,
        f"{subtitle}   |   Exact WGS84: {lat:.5f}°, {lon:.5f}°   |   Google Maps: https://maps.google.com/?q={lat:.5f},{lon:.5f}",
        color="#93c5fd",
        fontsize=9.5,
    )
    fig.text(
        0.96,
        0.962,
        verdict_badge,
        color="#0b1120",
        fontsize=9.8,
        fontweight="bold",
        ha="right",
        va="center",
        bbox=dict(boxstyle="round,pad=0.42", facecolor=badge_color, edgecolor="white", lw=1.2),
    )

    out_path = STORE_OUT / out_name
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    sync_out(out_path)
    print(f"Wrote {out_path} ({out_path.stat().st_size / 1e6:.2f} MB)", flush=True)
    return out_path


def main() -> None:
    render_locator_atlas()

    laz414 = WORK / "amazon_eba/laz/NP_T-0414.laz"

    # 1. Quirigua Acropolis Control
    b_qc = vst.fetch_motagua_bundle("QC-Acropolis", "quirigua_museo_dtm_cog.tif", 7, 15.26885, -89.04025, half_m=190.0)
    render_diagnostic_sheet(
        b_qc,
        "sheet-qc-quirigua-acropolis-control.png",
        "POSITIVE CONTROL QC-1: Quiriguá UNESCO Maya Acropolis & Ballcourt Plaza (Izabal, Guatemala)",
        "Missouri S&T UAV LiDAR (OTDS.072026.32616.7 • 166 pts/m² • 0.25 m DTM)",
        "KNOWN UNESCO WORLD HERITAGE CONTROL",
        "#4ade80",
        transect_p0=(-120.0, 125.0),
        transect_p1=(35.0, -135.0),
        lrm_lim=5.0,
        feature_note="Great Plaza (~67.5 m) rising +11 m onto the masonry Southern Acropolis temples under 30–38 m rainforest canopy",
    )

    # 2. Acre Geoglyph Control (acrc33 + mouv65)
    b_ac = vst.fetch_eba_acre_bundle("AC-Geoglyph", laz414, -9.85446, -67.27142, half_m=220.0)
    render_diagnostic_sheet(
        b_ac,
        "sheet-ac-acre-ring-geoglyph-control.png",
        "POSITIVE CONTROL AC-1: 267 m Circular Ditch Geoglyph & Mounded Ring Village (Acre, Brazil)",
        "INPE EBA Transect NP_T-0414 • 1.0 m Circular-Disk PMF Bare-Earth DTM (Jacobs acrc33 / mouv65)",
        "KNOWN AMAZON GEOGLYPH CONTROL",
        "#4ade80",
        transect_p0=(-185.0, 20.0),
        transect_p1=(185.0, -25.0),
        lrm_lim=1.8,
        feature_note="Classic -2.2 m excavated circular geoglyph ditch + outer embankment & radiating causeways on flat 143.5 m plateau",
    )

    # 3. AM-05 Forensic Audit Sheet
    b_am05 = vst.fetch_eba_acre_bundle("AM-05-Audit", laz414, -9.89700, -67.30296, half_m=200.0)
    render_diagnostic_sheet(
        b_am05,
        "sheet-am05-acre-canopy-audit.png",
        "FORENSIC AUDIT AM-05: Why Pass-1 Looked Rectangular — LAZ Class-2 Canopy Boundary Artifact (Acre, Brazil)",
        "INPE EBA Transect NP_T-0414 • 2.4 km SW of Known Quadrangle Geoglyph",
        "AUDIT: CANOPY-EDGE ARTIFACT (REJECTED)",
        "#fb923c",
        transect_p0=(-150.0, 0.0),
        transect_p1=(150.0, 0.0),
        lrm_lim=1.5,
        feature_note="West (-40..+20 m) is cleared pasture (0 m canopy); East (+40..+150 m) has 22 m trees mislabeled Class 2 in raw LAZ! True ground is smooth.",
    )

    # 4. MG-Gualan A (Strongest Motagua Candidate)
    b_ga = vst.fetch_motagua_bundle("MG-Gualan-A", "gualan_dtm_cog.tif", 8, 15.13427, -89.31641, half_m=165.0)
    render_diagnostic_sheet(
        b_ga,
        "sheet-mg-gualan-a-terrace-complex.png",
        "CANDIDATE MG-Gualan A: Elevated River-Terrace Platform & Mound Group (Gualán, Zacapa, Guatemala)",
        "Missouri S&T Fault UAV LiDAR (OTDS.072026.32616.8 • 119 pts/m² • 0.50 m DTM) • 33 km WSW of Quiriguá",
        "STRONGEST UNCATALOGUED MOTAGUA LEAD",
        "#facc15",
        transect_p0=(-125.0, 35.0),
        transect_p1=(115.0, -25.0),
        lrm_lim=3.5,
        feature_note="Solid +3–6 m rectilinear terrace platform overlooking active stream channel; survives 3D sub-canopy tree removal",
    )

    # 5. MG-Gualan B (Roadside Audit)
    b_gb = vst.fetch_motagua_bundle("MG-Gualan-B", "gualan_dtm_cog.tif", 8, 15.12571, -89.32880, half_m=165.0)
    render_diagnostic_sheet(
        b_gb,
        "sheet-mg-gualan-b-roadside-audit.png",
        "CANDIDATE AUDIT MG-Gualan B: Low Terraced Plots South of Modern Roadway (Gualán, Guatemala)",
        "Missouri S&T Fault UAV LiDAR (OTDS.072026.32616.8 • 0.50 m DTM)",
        "DOWNGRADED: ROADSIDE + TREE LEAKAGE",
        "#fb923c",
        transect_p0=(-130.0, 65.0),
        transect_p1=(130.0, -45.0),
        lrm_lim=3.0,
        feature_note="Red dotted line shows +16 m leaked tree crown in delivered DTM; true bare-earth shows +1.5 m low roadside terrace",
    )

    # 6. MG-RioTambor (Jade Hinterland Dissected Spurs)
    b_rt = vst.fetch_motagua_bundle("MG-RioTambor", "rio_tambor_dtm_cog.tif", 9, 14.93932, -89.81672, half_m=165.0)
    render_diagnostic_sheet(
        b_rt,
        "sheet-mg-rio-tambor-jade-hinterland.png",
        "CANDIDATE AUDIT MG-RioTambor: Forested Ridge Knolls ~21 km E of Taube/Hruby Jade Workshops (Zacapa, Guatemala)",
        "Missouri S&T Fault UAV LiDAR (OTDS.072026.32616.9 • 165 pts/m² • 0.50 m DTM)",
        "AMBIGUOUS: ERODED SPURS / NEEDS FIELD CHECK",
        "#f87171",
        transect_p0=(-135.0, 95.0),
        transect_p1=(125.0, -95.0),
        lrm_lim=3.5,
        feature_note="Heavily dissected tectonic ridge spurs under 17 m dry-forest canopy; lacks crisp courtyard plaza geometry of NC-1..NC-6",
    )


if __name__ == "__main__":
    main()
