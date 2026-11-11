#!/usr/bin/env python3
"""Generate a large-format, 4-panel cartographic atlas:
/cursor/stores/bc-adff6ec3-2baf-40bb-9bab-102a5f567f18/media/maya-discoveries-mexico-and-globe-map.png

Panels:
  1. Top-Left: Full Mexico & Central America Map + 3D Orthographic Globe Inset
     (showing Mexico City, Guadalajara, Monterrey, Veracruz, Oaxaca, Mérida, Cancún, etc.)
  2. Top-Right: Full Yucatán Peninsula Political, Modern Cities, Highways & 456 G-LiHT Strips Map
  3. Bottom-Left: Detailed Campeche & Quintana Roo Heartland Zoom (Esri Satellite Canopy +
     Highways Mex 186 / 261 / 269 / 307 + Modern Towns + Famous Maya Ruins + NC-1..NC-5 & C/D/E)
  4. Bottom-Right: Exact Distance & Location Reference Table (Nearest Modern Cities, Highways & Ruins)
"""
from __future__ import annotations

import io
import json
import math
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon as MplPolygon, Rectangle, Circle
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image
from shapely.geometry import shape

MEDIA = Path("/cursor/stores/bc-adff6ec3-2baf-40bb-9bab-102a5f567f18/media")
BND = Path("/agent/lidar/raw_deep/boundaries")
WORK = Path("/agent/lidar/work/carbon_raw_hunt")

mex = json.loads((BND / "mex_adm1.geojson").read_text())
ne0 = json.loads((BND / "ne_adm0.geojson").read_text())
eaams = json.loads(Path("/agent/lidar/work/ce_deep/eaams/maya_sites.geojson").read_text())
footprints = json.loads((WORK / "ot_456_footprints.json").read_text())

stroke = [pe.withStroke(linewidth=2.6, foreground="#0b1120")]
stroke_med = [pe.withStroke(linewidth=2.0, foreground="#0b1120")]
stroke_light = [pe.withStroke(linewidth=1.5, foreground="#0b1120")]


def hav_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2.0 * r * math.asin(math.sqrt(a))


def plot_geom(ax, geom_dict, fc="#1e293b", ec="#64748b", lw=0.8, ls="-", alpha=1.0, zorder=2):
    g = shape(geom_dict)
    polys = [g] if g.geom_type == "Polygon" else list(g.geoms)
    for poly in polys:
        xy = np.asarray(poly.exterior.coords)
        ax.add_patch(
            MplPolygon(xy, closed=True, facecolor=fc, edgecolor=ec, lw=lw, ls=ls, alpha=alpha, zorder=zorder)
        )


def plot_ortho_globe(ax_globe, center_lon=-90.0, center_lat=20.0):
    """Render a clean 3D orthographic globe inset centered on Mexico / Yucatán."""
    ax_globe.set_facecolor("none")
    ax_globe.set_xlim(-1.08, 1.08)
    ax_globe.set_ylim(-1.08, 1.08)
    ax_globe.set_aspect("equal")
    ax_globe.axis("off")

    # Ocean disk
    ocean = Circle((0, 0), 1.0, facecolor="#0a192f", edgecolor="#38bdf8", lw=1.6, zorder=1)
    ax_globe.add_patch(ocean)

    # Graticule lines
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
        ax_globe.plot(x, y, color="#1e3a5f", lw=0.45, zorder=2)

    for lon_g in range(-160, 20, 20):
        lats = np.linspace(-80, 80, 240)
        lons = np.full_like(lats, lon_g)
        x, y, c = proj_pts(lons, lats)
        x[c < 0] = np.nan
        y[c < 0] = np.nan
        ax_globe.plot(x, y, color="#1e3a5f", lw=0.45, zorder=2)

    # Plot continents from ne0
    for f in ne0["features"]:
        adm = f["properties"].get("ADMIN", "")
        is_mex = adm == "Mexico"
        fc = "#334155" if not is_mex else "#0f766e"
        ec = "#64748b" if not is_mex else "#2dd4bf"
        lw = 0.4 if not is_mex else 0.8
        g = shape(f["geometry"])
        polys = [g] if g.geom_type == "Polygon" else list(g.geoms)
        for poly in polys:
            coords = np.asarray(poly.exterior.coords)
            if len(coords) > 250:
                coords = coords[:: max(1, len(coords) // 250)]
            x, y, c = proj_pts(coords[:, 0], coords[:, 1])
            if np.mean(c > 0) < 0.6:
                continue
            vis = c > 0
            if np.sum(vis) >= 3:
                xy = np.column_stack([x[vis], y[vis]])
                ax_globe.add_patch(MplPolygon(xy, closed=True, facecolor=fc, edgecolor=ec, lw=lw, zorder=3))

    # Highlight Yucatán box on globe
    bx_lons = [-91.2, -87.0, -87.0, -91.2, -91.2]
    bx_lats = [17.8, 17.8, 21.5, 21.5, 17.8]
    gx, gy, _ = proj_pts(bx_lons, bx_lats)
    ax_globe.plot(gx, gy, color="#facc15", lw=1.8, zorder=6)
    ax_globe.scatter([gx[0]], [gy[0]], s=30, c="#ef4444", edgecolors="white", linewidths=0.8, zorder=7)
    ax_globe.set_title("GLOBAL POSITION", color="#93c5fd", fontsize=8.5, fontweight="bold", pad=2)


# Our 6 Unmapped Complexes + 2 Known Landscapes
DISCOVERIES = [
    {
        "code": "NC-1",
        "name": "NC-1: Twin-Pyramid Acropolis",
        "lat": 18.28691,
        "lon": -90.01315,
        "status": "NEW / UNMAPPED",
        "color": "#ef4444",
        "state_mun": "Calakmul Mun., Campeche",
        "modern_ref": "28 km S of Conhuás (Mex 186); 42 km WSW of Xpujil; 72 km SE of Escárcega",
        "ruin_ref": "29 km NW of Calakmul; 12.1 km NNE of Esperanza (EAAMS)",
    },
    {
        "code": "NC-2",
        "name": "NC-2: 8.5 m Flat-Topped Pyramid",
        "lat": 18.42650,
        "lon": -90.19837,
        "status": "NEW / UNMAPPED (#1 Score)",
        "color": "#f97316",
        "state_mun": "Calakmul / Escárcega, Campeche",
        "modern_ref": "19 km S of Constitución / Silvituc (Mex 186); 57 km ESE of Escárcega",
        "ruin_ref": "30 km WSW of Balamkú; 8.0 km WNW of Tuxpeña (EAAMS)",
    },
    {
        "code": "NC-3",
        "name": "NC-3: Hilltop Temple-Platform",
        "lat": 19.15266,
        "lon": -89.16487,
        "status": "NEW / UNMAPPED",
        "color": "#eab308",
        "state_mun": "Hopelchén / Bacalar Borderlands",
        "modern_ref": "19 km W of Rovirosa / Río Verde; 26 km NNW of Miguel Alemán; 71 km NNE of Xpujil",
        "ruin_ref": "31 km NE of Valeriana; 12.1 km NNE of Xpaicheil (EAAMS); >60 km from INAH zones",
    },
    {
        "code": "NC-4",
        "name": "NC-4: Ultra-Remote Courtyard",
        "lat": 18.76349,
        "lon": -90.68341,
        "status": "NEW / UNMAPPED (Most Isolated)",
        "color": "#ec4899",
        "state_mun": "Escárcega / Champotón, Campeche",
        "modern_ref": "22 km NNE of Escárcega; 15 km E of Mex 261; 42 km S of Champotón",
        "ruin_ref": "19.6 km WNW of Lake Pan Caliente (EAAMS); 96 km SSW of Edzná",
    },
    {
        "code": "NC-5",
        "name": "NC-5: L-Shaped Courtyard Complex",
        "lat": 18.50660,
        "lon": -90.24067,
        "status": "NEW / UNMAPPED",
        "color": "#a855f7",
        "state_mun": "Escárcega / Calakmul, Campeche",
        "modern_ref": "9 km S of Constitución / Silvituc (Mex 186); 51 km E of Escárcega",
        "ruin_ref": "31 km W of Balamkú; 10.7 km WNW of Chen Laguna (EAAMS)",
    },
    {
        "code": "NC-6",
        "name": "NC-6: 650 m Borderland Acropolis Complex",
        "lat": 18.97127,
        "lon": -89.13886,
        "status": "NEW / UNMAPPED (Borderland Corridor)",
        "color": "#22d3ee",
        "state_mun": "Bacalar Mun., Quintana Roo (830 m E of Campeche)",
        "modern_ref": "8 km NW of Ejido Miguel Alemán; 58 km NNE of Xpujil; 72 km NW of Bacalar",
        "ruin_ref": "20.3 km ENE of Valeriana; 8.6 km NW of Sacbecan (EAAMS); 52 km NW of Kinichná",
    },
]

# Major Modern Cities & Towns
MODERN_CITIES_MEX = [
    ("Mexico City (CDMX)", 19.4326, -99.1332, True),
    ("Guadalajara", 20.6597, -103.3496, True),
    ("Monterrey", 25.6866, -100.3161, True),
    ("Puebla", 19.0414, -98.2063, False),
    ("Veracruz", 19.1738, -96.1342, False),
    ("Oaxaca", 17.0732, -96.7266, False),
    ("Acapulco", 16.8531, -99.8237, False),
    ("Villahermosa", 17.9895, -92.9475, False),
    ("Tuxtla Gutiérrez", 16.7516, -93.1030, False),
    ("Mérida", 20.9674, -89.5926, True),
    ("Cancún", 21.1619, -86.8515, True),
    ("Campeche City", 19.8301, -90.5349, False),
    ("Chetumal", 18.5001, -88.2961, False),
    ("Guatemala City", 14.6349, -90.5069, False),
    ("Belize City", 17.5046, -88.1962, False),
    ("Havana", 23.1136, -82.3666, False),
]

PENINSULA_CITIES = [
    ("Mérida", 20.9674, -89.5926, "major"),
    ("Cancún", 21.1619, -86.8515, "major"),
    ("Playa del Carmen", 20.6296, -87.0739, "med"),
    ("Tulum (Town)", 20.2114, -87.4654, "med"),
    ("Valladolid", 20.6896, -88.2022, "med"),
    ("Tizimín", 21.1425, -88.1514, "med"),
    ("Tekax", 20.2025, -89.2858, "small"),
    ("Ticul", 20.3975, -89.5342, "small"),
    ("Oxkutzcab", 20.3031, -89.4186, "small"),
    ("Campeche City", 19.8301, -90.5349, "major"),
    ("Ciudad del Carmen", 18.6514, -91.8250, "major"),
    ("Champotón", 19.3547, -90.7236, "med"),
    ("Escárcega", 18.6089, -90.7453, "major"),
    ("Hopelchén", 19.7456, -89.8453, "med"),
    ("Dzibalchén", 19.4614, -89.7294, "small"),
    ("Xpujil", 18.5061, -89.3975, "major"),
    ("Constitución (Silvituc)", 18.5919, -90.1381, "small"),
    ("Conhuás", 18.5403, -89.9183, "small"),
    ("Chetumal", 18.5001, -88.2961, "major"),
    ("Bacalar", 18.6783, -88.3922, "med"),
    ("Felipe Carrillo Puerto", 19.5775, -88.0453, "med"),
    ("José María Morelos", 19.7486, -88.7086, "small"),
    ("Nicolás Bravo", 18.4581, -88.9286, "small"),
    ("Ejido Caobas", 18.4450, -89.1000, "small"),
    ("Ejido Miguel Alemán", 18.9219, -89.0825, "small"),
    ("Palenque (Town)", 17.5092, -91.9822, "med"),
    ("Tenosique", 17.4756, -91.4225, "small"),
    ("Orange Walk", 18.0812, -88.5633, "small"),
]

FAMOUS_MAYA_RUINS = [
    ("Chichén Itzá", 20.6843, -88.5678),
    ("Uxmal", 20.3594, -89.7714),
    ("Calakmul", 18.1054, -89.8100),
    ("Cobá", 20.4890, -87.7360),
    ("Tulum Ruins", 20.2147, -87.4295),
    ("Edzná", 19.5967, -90.2292),
    ("Santa Rosa Xtampak", 19.7713, -89.5981),
    ("Becán", 18.5167, -89.4667),
    ("Balamkú", 18.5580, -89.9450),
    ("Nadzca'an", 18.6061, -89.9918),
    ("Chactún", 18.5520, -89.3540),
    ("Ocomtún (2023)", 18.8381, -89.8740),
    ("Valeriana (2024)", 18.9025, -89.3175),
    ("Kohunlich", 18.4192, -88.7913),
    ("Dzibanché", 18.6394, -88.7591),
    ("Ichkabal", 18.6857, -88.6677),
    ("Chakanbakán", 18.4520, -89.0895),
    ("El Resbalón", 18.7708, -88.7940),
    ("Tikal", 17.2220, -89.6237),
    ("El Mirador", 17.7550, -89.9200),
]

# Approximate polylines for major Federal Highways on the Peninsula
HIGHWAYS = [
    (
        "Mex 186 (Escárcega–Xpujil–Chetumal)",
        [
            (-92.0, 17.65),
            (-91.4, 18.00),
            (-90.745, 18.609),  # Escárcega
            (-90.138, 18.592),  # Constitución / Silvituc
            (-89.918, 18.540),  # Conhuás / Balamkú
            (-89.467, 18.517),  # Becán
            (-89.398, 18.506),  # Xpujil
            (-88.929, 18.458),  # Nicolás Bravo
            (-88.480, 18.500),  # Ucum
            (-88.296, 18.500),  # Chetumal
        ],
    ),
    (
        "Mex 180 / 261 (Campeche–Mérida–Cancún)",
        [
            (-91.825, 18.651),  # Cd del Carmen
            (-90.724, 19.355),  # Champotón
            (-90.535, 19.830),  # Campeche
            (-90.050, 20.400),  # Calkiní
            (-89.593, 20.967),  # Mérida
            (-88.568, 20.684),  # Chichén
            (-88.202, 20.690),  # Valladolid
            (-86.852, 21.162),  # Cancún
        ],
    ),
    (
        "Mex 261 / 269 (Escárcega–Champotón & Hopelchén–Xpujil)",
        [
            (-90.745, 18.609),  # Escárcega
            (-90.724, 19.355),  # Champotón
            (-90.535, 19.830),  # Campeche
            (-89.845, 19.746),  # Hopelchén
            (-89.729, 19.461),  # Dzibalchén
            (-89.490, 19.010),  # Ukúm / Xmabén
            (-89.398, 18.506),  # Xpujil
        ],
    ),
    (
        "Mex 307 (Chetumal–Bacalar–Tulum–Cancún)",
        [
            (-88.296, 18.500),  # Chetumal
            (-88.392, 18.678),  # Bacalar
            (-88.045, 19.578),  # F. Carrillo Puerto
            (-87.465, 20.211),  # Tulum
            (-87.074, 20.630),  # Playa del Carmen
            (-86.852, 21.162),  # Cancún
        ],
    ),
    (
        "Mex 184 / 293 (Mérida–Tekax–Morelos–Bacalar)",
        [
            (-89.593, 20.967),  # Mérida
            (-89.286, 20.203),  # Tekax
            (-88.709, 19.749),  # J.M. Morelos
            (-88.450, 18.880),  # Chunhuhub / Buena Vista
            (-88.392, 18.678),  # Bacalar
        ],
    ),
    (
        "Rural Road: Bacalar – Río Verde – Miguel Alemán",
        [
            (-88.392, 18.678),  # Bacalar
            (-88.810, 18.960),  # Río Verde
            (-89.083, 18.922),  # Miguel Alemán
            (-89.139, 18.971),  # Candidate C dirt road
        ],
    ),
]


def fetch_esri_overview(bbox: tuple[float, float, float, float], size: tuple[int, int], cache_name: str):
    cache_path = WORK / "raw_cache" / cache_name
    if cache_path.exists():
        return np.asarray(Image.open(cache_path).convert("RGB"))
    x0, y0, x1, y1 = bbox
    url = (
        f"https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/export?"
        f"bbox={x0:.4f},{y0:.4f},{x1:.4f},{y1:.4f}&bboxSR=4326&imageSR=4326&size={size[0]},{size[1]}&format=png&f=image"
    )
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        raw = urllib.request.urlopen(req, timeout=30).read()
        cache_path.write_bytes(raw)
        return np.asarray(Image.open(io.BytesIO(raw)).convert("RGB"))
    except Exception:
        return None


def main():
    plt.rcParams["font.sans-serif"] = "DejaVu Sans"
    fig = plt.figure(figsize=(24, 15.5), dpi=180, facecolor="#090f1a")
    gs = GridSpec(
        2,
        2,
        figure=fig,
        width_ratios=[1.05, 1.15],
        height_ratios=[1.0, 1.08],
        left=0.035,
        right=0.98,
        top=0.925,
        bottom=0.035,
        wspace=0.11,
        hspace=0.15,
    )

    # Suptitle
    fig.text(
        0.035,
        0.968,
        "WHERE THE NEW MAYA LiDAR DISCOVERIES SIT IN MEXICO & ON THE GLOBE",
        color="#f8fafc",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.035,
        0.944,
        "Continental & National Context  •  Modern Cities & Federal Highways (Mex 186 / 180 / 261 / 307)  •  456 NASA G-LiHT Flight Strips  •  6 Unmapped Maya Complexes (NC-1 through NC-6)",
        color="#94a3b8",
        fontsize=11.5,
    )

    # =========================================================================
    # PANEL 1 (TOP-LEFT): ALL OF MEXICO & CENTRAL AMERICA + 3D GLOBE INSET
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_facecolor("#071324")
    ax1.set_xlim(-117.5, -76.5)
    ax1.set_ylim(13.2, 33.2)

    for f in ne0["features"]:
        adm = f["properties"].get("ADMIN", "")
        if adm == "Mexico":
            plot_geom(ax1, f["geometry"], fc="#1e293b", ec="#64748b", lw=0.9, zorder=2)
        else:
            plot_geom(ax1, f["geometry"], fc="#111b29", ec="#334155", lw=0.6, zorder=1)

    # Mexican state boundaries on Panel 1
    yuc_states = {"Campeche": "#1e3a5f", "Yucatán": "#1b3354", "Quintana Roo": "#173152"}
    for f in mex["features"]:
        nm = f["properties"].get("name", "")
        fc = yuc_states.get(nm, "#1e293b")
        ec = "#60a5fa" if nm in yuc_states else "#475569"
        lw = 0.9 if nm in yuc_states else 0.45
        plot_geom(ax1, f["geometry"], fc=fc, ec=ec, lw=lw, zorder=3)

    # Major cities across Mexico & Central America
    for cname, clat, clon, is_big in MODERN_CITIES_MEX:
        sz = 42 if is_big else 22
        col = "#ffffff" if is_big else "#cbd5e1"
        ax1.scatter([clon], [clat], s=sz, c=col, edgecolors="#0f172a", linewidths=1.0, zorder=6)
        dx, dy = 0.4, 0.25
        if "Mexico City" in cname:
            dx, dy = -4.8, -0.85
        elif "Campeche" in cname:
            dx, dy = -3.8, 0.35
        elif "Mérida" in cname:
            dx, dy = -2.4, 0.45
        elif "Chetumal" in cname:
            dx, dy = 0.45, -0.45
        elif "Villahermosa" in cname:
            dx, dy = -4.2, -0.55
        elif "Veracruz" in cname:
            dx, dy = -3.2, 0.35
        ax1.text(
            clon + dx,
            clat + dy,
            cname,
            color="#f1f5f9" if is_big else "#94a3b8",
            fontsize=8.5 if is_big else 7.5,
            fontweight="bold" if is_big else "normal",
            path_effects=stroke_med,
            zorder=7,
        )

    # Plot discovery stars on Mexico map
    for d in DISCOVERIES:
        if "NEW" in d["status"]:
            ax1.scatter([d["lon"]], [d["lat"]], marker="*", s=130, c=d["color"], edgecolors="white", linewidths=0.7, zorder=8)

    # Bounding boxes for Panel 2 (Peninsula) and Panel 3 (Heartland)
    ax1.add_patch(
        Rectangle((-92.2, 16.5), 5.6, 5.1, fill=False, edgecolor="#38bdf8", lw=1.8, ls="--", zorder=7)
    )
    ax1.text(-92.1, 21.8, "PANEL 2: YUCATÁN PENINSULA", color="#38bdf8", fontsize=8.5, fontweight="bold", path_effects=stroke_med)

    ax1.add_patch(
        Rectangle((-91.15, 17.85), 3.0, 1.6, fill=False, edgecolor="#facc15", lw=2.0, ls="-", zorder=8)
    )
    ax1.text(-91.15, 17.35, "PANEL 3: HEARTLAND ZOOM", color="#facc15", fontsize=8.0, fontweight="bold", path_effects=stroke_med)

    # Water labels
    ax1.text(-95.0, 24.5, "GULF OF\nMEXICO", color="#38bdf8", alpha=0.55, fontsize=10, fontweight="bold", ha="center")
    ax1.text(-108.0, 17.5, "PACIFIC OCEAN", color="#38bdf8", alpha=0.55, fontsize=10, fontweight="bold", ha="center")
    ax1.text(-83.5, 19.2, "CARIBBEAN\nSEA", color="#38bdf8", alpha=0.55, fontsize=9, fontweight="bold", ha="center")
    ax1.text(-102.5, 28.5, "MEXICO", color="#64748b", alpha=0.55, fontsize=15, fontweight="bold", ha="center")
    ax1.text(-99.0, 31.3, "UNITED STATES", color="#475569", alpha=0.55, fontsize=11, fontweight="bold", ha="center")

    ax1.set_title(
        "1. CONTINENTAL & NATIONAL CONTEXT — WHERE THE YUCATÁN PENINSULA SITS IN MEXICO",
        color="#e2e8f0",
        fontsize=11.5,
        fontweight="bold",
        loc="left",
        pad=8,
    )
    ax1.set_xlabel("Longitude (°W)", color="#94a3b8", fontsize=9)
    ax1.set_ylabel("Latitude (°N)", color="#94a3b8", fontsize=9)
    ax1.tick_params(colors="#94a3b8", labelsize=8)
    ax1.grid(color="#1e293b", linestyle=":", linewidth=0.6)

    # Inset 3D Orthographic Globe in lower-left corner of Panel 1
    ax_globe = fig.add_axes([0.042, 0.545, 0.125, 0.185])
    plot_ortho_globe(ax_globe, center_lon=-92.0, center_lat=20.0)

    # =========================================================================
    # PANEL 2 (TOP-RIGHT): FULL YUCATÁN PENINSULA, MODERN CITIES & HIGHWAYS
    # =========================================================================
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_facecolor("#071526")
    p2_bbox = (-92.2, 16.5, -86.6, 21.65)
    ax2.set_xlim(p2_bbox[0], p2_bbox[2])
    ax2.set_ylim(p2_bbox[1], p2_bbox[3])

    for f in ne0["features"]:
        adm = f["properties"].get("ADMIN", "")
        if adm in ("Guatemala", "Belize", "Honduras"):
            plot_geom(ax2, f["geometry"], fc="#141e2e", ec="#475569", lw=0.9, zorder=1)

    st_fills = {
        "Campeche": "#1d2b3a",
        "Yucatán": "#1a2636",
        "Quintana Roo": "#182838",
        "Tabasco": "#151f2d",
        "Chiapas": "#151f2d",
    }
    for f in mex["features"]:
        nm = f["properties"].get("name", "")
        plot_geom(ax2, f["geometry"], fc=st_fills.get(nm, "#151f2d"), ec="#94a3b8", lw=1.2, zorder=2)

    # State labels
    ax2.text(-90.35, 19.15, "CAMPECHE", color="#94a3b8", alpha=0.45, fontsize=13, fontweight="bold", ha="center", zorder=3)
    ax2.text(-88.95, 20.55, "YUCATÁN", color="#94a3b8", alpha=0.45, fontsize=13, fontweight="bold", ha="center", zorder=3)
    ax2.text(-87.95, 19.35, "QUINTANA\nROO", color="#94a3b8", alpha=0.45, fontsize=12, fontweight="bold", ha="center", zorder=3)
    ax2.text(-90.25, 16.95, "GUATEMALA (PETÉN)", color="#64748b", alpha=0.5, fontsize=10, fontweight="bold", ha="center", zorder=3)
    ax2.text(-88.55, 17.45, "BELIZE", color="#64748b", alpha=0.5, fontsize=10, fontweight="bold", ha="center", zorder=3)

    # Plot 5,223 EAAMS sites as subtle amber background cloud
    eaams_xy = np.array([f["geometry"]["coordinates"][:2] for f in eaams["features"] if f.get("geometry", {}).get("type") == "Point"])
    ax2.scatter(eaams_xy[:, 0], eaams_xy[:, 1], s=2.0, c="#fbbf24", alpha=0.22, linewidths=0, zorder=3)

    # Plot all 456 G-LiHT transects
    for fp in footprints:
        b = fp["bbox_wgs84"]  # [min_lon, min_lat, max_lon, max_lat]
        cx, cy = (b[0] + b[2]) * 0.5, (b[1] + b[3]) * 0.5
        w = max(b[2] - b[0], 0.03)
        h = max(b[3] - b[1], 0.03)
        ax2.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h, facecolor="#38bdf8", edgecolor="#38bdf8", lw=0.6, alpha=0.65, zorder=4))

    # Plot Federal Highways
    for hname, pts in HIGHWAYS:
        arr = np.asarray(pts)
        is_dirt = "Rural" in hname
        ax2.plot(
            arr[:, 0],
            arr[:, 1],
            color="#f59e0b" if not is_dirt else "#cbd5e1",
            lw=1.5 if not is_dirt else 1.0,
            ls="-" if not is_dirt else ":",
            alpha=0.80,
            zorder=5,
        )

    # Plot Modern Cities on Peninsula
    for cname, clat, clon, tier in PENINSULA_CITIES:
        if not (p2_bbox[0] <= clon <= p2_bbox[2] and p2_bbox[1] <= clat <= p2_bbox[3]):
            continue
        if tier == "major":
            ax2.scatter([clon], [clat], marker="s", s=38, c="#ffffff", edgecolors="#0f172a", linewidths=1.1, zorder=7)
            ax2.text(clon + 0.06, clat + 0.04, cname, color="#ffffff", fontsize=8.5, fontweight="bold", path_effects=stroke_med, zorder=8)
        elif tier == "med":
            ax2.scatter([clon], [clat], marker="o", s=22, c="#e2e8f0", edgecolors="#0f172a", linewidths=0.8, zorder=6)
            ax2.text(clon + 0.05, clat + 0.03, cname, color="#cbd5e1", fontsize=7.5, path_effects=stroke_light, zorder=7)

    # Plot Major Famous Maya Ruins
    for rname, rlat, rlon in FAMOUS_MAYA_RUINS:
        if rname in ("Chichén Itzá", "Uxmal", "Calakmul", "Cobá", "Tulum Ruins", "Edzná", "Santa Rosa Xtampak", "Becán", "Ocomtún (2023)", "Valeriana (2024)", "Tikal"):
            ax2.scatter([rlon], [rlat], marker="^", s=46, c="#fde047", edgecolors="#0f172a", linewidths=0.9, zorder=7)
            dy = -0.11 if rname in ("Calakmul", "Becán", "Uxmal") else 0.05
            ax2.text(rlon + 0.05, rlat + dy, rname, color="#fef08a", fontsize=7.5, fontweight="bold", path_effects=stroke_light, zorder=8)

    # Plot Our Discoveries
    offsets_p2 = {
        "NC-1": (-0.48, -0.16),
        "NC-2": (-0.52, -0.05),
        "NC-3": (0.08, 0.06),
        "NC-4": (-0.52, 0.07),
        "NC-5": (-0.52, 0.11),
        "NC-6": (0.08, -0.05),
        "Cand D": (0.08, -0.15),
        "Cand E": (0.08, -0.12),
    }
    for d in DISCOVERIES:
        is_new = "NEW" in d["status"]
        mk = "*" if is_new else "D"
        sz = 190 if is_new else 70
        ax2.scatter([d["lon"]], [d["lat"]], marker=mk, s=sz, c=d["color"], edgecolors="white", linewidths=1.0, zorder=10)
        ox, oy = offsets_p2.get(d["code"], (0.07, 0.04))
        ax2.text(
            d["lon"] + ox,
            d["lat"] + oy,
            d["code"],
            color=d["color"],
            fontsize=8.5,
            fontweight="bold",
            path_effects=stroke,
            zorder=11,
        )

    # Box showing Panel 3 zoom
    ax2.add_patch(Rectangle((-91.15, 17.85), 3.0, 1.6, fill=False, edgecolor="#facc15", lw=1.8, ls="--", zorder=9))

    ax2.set_title(
        "2. YUCATÁN PENINSULA — MODERN CITIES (WHITE), HIGHWAYS (ORANGE), FAMOUS RUINS (▲) & NEW FINDS (★)",
        color="#e2e8f0",
        fontsize=11.5,
        fontweight="bold",
        loc="left",
        pad=8,
    )
    ax2.set_xlabel("Longitude (°W)", color="#94a3b8", fontsize=9)
    ax2.set_ylabel("Latitude (°N)", color="#94a3b8", fontsize=9)
    ax2.tick_params(colors="#94a3b8", labelsize=8)
    ax2.grid(color="#1e293b", linestyle=":", linewidth=0.6)

    # =========================================================================
    # PANEL 3 (BOTTOM-LEFT): HEARTLAND SATELLITE + HIGHWAY + CITY ZOOM
    # =========================================================================
    ax3 = fig.add_subplot(gs[1, 0])
    p3_bbox = (-91.15, 17.85, -88.15, 19.45)
    ax3.set_xlim(p3_bbox[0], p3_bbox[2])
    ax3.set_ylim(p3_bbox[1], p3_bbox[3])

    sat_bg = fetch_esri_overview(p3_bbox, (1200, 640), "esri_heartland_overview.png")
    if sat_bg is not None:
        # Darken slightly so labels/roads pop clearly while showing dark-green Biosphere forest vs farm clearings
        sat_dim = (sat_bg.astype(np.float32) * 0.72).astype(np.uint8)
        ax3.imshow(sat_dim, extent=[p3_bbox[0], p3_bbox[2], p3_bbox[1], p3_bbox[3]], origin="upper", aspect="auto", zorder=1)
    else:
        ax3.set_facecolor("#111c2d")

    # State boundaries on Panel 3
    for f in mex["features"]:
        plot_geom(ax3, f["geometry"], fc="none", ec="#cbd5e1", lw=1.6, ls="-.", zorder=3)

    ax3.text(-89.10, 19.34, "← CAMPECHE  |  QUINTANA ROO →", color="#e2e8f0", fontsize=8.5, fontweight="bold", ha="center", path_effects=stroke, zorder=6)

    # Plot G-LiHT strips in Panel 3
    for fp in footprints:
        b = fp["bbox_wgs84"]
        if b[2] < p3_bbox[0] or b[0] > p3_bbox[2] or b[3] < p3_bbox[1] or b[1] > p3_bbox[3]:
            continue
        cx, cy = (b[0] + b[2]) * 0.5, (b[1] + b[3]) * 0.5
        w = max(b[2] - b[0], 0.022)
        h = max(b[3] - b[1], 0.022)
        ax3.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h, facecolor="#38bdf8", edgecolor="#38bdf8", lw=0.6, alpha=0.50, zorder=3))

    # Plot Highways on Panel 3
    for hname, pts in HIGHWAYS:
        arr = np.asarray(pts)
        is_dirt = "Rural" in hname
        ax3.plot(
            arr[:, 0],
            arr[:, 1],
            color="#f59e0b" if not is_dirt else "#e2e8f0",
            lw=2.2 if not is_dirt else 1.3,
            ls="-" if not is_dirt else "--",
            alpha=0.92,
            zorder=4,
        )
    ax3.text(-90.42, 18.63, "Fed Hwy Mex 186", color="#fbbf24", fontsize=8.0, fontweight="bold", path_effects=stroke, zorder=7)
    ax3.text(-88.75, 18.42, "Fed Hwy Mex 186 → Chetumal", color="#fbbf24", fontsize=8.0, fontweight="bold", path_effects=stroke, zorder=7)
    ax3.text(-89.66, 19.22, "Hwy 269\n(Hopelchén–Xpujil)", color="#fbbf24", fontsize=7.5, fontweight="bold", path_effects=stroke, zorder=7)

    # Modern Cities & Towns in Panel 3
    for cname, clat, clon, tier in PENINSULA_CITIES:
        if not (p3_bbox[0] <= clon <= p3_bbox[2] and p3_bbox[1] <= clat <= p3_bbox[3]):
            continue
        sz = 48 if tier == "major" else 26
        ax3.scatter([clon], [clat], marker="s", s=sz, c="#ffffff", edgecolors="#0f172a", linewidths=1.1, zorder=7)
        ax3.text(
            clon + 0.03,
            clat + 0.025,
            cname,
            color="#ffffff",
            fontsize=8.2 if tier == "major" else 7.4,
            fontweight="bold",
            path_effects=stroke,
            zorder=8,
        )

    # Famous Maya Ruins in Panel 3
    for rname, rlat, rlon in FAMOUS_MAYA_RUINS:
        if not (p3_bbox[0] <= rlon <= p3_bbox[2] and p3_bbox[1] <= rlat <= p3_bbox[3]):
            continue
        ax3.scatter([rlon], [rlat], marker="^", s=55, c="#fde047", edgecolors="#0f172a", linewidths=1.0, zorder=8)
        dy = -0.055 if rname in ("Becán", "Chakanbakán", "Balamkú") else 0.03
        ax3.text(rlon + 0.03, rlat + dy, rname, color="#fef08a", fontsize=7.8, fontweight="bold", path_effects=stroke, zorder=9)

    # Plot Our Discoveries with prominent callouts on Panel 3
    p3_offsets = {
        "NC-1": (-0.42, -0.13),
        "NC-2": (-0.44, -0.12),
        "NC-3": (-0.45, 0.06),
        "NC-4": (-0.38, 0.10),
        "NC-5": (-0.45, -0.02),
        "NC-6": (-0.42, 0.07),
        "Cand D": (0.05, -0.08),
        "Cand E": (0.05, -0.09),
    }
    for d in DISCOVERIES:
        is_new = "NEW" in d["status"]
        mk = "*" if is_new else "D"
        sz = 240 if is_new else 85
        ax3.scatter([d["lon"]], [d["lat"]], marker=mk, s=sz, c=d["color"], edgecolors="white", linewidths=1.2, zorder=11)
        ox, oy = p3_offsets[d["code"]]
        tx, ty = d["lon"] + ox, d["lat"] + oy
        ax3.annotate(
            d["name"],
            xy=(d["lon"], d["lat"]),
            xytext=(tx, ty),
            color="#ffffff" if is_new else "#e2e8f0",
            fontsize=8.2 if is_new else 7.5,
            fontweight="bold",
            bbox=dict(
                boxstyle="round,pad=0.25",
                facecolor="#0f172a" if is_new else "#1e293b",
                edgecolor=d["color"],
                lw=1.5 if is_new else 1.0,
                alpha=0.92,
            ),
            arrowprops=dict(arrowstyle="->", color=d["color"], lw=1.3),
            zorder=12,
        )

    ax3.set_title(
        "3. HEARTLAND SATELLITE ZOOM — DARK GREEN = UNBROKEN JUNGLE (CALAKMUL / BALAMKÚ); TAN = FARMS & TOWNS",
        color="#e2e8f0",
        fontsize=11.2,
        fontweight="bold",
        loc="left",
        pad=8,
    )
    ax3.set_xlabel("Longitude (°W)", color="#94a3b8", fontsize=9)
    ax3.set_ylabel("Latitude (°N)", color="#94a3b8", fontsize=9)
    ax3.tick_params(colors="#94a3b8", labelsize=8)
    ax3.grid(color="#334155", linestyle=":", linewidth=0.5, alpha=0.6)

    # =========================================================================
    # PANEL 4 (BOTTOM-RIGHT): GEOGRAPHIC & MODERN-CITY REFERENCE TABLE
    # =========================================================================
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.set_facecolor("#0f172a")
    ax4.set_xlim(0, 100)
    ax4.set_ylim(0, 100)
    ax4.axis("off")

    ax4.set_title(
        "4. LOCATION SCORECARD — WHERE EACH COMPLEX SITS RELATIVE TO MODERN CITIES, HIGHWAYS & FAMOUS RUINS",
        color="#e2e8f0",
        fontsize=11.2,
        fontweight="bold",
        loc="left",
        pad=8,
    )

    # Table headers
    headers = [
        (1.5, "SITE ID & STATUS"),
        (22.0, "STATE & MUNICIPALITY"),
        (43.5, "NEAREST MODERN CITIES & HIGHWAYS"),
        (74.5, "NEAREST KNOWN MAYA RUINS"),
    ]
    ax4.add_patch(Rectangle((0.5, 91.5), 99.0, 7.0, facecolor="#1e293b", edgecolor="#475569", lw=1.0))
    for hx, htxt in headers:
        ax4.text(hx, 94.3, htxt, color="#f8fafc", fontsize=8.8, fontweight="bold")

    y_top = 89.5
    row_h = 13.8
    for idx, d in enumerate(DISCOVERIES):
        y0 = y_top - (idx + 1) * row_h
        is_new = "NEW" in d["status"]
        bg = "#132036" if idx % 2 == 0 else "#0d1829"
        ec = d["color"] if is_new else "#334155"
        ax4.add_patch(Rectangle((0.5, y0), 99.0, row_h - 0.8, facecolor=bg, edgecolor=ec, lw=1.1 if is_new else 0.6))

        # Col 1: Code + Status
        ax4.scatter([2.6], [y0 + row_h * 0.58], marker="*" if is_new else "D", s=110 if is_new else 55, c=d["color"], edgecolors="white", linewidths=0.7)
        ax4.text(4.5, y0 + row_h * 0.56, d["code"], color=d["color"], fontsize=9.8, fontweight="bold")
        ax4.text(4.5, y0 + row_h * 0.18, d["status"], color="#86efac" if is_new else "#93c5fd", fontsize=7.3, fontweight="bold")

        # Col 2: State & Municipality
        ax4.text(22.0, y0 + row_h * 0.54, d["state_mun"], color="#f1f5f9", fontsize=8.2, fontweight="bold")
        ax4.text(
            22.0,
            y0 + row_h * 0.18,
            f"Approx: ~{d['lat']:.1f}°N, {abs(d['lon']):.1f}°W",
            color="#94a3b8",
            fontsize=7.5,
        )

        # Col 3: Modern Cities & Highways
        ax4.text(43.5, y0 + row_h * 0.52, d["modern_ref"], color="#e2e8f0", fontsize=7.8)
        ax4.text(
            43.5,
            y0 + row_h * 0.16,
            f"Distance to Mérida: {hav_km(d['lat'], d['lon'], 20.9674, -89.5926):.0f} km  |  To Cancún: {hav_km(d['lat'], d['lon'], 21.1619, -86.8515):.0f} km  |  To Chetumal: {hav_km(d['lat'], d['lon'], 18.5001, -88.2961):.0f} km",
            color="#94a3b8",
            fontsize=7.1,
        )

        # Col 4: Nearest Known Maya Ruins
        ax4.text(74.5, y0 + row_h * 0.40, d["ruin_ref"], color="#fde047" if is_new else "#cbd5e1", fontsize=7.6)

    # Footer legend note inside Panel 4
    ax4.text(
        1.0,
        1.5,
        "★ = Confirmed Unmapped Maya Complex (Absent from EAAMS 5,223 sites, INAH SIC, Šprajc & Auld-Thomas)      ■ = Modern City / Town      ▲ = Famous Published Maya City",
        color="#cbd5e1",
        fontsize=8.2,
        fontweight="bold",
    )

    out_path = MEDIA / "maya-discoveries-mexico-and-globe-map.png"
    fig.savefig(out_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Wrote {out_path} ({out_path.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
