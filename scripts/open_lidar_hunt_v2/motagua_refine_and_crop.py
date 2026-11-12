#!/usr/bin/env python3
"""Rescore Motagua candidates for architectural (not tectonic) relief and crop zooms."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image, ImageDraw
from rasterio.transform import rowcol, xy as rio_xy
from rasterio.warp import transform as rio_transform
from scipy import ndimage as ndi

WORK = Path("/agent/lidar/work/open_lidar_hunt_v2")
DTM_DIR = WORK / "motagua/dtm"
OUT = WORK / "out/motagua"
PREV = WORK / "previews/motagua/crops"
PREV.mkdir(parents=True, exist_ok=True)

KNOWN = [
    {"name": "Quirigua", "lat": 15.2694, "lon": -89.0404},
]


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(min(1.0, a)))


def lrm_rgb(lrm, valid):
    v = np.nan_to_num(lrm, nan=0.0)
    lo, hi = np.percentile(v[valid], [5, 95]) if valid.any() else (-1, 1)
    t = np.clip((v - lo) / max(hi - lo, 1e-6), 0, 1)
    rgb = np.zeros(v.shape + (3,), dtype=np.uint8)
    rgb[..., 0] = (30 + 200 * t).astype(np.uint8)
    rgb[..., 1] = (40 + 160 * (1 - np.abs(t - 0.5) * 2)).astype(np.uint8)
    rgb[..., 2] = (40 + 200 * (1 - t)).astype(np.uint8)
    rgb[~valid] = 15
    return rgb


def norm_lrm(dem, valid, sigma_px):
    step = 4 if sigma_px >= 12 else 1
    h, w = dem.shape
    sub_v = valid[::step, ::step].astype(np.float32)
    sub_d = np.where(valid[::step, ::step], dem[::step, ::step], 0.0).astype(np.float32)
    num = ndi.gaussian_filter(sub_d, sigma=sigma_px / step, mode="constant", cval=0.0)
    den = ndi.gaussian_filter(sub_v, sigma=sigma_px / step, mode="constant", cval=0.0)
    sub_tr = np.where(den > 1e-3, num / np.maximum(den, 1e-3), 0.0).astype(np.float32)
    trend = ndi.zoom(sub_tr, (h / sub_tr.shape[0], w / sub_tr.shape[1]), order=1)[:h, :w]
    lrm = dem - trend
    lrm[~valid] = np.nan
    return lrm


def arch_score(c):
    """Down-weight tectonic hills; prefer 1–8 m orthogonal plaza groups."""
    h = c["h_max"]
    if h > 12.0:
        return -1.0  # natural ridgeline / fault scarp
    if h < 0.9:
        return -1.0
    # sweet spot 1.2–6.5 m platforms
    h_term = 6.0 * h if h <= 6.5 else 6.0 * 6.5 - 4.0 * (h - 6.5)
    s = (
        h_term
        + 4.0 * c["n_mounds"]
        + 10.0 * c.get("rect_mean", 0)
        + 12.0 * c.get("ra_mean", 0)
        + 3.0 * c.get("quads", 1)
        + min(8.0, c.get("plaza_m2", 0) / 100.0)
    )
    if c.get("ra_mean", 0) < 0.18 and c.get("rect_mean", 0) < 0.45 and c["n_mounds"] < 4:
        s *= 0.35
    if c.get("solitary") and c.get("ra_mean", 0) < 0.25:
        s *= 0.5
    return float(s)


def crop_candidate(ds_path, lat, lon, tag, half_m=90):
    with rasterio.open(ds_path) as ds:
        xs, ys = rio_transform("EPSG:4326", ds.crs, [lon], [lat])
        r, c = rowcol(ds.transform, xs[0], ys[0])
        res = abs(ds.res[0])
        half = int(half_m / res)
        r0, r1 = max(0, r - half), min(ds.height, r + half)
        c0, c1 = max(0, c - half), min(ds.width, c + half)
        dem = ds.read(1, window=((r0, r1), (c0, c1))).astype(np.float32)
        nodata = ds.nodata
        transform = ds.window_transform(((r0, r1), (c0, c1)))
    valid = np.isfinite(dem)
    if nodata is not None:
        valid &= dem != nodata
    if valid.sum() < 100:
        return None
    dem = np.where(valid, dem, np.nan)
    lrm = norm_lrm(dem, valid, 18.0 / abs(transform.a))
    rgb = lrm_rgb(lrm, valid)
    # hillshade overlay
    gy, gx = np.gradient(np.nan_to_num(dem, nan=np.nanmedian(dem[valid])))
    hs = gx * 0.5 + gy * 0.7
    hs = (hs - np.nanmin(hs)) / max(np.nanmax(hs) - np.nanmin(hs), 1e-6)
    for i in range(3):
        rgb[..., i] = (rgb[..., i] * (0.55 + 0.45 * hs)).astype(np.uint8)
    im = Image.fromarray(rgb)
    draw = ImageDraw.Draw(im)
    cy, cx = dem.shape[0] // 2, dem.shape[1] // 2
    draw.ellipse([cx - 8, cy - 8, cx + 8, cy + 8], outline=(255, 255, 0), width=2)
    draw.text((4, 4), tag, fill=(255, 255, 0))
    out = PREV / f"{tag}.png"
    im.save(out)
    return str(out)


def main():
    raw = json.loads((OUT / "motagua_candidates.json").read_text())
    scored = []
    for c in raw:
        # skip military base for public archaeological claims
        if "military" in c["dataset"]:
            continue
        s = arch_score(c)
        if s < 12.0:
            continue
        nearest = min(KNOWN, key=lambda k: haversine_km(c["lat"], c["lon"], k["lat"], k["lon"]))
        d = haversine_km(c["lat"], c["lon"], nearest["lat"], nearest["lon"])
        scored.append({**c, "arch_score": s, "dist_quirigua_km": d})
    scored.sort(key=lambda x: x["arch_score"], reverse=True)
    # NMS 80 m
    kept = []
    for c in scored:
        if any(haversine_km(c["lat"], c["lon"], k["lat"], k["lon"]) * 1000 < 80 for k in kept):
            continue
        kept.append(c)
    # split control vs new
    control = [c for c in kept if c["dist_quirigua_km"] < 0.6]
    novel = [c for c in kept if c["dist_quirigua_km"] >= 0.6]
    # crop top novel + control
    crops = []
    for i, c in enumerate(control[:4] + novel[:16], start=1):
        tag = f"{'QC' if c['dist_quirigua_km'] < 0.6 else 'MG'}-{i:02d}_{c['dataset'].replace('_dtm_cog.tif','')}"
        path = DTM_DIR / c["dataset"]
        crop = crop_candidate(path, c["lat"], c["lon"], tag)
        c["crop"] = crop
        c["id"] = tag
        crops.append(c)
        print(f"{tag} arch={c['arch_score']:.1f} h={c['h_max']:.1f} n={c['n_mounds']} "
              f"{c['lat']:.5f},{c['lon']:.5f} quirigua={c['dist_quirigua_km']:.2f}km")
    (OUT / "motagua_architectural.json").write_text(json.dumps(crops, indent=2))
    (OUT / "motagua_novel_top.json").write_text(json.dumps(novel[:30], indent=2))
    print("control", len(control), "novel kept", len(novel), "cropped", len(crops))


if __name__ == "__main__":
    main()
