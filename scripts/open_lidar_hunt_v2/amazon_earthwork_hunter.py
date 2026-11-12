#!/usr/bin/env python3
"""Amazon earthwork hunter for SFB DTMs and EBA ground-classified LAZ strips.

Targets geometric ditched enclosures (geoglyphs), ring ditches, and mound
clusters using local-relief + circular/rectangular ditch detectors — distinct
from the Maya plaza grammar used on Motagua / Yucatán.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import laspy
import numpy as np
import rasterio
from PIL import Image, ImageDraw
from pyproj import Transformer
from rasterio.transform import Affine, from_origin, xy as rio_xy
from rasterio.warp import transform as rio_transform
from scipy import ndimage as ndi

WORK = Path("/agent/lidar/work/open_lidar_hunt_v2")
OUT = WORK / "out/amazon"
PREV = WORK / "previews/amazon"
DTM_OUT = WORK / "amazon_eba/dtm"
OUT.mkdir(parents=True, exist_ok=True)
PREV.mkdir(parents=True, exist_ok=True)
DTM_OUT.mkdir(parents=True, exist_ok=True)

GEO = []
with open(WORK / "known_sites/amazon_geoglyphs_points.csv") as f:
    for row in csv.DictReader(f):
        GEO.append((float(row["lat"]), float(row["lon"]), row["name"]))


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(min(1.0, a)))


def nearest_geo(lat, lon):
    best = min(GEO, key=lambda g: haversine_km(lat, lon, g[0], g[1]))
    return best[2], haversine_km(lat, lon, best[0], best[1]), best[0], best[1]


def norm_lrm(dem, valid, sigma_px):
    h, w = dem.shape
    step = 4 if sigma_px >= 10 else 1
    sub_v = valid[::step, ::step].astype(np.float32)
    sub_d = np.where(valid[::step, ::step], dem[::step, ::step], 0.0).astype(np.float32)
    num = ndi.gaussian_filter(sub_d, sigma=sigma_px / step, mode="constant", cval=0.0)
    den = ndi.gaussian_filter(sub_v, sigma=sigma_px / step, mode="constant", cval=0.0)
    sub_tr = np.where(den > 1e-3, num / np.maximum(den, 1e-3), 0.0).astype(np.float32)
    trend = ndi.zoom(sub_tr, (h / sub_tr.shape[0], w / sub_tr.shape[1]), order=1)[:h, :w]
    lrm = dem - trend
    lrm[~valid] = np.nan
    return lrm


def laz_to_dtm(laz_path: Path, res_m: float = 1.0) -> Path:
    out = DTM_OUT / f"{laz_path.stem}_dtm_{res_m:.0f}m.tif"
    if out.exists():
        return out
    print(f"GRID {laz_path.name}", flush=True)
    las = laspy.read(laz_path)
    cls = np.array(las.classification)
    ground = cls == 2
    if ground.sum() < 5000:
        # fallback: lowest 15% of points in cells later; use all for now
        x = np.array(las.x)
        y = np.array(las.y)
        z = np.array(las.z)
    else:
        x = np.array(las.x)[ground]
        y = np.array(las.y)[ground]
        z = np.array(las.z)[ground]
    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())
    # slightly pad
    cols = int(math.ceil((xmax - xmin) / res_m)) + 1
    rows = int(math.ceil((ymax - ymin) / res_m)) + 1
    # bin min-Z
    ci = np.clip(((x - xmin) / res_m).astype(np.int32), 0, cols - 1)
    ri = np.clip(((ymax - y) / res_m).astype(np.int32), 0, rows - 1)
    dem = np.full((rows, cols), np.inf, dtype=np.float32)
    flat = ri * cols + ci
    # scatter min via sort
    order = np.argsort(flat)
    flat_s = flat[order]
    z_s = z[order]
    # first occurrence after sort of keys with z ascending within key
    # better: sort by flat then z
    order2 = np.lexsort((z_s, flat_s))
    flat_s = flat_s[order2]
    z_s = z_s[order2]
    uniq, idx = np.unique(flat_s, return_index=True)
    dem.ravel()[uniq] = z_s[idx]
    dem[np.isinf(dem)] = np.nan
    # fill small holes
    mask = np.isfinite(dem)
    if mask.any():
        # nearest fill for gaps < 5 cells via binary closing on mask then IDW-ish
        for _ in range(2):
            filled = ndi.generic_filter(
                np.nan_to_num(dem, nan=0.0),
                function=lambda a: np.nanmean(a[a != 0]) if np.any(a != 0) else np.nan,
                size=3,
                mode="constant",
                cval=0.0,
            )
            dem = np.where(mask, dem, filled)
            mask = np.isfinite(dem)
    transform = from_origin(xmin, ymax, res_m, res_m)
    crs = las.header.parse_crs()
    with rasterio.open(
        out,
        "w",
        driver="GTiff",
        height=rows,
        width=cols,
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        compress="deflate",
        nodata=np.nan,
    ) as dst:
        dst.write(dem.astype(np.float32), 1)
    print(f"  wrote {out} {cols}x{rows}", flush=True)
    return out


def circularity(comp_mask, res_m):
    ys, xs = np.where(comp_mask)
    if len(xs) < 40:
        return 0.0, 0.0, 0.0
    cy, cx = ys.mean(), xs.mean()
    dist = np.hypot(ys - cy, xs - cx) * res_m
    r_med = float(np.median(dist))
    if r_med < 15:
        return 0.0, r_med, 0.0
    # ring thickness score: points near median radius
    ring = np.abs(dist - r_med) < max(8.0, 0.18 * r_med)
    circ = float(ring.mean())
    # perimeter regularity
    return circ, r_med, float(dist.std())


def detect_earthworks(dem_path: Path, source: str):
    print(f"SCAN {dem_path}", flush=True)
    with rasterio.open(dem_path) as ds:
        dem = ds.read(1).astype(np.float32)
        transform = ds.transform
        crs = ds.crs
        res_m = float(abs(ds.res[0]))
        nodata = ds.nodata
    valid = np.isfinite(dem)
    if nodata is not None and np.isfinite(nodata):
        valid &= dem != nodata
    if valid.sum() < 1000:
        return []
    dem = np.where(valid, dem, np.nan)
    interior = ndi.binary_erosion(valid, iterations=max(2, int(3 / res_m)))
    lrm = norm_lrm(dem, valid, max(25.0 / res_m, 8.0))
    gy, gx = np.gradient(np.nan_to_num(dem, nan=np.nanmedian(dem[valid])))
    slope = np.degrees(np.arctan(np.hypot(gx, gy) / res_m))
    # ditches = negative LRM rings
    ditch = (lrm <= -0.55) & interior
    ditch = ndi.binary_closing(ditch, iterations=1)
    ditch = ndi.binary_opening(ditch, iterations=1)
    lab, n = ndi.label(ditch)
    hits = []
    for i, slc in enumerate(ndi.find_objects(lab), start=1):
        if slc is None:
            continue
        comp = lab[slc] == i
        area = float(comp.sum() * res_m * res_m)
        if area < 800 or area > 120000:
            continue
        circ, r_med, r_std = circularity(comp, res_m)
        # diameter ~ 2*r
        diam = 2 * r_med
        if diam < 40 or diam > 450:
            continue
        # geoglyphs often have circ>0.25 for ring-like masks; also allow filled ditch squares
        ys, xs = np.where(comp)
        cy = slc[0].start + ys.mean()
        cx = slc[1].start + xs.mean()
        depth = float(-np.nanpercentile(lrm[slc][comp], 10))
        # rectangularity of bbox
        h_m = (slc[0].stop - slc[0].start) * res_m
        w_m = (slc[1].stop - slc[1].start) * res_m
        rect_fill = area / max(h_m * w_m, 1.0)
        # score
        score = 8.0 * depth + 6.0 * circ + min(8.0, diam / 40.0) + 3.0 * rect_fill
        if circ < 0.12 and rect_fill < 0.25:
            continue
        if score < 8.0:
            continue
        x, y = rio_xy(transform, cy, cx)
        lon, lat = rio_transform(crs, "EPSG:4326", [x], [y])
        lat, lon = float(lat[0]), float(lon[0])
        gname, gdkm, glat, glon = nearest_geo(lat, lon)
        hits.append(
            {
                "lat": lat,
                "lon": lon,
                "score": float(score),
                "depth_m": depth,
                "diameter_m": float(diam),
                "circ": float(circ),
                "area_m2": area,
                "rect_fill": float(rect_fill),
                "source": source,
                "dataset": dem_path.name,
                "nearest_geoglyph": gname,
                "dist_known_km": float(gdkm),
                "maps": f"https://www.google.com/maps?q={lat:.6f},{lon:.6f}&t=k&z=17",
            }
        )
    # also positive mounds (mound villages)
    mound = (lrm >= 0.7) & interior
    mound = ndi.binary_opening(mound, iterations=1)
    mlab, mn = ndi.label(mound)
    for i, slc in enumerate(ndi.find_objects(mlab), start=1):
        if slc is None:
            continue
        comp = mlab[slc] == i
        area = float(comp.sum() * res_m * res_m)
        if area < 200 or area > 25000:
            continue
        h = float(np.nanpercentile(lrm[slc][comp], 90))
        if h < 1.0 or h > 8.0:
            continue
        ys, xs = np.where(comp)
        cy = slc[0].start + ys.mean()
        cx = slc[1].start + xs.mean()
        x, y = rio_xy(transform, cy, cx)
        lon, lat = rio_transform(crs, "EPSG:4326", [x], [y])
        lat, lon = float(lat[0]), float(lon[0])
        gname, gdkm, _, _ = nearest_geo(lat, lon)
        # only keep if clustered: count neighbors later
        hits.append(
            {
                "lat": lat,
                "lon": lon,
                "score": float(5.0 * h + area / 2000.0),
                "depth_m": 0.0,
                "height_m": h,
                "diameter_m": math.sqrt(area / math.pi) * 2,
                "circ": 0.0,
                "area_m2": area,
                "rect_fill": 0.0,
                "source": source,
                "dataset": dem_path.name,
                "type": "mound",
                "nearest_geoglyph": gname,
                "dist_known_km": float(gdkm),
                "maps": f"https://www.google.com/maps?q={lat:.6f},{lon:.6f}&t=k&z=17",
            }
        )
    # preview
    v = np.nan_to_num(lrm, nan=0.0)
    lo, hi = np.percentile(v[valid], [3, 97])
    t = np.clip((v - lo) / max(hi - lo, 1e-6), 0, 1)
    rgb = np.zeros(v.shape + (3,), dtype=np.uint8)
    rgb[..., 0] = (30 + 200 * t).astype(np.uint8)
    rgb[..., 1] = (50 + 140 * (1 - np.abs(t - 0.5) * 2)).astype(np.uint8)
    rgb[..., 2] = (40 + 200 * (1 - t)).astype(np.uint8)
    rgb[~valid] = 10
    step = max(1, int(round(dem.shape[1] / 1800)))
    im = Image.fromarray(rgb[::step, ::step])
    draw = ImageDraw.Draw(im)
    for j, h in enumerate(sorted(hits, key=lambda z: -z["score"])[:15], 1):
        # approximate pixel from latlon via reverse is heavy; skip if no rowcol stored
        pass
    im.save(PREV / f"{dem_path.stem}_lrm.png")
    hits.sort(key=lambda h: h["score"], reverse=True)
    # NMS 60m
    kept = []
    for h in hits:
        if any(haversine_km(h["lat"], h["lon"], k["lat"], k["lon"]) * 1000 < 60 for k in kept):
            continue
        kept.append(h)
    print(f"  hits={len(kept)}", flush=True)
    return kept


def main():
    all_hits = []
    # SFB Jacunda DTMs
    for p in sorted((WORK / "amazon_sfb").rglob("dtm.tif")):
        all_hits.extend(detect_earthworks(p, source="SFB_Jacunda"))
    # Jamari if extracted
    for p in sorted((WORK / "amazon_sfb").rglob("**/dtm.tif")):
        pass  # already covered by rglob above
    # EBA LAZ -> DTM
    for laz in sorted((WORK / "amazon_eba/laz").glob("*.laz")):
        dtm = laz_to_dtm(laz, res_m=1.0)
        all_hits.extend(detect_earthworks(dtm, source="EBA_LAZ"))
    all_hits.sort(key=lambda h: h["score"], reverse=True)
    novel = [h for h in all_hits if h["dist_known_km"] > 2.0]
    known_overlap = [h for h in all_hits if h["dist_known_km"] <= 2.0]
    (OUT / "amazon_all_hits.json").write_text(json.dumps(all_hits[:200], indent=2))
    (OUT / "amazon_novel.json").write_text(json.dumps(novel[:80], indent=2))
    (OUT / "amazon_known_overlap.json").write_text(json.dumps(known_overlap[:40], indent=2))
    print("TOTAL", len(all_hits), "novel>2km", len(novel), "overlap", len(known_overlap))
    print("TOP novel:")
    for i, h in enumerate(novel[:15], 1):
        print(
            f"  {i}. score={h['score']:.1f} diam={h.get('diameter_m',0):.0f}m "
            f"{h['lat']:.5f},{h['lon']:.5f} {h['source']} known={h['dist_known_km']:.1f}km"
        )


if __name__ == "__main__":
    main()
