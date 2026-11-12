#!/usr/bin/env python3
"""Motagua–Polochic UAV LiDAR structure hunter.

Scans OpenTopography Community Dataspace 0.5 m (Quirigua 0.25 m) DTMs for
multi-mound orthogonal / plaza-like complexes and isolated monumental mounds,
using the same LRM + architectural-grammar approach as the Yucatán G-LiHT hunt.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image, ImageDraw, ImageFont
from rasterio.transform import xy as rio_xy
from rasterio.warp import transform as rio_transform
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

WORK = Path("/agent/lidar/work/open_lidar_hunt_v2")
DTM_DIR = WORK / "motagua/dtm"
OUT = WORK / "out/motagua"
PREV = WORK / "previews/motagua"
OUT.mkdir(parents=True, exist_ok=True)
PREV.mkdir(parents=True, exist_ok=True)

# Known archaeological landmarks for cross-check (approx centers)
KNOWN = [
    {"name": "Quirigua (Museo / Acropolis)", "lat": 15.2694, "lon": -89.0404, "note": "UNESCO Maya site; positive-control footprint"},
    {"name": "Quirigua Great Plaza vicinity", "lat": 15.2715, "lon": -89.0400, "note": "Classic Maya"},
]


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(min(1.0, a)))


def norm_gaussian_lrm(dem, valid, sigma_px):
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


def extract_mounds(lrm, slope_deg, interior, res_m):
    base = (lrm >= 0.55) & interior
    base = ndi.binary_opening(base, iterations=1)
    base = ndi.binary_closing(base, iterations=1)
    labeled, nlab = ndi.label(base)
    if nlab == 0:
        return []
    mounds = []
    px_area = res_m * res_m
    for lab_idx, slc in enumerate(ndi.find_objects(labeled), start=1):
        if slc is None:
            continue
        comp = labeled[slc] == lab_idx
        n_px = int(comp.sum())
        area = n_px * px_area
        if area < 35.0 or area > 12000.0:
            continue
        vals = lrm[slc][comp]
        vals = vals[np.isfinite(vals)]
        if len(vals) < 25:
            continue
        h_max = float(np.max(vals))
        h_p90 = float(np.percentile(vals, 90))
        if h_max < 0.80 or h_p90 < 0.60:
            continue
        slp = slope_deg[slc][comp]
        slp_p80 = float(np.nanpercentile(slp, 80)) if np.isfinite(slp).any() else 0.0
        if slp_p80 < 6.0:
            continue
        ys_rel, xs_rel = np.where(comp)
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
        if w_m < 6.0 or w_m > 160.0 or h_m < 4.0:
            continue
        aspect = w_m / h_m
        if aspect > 6.0:
            continue
        rect = min(1.0, area / max(w_m * h_m, 1.0))
        theta = float(math.atan2(evecs[1, 0], evecs[0, 0]) % math.pi)
        # right-angle score
        pad = (
            slice(max(0, slc[0].start - 3), min(lrm.shape[0], slc[0].stop + 3)),
            slice(max(0, slc[1].start - 3), min(lrm.shape[1], slc[1].stop + 3)),
        )
        c_pad = labeled[pad] == lab_idx
        eroded = ndi.binary_erosion(c_pad, iterations=1)
        ey, ex = np.where(c_pad & ~eroded)
        ra = 0.0
        if len(ex) >= 14:
            dist_f = ndi.distance_transform_edt(~c_pad)
            gy, gx = np.gradient(dist_f.astype(np.float32))
            ang = np.mod(np.arctan2(gy[ey, ex], gx[ey, ex]), np.pi)
            hist, bins = np.histogram(ang, bins=18, range=(0, np.pi))
            top2 = np.argsort(hist)[-2:]
            if hist[top2].sum() >= 8:
                a0 = (bins[top2[0]] + bins[top2[0] + 1]) / 2
                a1 = (bins[top2[1]] + bins[top2[1] + 1]) / 2
                diff = abs(a0 - a1)
                diff = min(diff, np.pi - diff)
                ra = float(max(0.0, 1.0 - abs(diff - (np.pi / 2)) / (np.pi / 4)))
                ra *= 0.5 + 0.5 * float(hist[top2].sum() / max(hist.sum(), 1))
        mounds.append(
            {
                "cy": float(mean[1]),
                "cx": float(mean[0]),
                "area_m2": area,
                "length_m": w_m,
                "width_m": h_m,
                "aspect": aspect,
                "h_max": h_max,
                "h_p90": h_p90,
                "slp_p80": slp_p80,
                "rect": rect,
                "ra": ra,
                "theta": theta,
                "ys": ys,
                "xs": xs,
            }
        )
    return mounds


def score_complexes(mounds, lrm, slope_deg, interior, res_m):
    if len(mounds) < 2:
        return []
    centers = np.array([[m["cy"], m["cx"]] for m in mounds])
    tree = cKDTree(centers * res_m)
    complexes = []
    used = np.zeros(len(mounds), dtype=bool)
    order = sorted(range(len(mounds)), key=lambda i: mounds[i]["h_max"], reverse=True)
    for i in order:
        if used[i]:
            continue
        if mounds[i]["h_max"] < 1.05:
            break
        members = tree.query_ball_point(centers[i] * res_m, r=95.0)
        if len(members) < 2:
            continue
        # require at least 2 solid mounds
        mem = [mounds[j] for j in members if mounds[j]["h_max"] >= 0.85]
        if len(mem) < 2:
            continue
        for j in members:
            used[j] = True
        hs = sorted([m["h_max"] for m in mem], reverse=True)
        rect_mean = float(np.mean([m["rect"] for m in mem]))
        ra_mean = float(np.mean([m["ra"] for m in mem]))
        # quadrant diversity
        cy = float(np.mean([m["cy"] for m in mem]))
        cx = float(np.mean([m["cx"] for m in mem]))
        quads = set()
        for m in mem:
            quads.add((0 if m["cy"] >= cy else 1, 0 if m["cx"] >= cx else 1))
        # plaza: flat patch near center
        r_px = int(45.0 / res_m)
        y0, y1 = max(0, int(cy) - r_px), min(lrm.shape[0], int(cy) + r_px)
        x0, x1 = max(0, int(cx) - r_px), min(lrm.shape[1], int(cx) + r_px)
        sub_l = lrm[y0:y1, x0:x1]
        sub_s = slope_deg[y0:y1, x0:x1]
        sub_i = interior[y0:y1, x0:x1]
        plaza = sub_i & np.isfinite(sub_l) & (np.abs(sub_l) <= 0.45) & (sub_s <= 5.0)
        plaza_m2 = float(plaza.sum() * res_m * res_m)
        score = (
            8.0 * hs[0]
            + 4.0 * (hs[1] if len(hs) > 1 else 0)
            + 3.5 * len(mem)
            + 6.0 * rect_mean
            + 5.0 * ra_mean
            + 2.5 * len(quads)
            + min(6.0, plaza_m2 / 120.0)
        )
        # reject very elongated road embankments disguised as complexes
        if np.mean([m["aspect"] for m in mem]) > 4.2 and ra_mean < 0.25:
            continue
        if score < 18.0:
            continue
        complexes.append(
            {
                "cy": cy,
                "cx": cx,
                "n_mounds": len(mem),
                "h_max": hs[0],
                "h2": hs[1] if len(hs) > 1 else 0.0,
                "rect_mean": rect_mean,
                "ra_mean": ra_mean,
                "quads": len(quads),
                "plaza_m2": plaza_m2,
                "score": float(score),
                "mounds": [
                    {
                        "h_max": m["h_max"],
                        "area_m2": m["area_m2"],
                        "length_m": m["length_m"],
                        "width_m": m["width_m"],
                        "rect": m["rect"],
                        "ra": m["ra"],
                    }
                    for m in mem
                ],
            }
        )
    # also keep solitary monumental mounds
    for m in mounds:
        if m["h_max"] >= 2.4 and m["area_m2"] >= 180 and m["rect"] >= 0.45:
            complexes.append(
                {
                    "cy": m["cy"],
                    "cx": m["cx"],
                    "n_mounds": 1,
                    "h_max": m["h_max"],
                    "h2": 0.0,
                    "rect_mean": m["rect"],
                    "ra_mean": m["ra"],
                    "quads": 1,
                    "plaza_m2": 0.0,
                    "score": float(10.0 * m["h_max"] + 8.0 * m["rect"] + 4.0 * m["ra"]),
                    "mounds": [
                        {
                            "h_max": m["h_max"],
                            "area_m2": m["area_m2"],
                            "length_m": m["length_m"],
                            "width_m": m["width_m"],
                            "rect": m["rect"],
                            "ra": m["ra"],
                        }
                    ],
                    "solitary": True,
                }
            )
    # NMS
    complexes.sort(key=lambda c: c["score"], reverse=True)
    kept = []
    for c in complexes:
        if any(
            math.hypot((c["cy"] - k["cy"]) * res_m, (c["cx"] - k["cx"]) * res_m) < 70.0
            for k in kept
        ):
            continue
        kept.append(c)
    return kept


def lrm_to_rgb(lrm, valid):
    v = np.nan_to_num(lrm, nan=0.0)
    lo, hi = np.percentile(v[valid], [2, 98]) if valid.any() else (-1, 1)
    t = np.clip((v - lo) / max(hi - lo, 1e-6), 0, 1)
    rgb = np.zeros(v.shape + (3,), dtype=np.uint8)
    # teal-brown diverging-ish
    rgb[..., 0] = (40 + 180 * t).astype(np.uint8)
    rgb[..., 1] = (50 + 140 * (1 - np.abs(t - 0.5) * 2)).astype(np.uint8)
    rgb[..., 2] = (40 + 180 * (1 - t)).astype(np.uint8)
    rgb[~valid] = 20
    return rgb


def render_preview(stem, dem, lrm, valid, transform, crs, cands, res_m):
    rgb = lrm_to_rgb(lrm, valid)
    # downsample for preview
    step = max(1, int(round(2.0 / res_m)))
    small = rgb[::step, ::step]
    im = Image.fromarray(small)
    draw = ImageDraw.Draw(im)
    for i, c in enumerate(cands[:12], start=1):
        x = c["cx"] / step
        y = c["cy"] / step
        r = max(6, 25 / step)
        draw.ellipse([x - r, y - r, x + r, y + r], outline=(255, 255, 0), width=2)
        draw.text((x + r + 2, y - 6), f"{i}:{c['score']:.0f}", fill=(255, 255, 0))
    out = PREV / f"{stem}_lrm_candidates.png"
    im.save(out)
    return str(out)


def process_dtm(path: Path):
    print(f"SCAN {path.name}", flush=True)
    with rasterio.open(path) as ds:
        dem = ds.read(1).astype(np.float32)
        transform = ds.transform
        crs = ds.crs
        nodata = ds.nodata
        res_m = float(abs(ds.res[0]))
    valid = np.isfinite(dem)
    if nodata is not None:
        valid &= dem != nodata
    # COGs sometimes use very large nodata
    if valid.any():
        p1, p99 = np.percentile(dem[valid], [0.5, 99.5])
        valid &= (dem > p1 - 50) & (dem < p99 + 50)
    dem = np.where(valid, dem, np.nan)
    # interior buffer away from edges / voids
    interior = ndi.binary_erosion(valid, iterations=max(3, int(4 / res_m)))
    gy, gx = np.gradient(np.nan_to_num(dem, nan=np.nanmedian(dem[valid]) if valid.any() else 0))
    slope_deg = np.degrees(np.arctan(np.hypot(gx, gy) / res_m))
    slope_deg[~valid] = np.nan
    sigma18 = 18.0 / res_m
    lrm = norm_gaussian_lrm(dem, valid, sigma18)
    mounds = extract_mounds(lrm, slope_deg, interior, res_m)
    complexes = score_complexes(mounds, lrm, slope_deg, interior, res_m)
    results = []
    for c in complexes:
        x, y = rio_xy(transform, c["cy"], c["cx"])
        lon, lat = rio_transform(crs, "EPSG:4326", [x], [y])
        lat, lon = float(lat[0]), float(lon[0])
        nearest = min(KNOWN, key=lambda k: haversine_km(lat, lon, k["lat"], k["lon"]))
        d_km = haversine_km(lat, lon, nearest["lat"], nearest["lon"])
        rec = {
            **{k: v for k, v in c.items() if k != "mounds"},
            "mounds": c["mounds"],
            "lat": lat,
            "lon": lon,
            "dataset": path.name,
            "res_m": res_m,
            "nearest_known": nearest["name"],
            "dist_known_km": d_km,
            "maps": f"https://www.google.com/maps?q={lat:.6f},{lon:.6f}&t=k&z=18",
        }
        results.append(rec)
    preview = render_preview(path.stem, dem, lrm, valid, transform, crs, results, res_m)
    # also save full-tile hillshade-ish overview
    print(f"  mounds={len(mounds)} complexes={len(results)} preview={preview}", flush=True)
    return results


def main():
    all_hits = []
    for p in sorted(DTM_DIR.glob("*.tif")):
        # Skip military base for ethics / sensitivity? Still scan but flag
        hits = process_dtm(p)
        all_hits.extend(hits)
    all_hits.sort(key=lambda h: h["score"], reverse=True)
    out_json = OUT / "motagua_candidates.json"
    out_json.write_text(json.dumps(all_hits, indent=2))
    # screened top list excluding known Quirigua < 0.4 km
    screened = [h for h in all_hits if h["dist_known_km"] > 0.35]
    (OUT / "motagua_screened.json").write_text(json.dumps(screened[:40], indent=2))
    summary = {
        "n_total": len(all_hits),
        "n_screened_away_from_quirigua": len(screened),
        "top10": screened[:10],
        "positive_control_hits_near_quirigua": [h for h in all_hits if h["dist_known_km"] <= 0.35][:5],
    }
    (OUT / "motagua_summary.json").write_text(json.dumps(summary, indent=2))
    print("WROTE", out_json)
    print("TOP screened:")
    for i, h in enumerate(screened[:12], 1):
        print(
            f"  {i}. score={h['score']:.1f} n={h['n_mounds']} h={h['h_max']:.1f}m "
            f"{h['lat']:.5f},{h['lon']:.5f} {h['dataset']} known={h['dist_known_km']:.2f}km"
        )


if __name__ == "__main__":
    main()
