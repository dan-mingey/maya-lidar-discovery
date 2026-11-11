#!/usr/bin/env python3
"""Part 1: Exact footprint intersection of all 456 G-LiHT transects against
the 5,223-site EAAMS catalog + INAH SIC directories, plus full cross-check
of the original 518 DEM hits (and top-300 hits.json) and Candidates A-E.
"""
from __future__ import annotations

import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import rasterio
from rasterio.features import shapes as rio_shapes
from rasterio.transform import rowcol, xy as rio_xy
from rasterio.warp import transform as rio_transform, transform_geom
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from shapely.geometry import Point, shape
from shapely.ops import unary_union

ROOT = Path("/agent/lidar")
WORK = ROOT / "work/carbon_raw_hunt"
WORK.mkdir(parents=True, exist_ok=True)


def haversine_km(lat1, lon1, lat2, lon2) -> float:
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(min(1.0, a)))


def latlon_to_xyz_km(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """Convert WGS84 lat/lon arrays to 3D chord coordinates (km) for exact KD-tree queries."""
    r = 6371.0088
    lat_r = np.radians(lats)
    lon_r = np.radians(lons)
    cos_lat = np.cos(lat_r)
    x = r * cos_lat * np.cos(lon_r)
    y = r * cos_lat * np.sin(lon_r)
    z = r * np.sin(lat_r)
    return np.column_stack([x, y, z])


def load_eaams() -> list[dict]:
    gj = json.loads((ROOT / "work/ce_deep/eaams/maya_sites.geojson").read_text())
    sites = []
    for idx, f in enumerate(gj["features"]):
        coords = f["geometry"]["coordinates"]
        lon, lat = float(coords[0]), float(coords[1])
        props = f.get("properties") or {}
        sites.append(
            {
                "id": idx,
                "name": props.get("name") or f"Unnamed_{idx}",
                "country": props.get("country") or "",
                "state": props.get("state") or "",
                "rank": props.get("rank"),
                "notes": props.get("notes") or "",
                "reference": props.get("reference") or "",
                "lat": lat,
                "lon": lon,
            }
        )
    return sites


def load_sic() -> list[dict]:
    sic_sites = []
    for p in [
        ROOT / "work/ce_deep/sic_campeche.csv",
        ROOT / "work/ce_deep/sic_qroo.csv",
    ]:
        if not p.exists():
            continue
        raw_bytes = p.read_bytes()
        try:
            txt = raw_bytes.decode("utf-8")
        except UnicodeDecodeError:
            txt = raw_bytes.decode("latin-1")
        reader = csv.DictReader(txt.splitlines())
        for row in reader:
            try:
                lat = float(row["gmaps_latitud"])
                lon = float(row["gmaps_longitud"])
            except Exception:
                continue
            sic_sites.append(
                {
                    "id": row.get("zona_arqueologica_id"),
                    "name": (row.get("zona_arqueologica_nombre") or "").strip(),
                    "state": (row.get("nom_ent") or "").strip(),
                    "municipality": (row.get("nom_mun") or "").strip(),
                    "lat": lat,
                    "lon": lon,
                    "url": row.get("pagina_web") or "",
                }
            )
    return sic_sites


def extract_tile_footprint_and_original_hits(args: tuple[str, str, bool]) -> dict:
    source, path_s, extract_footprint = args
    path = Path(path_s)
    ds = rasterio.open(path)
    dem = ds.read(1).astype(np.float32)
    nodata = ds.nodata
    if nodata is not None:
        dem = np.where(dem == nodata, np.nan, dem)
    dem = np.where((np.abs(dem) > 1e6) | (dem < -100), np.nan, dem)
    valid = np.isfinite(dem)
    valid_px = int(valid.sum())
    res_x, res_y = abs(ds.res[0]), abs(ds.res[1])
    area_km2 = valid_px * res_x * res_y / 1e6

    poly_wgs84_geojson = None
    bbox_wgs84 = None
    if extract_footprint and valid_px > 100:
        # Downsample valid mask by 10x (~10m) and close tiny scalloped edge pinholes
        step = 10
        sub_valid = ndi.binary_closing(valid[::step, ::step], iterations=2).astype(np.uint8)
        sub_transform = ds.transform * ds.transform.scale(
            dem.shape[1] / sub_valid.shape[1], dem.shape[0] / sub_valid.shape[0]
        )
        geoms = [
            shape(g)
            for g, val in rio_shapes(sub_valid, mask=(sub_valid == 1), transform=sub_transform)
            if val == 1
        ]
        if geoms:
            merged = unary_union(geoms).simplify(15.0, preserve_topology=True)
            wgs_geom = transform_geom(ds.crs, "EPSG:4326", merged.__geo_interface__)
            poly_wgs84_geojson = wgs_geom
            sh = shape(wgs_geom)
            bbox_wgs84 = list(sh.bounds)  # [minx, miny, maxx, maxy]

    # Also reproduce exact scan_structures.py hits so we get ALL 518 unique DEM hits (not truncated to 300)
    orig_hits = []
    if valid_px >= 1500:
        fill = dem.copy()
        ind = ndi.distance_transform_edt(~valid, return_distances=False, return_indices=True)
        fill = dem[tuple(ind)]
        sigma_px = max(20.0 / max(res_x, 0.5), 3.0)
        trend = ndi.gaussian_filter(fill, sigma=sigma_px, mode="nearest")
        lrm = dem - trend
        lrm[~valid] = np.nan
        finite = lrm[valid]
        if finite.size >= 1000:
            thr = min(2.2, max(0.55, float(np.percentile(finite, 92))))
            mask = (lrm >= thr) & valid
            mask = ndi.binary_opening(mask, iterations=1)
            mask = ndi.binary_closing(mask, iterations=2)
            labeled, nlab = ndi.label(mask)
            if nlab > 0:
                sizes = ndi.sum(mask, labeled, index=np.arange(1, nlab + 1))
                order = np.argsort(sizes)[::-1][:20]
                for idx in order:
                    lab = int(idx) + 1
                    area = float(sizes[idx] * res_x * res_y)
                    if area < 60 or area > 25000:
                        continue
                    comp = labeled == lab
                    ys, xs = np.where(comp)
                    if len(xs) < 40:
                        continue
                    vals = lrm[ys, xs]
                    vals = vals[np.isfinite(vals)]
                    if len(vals) < 40:
                        continue
                    height = float(np.percentile(vals, 90))
                    if height < 0.45:
                        continue
                    thr70 = np.percentile(vals, 70)
                    top = dem[ys, xs][lrm[ys, xs] >= thr70]
                    flat = float(1.0 / (1.0 + np.nanstd(top))) if len(top) > 5 else 0.0

                    # Rectangularity & right-angle score identical to scan_structures.py
                    pts = np.column_stack([xs.astype(np.float64), ys.astype(np.float64)])
                    mean = pts.mean(axis=0)
                    X = pts - mean
                    cov = X.T @ X / max(len(X) - 1, 1)
                    evals, evecs = np.linalg.eigh(cov)
                    evecs = evecs[:, np.argsort(evals)[::-1]]
                    rot = X @ evecs
                    w = max(rot[:, 0].max() - rot[:, 0].min(), 1e-6)
                    h = max(rot[:, 1].max() - rot[:, 1].min(), 1e-6)
                    rect_raw = float(len(xs) / (w * h))
                    aspect = max(w / h, h / w)
                    aspect_pen = 1.0 if aspect < 3.5 else max(0.0, 1.0 - (aspect - 3.5) / 4.0)
                    rect = rect_raw * aspect_pen

                    eroded = ndi.binary_erosion(comp, iterations=1)
                    ey_full, ex_full = np.where(comp & ~eroded)
                    ra = 0.0
                    if len(ex_full) >= 20:
                        by0 = max(0, int(ys.min()) - 6)
                        by1 = min(dem.shape[0], int(ys.max()) + 7)
                        bx0 = max(0, int(xs.min()) - 6)
                        bx1 = min(dem.shape[1], int(xs.max()) + 7)
                        comp_sub = comp[by0:by1, bx0:bx1]
                        ey = ey_full - by0
                        ex = ex_full - bx0
                        dist_f = ndi.distance_transform_edt(~comp_sub)
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
                    diam_px = float(math.sqrt(w * w + h * h))
                    diam_m = diam_px * (res_x + res_y) / 2
                    if diam_m < 8 or diam_m > 160:
                        scale = 0.15
                    elif 15 <= diam_m <= 90:
                        scale = 1.0
                    else:
                        scale = 0.6

                    cy, cx = int(ys.mean()), int(xs.mean())
                    r_w = max(int(diam_px * 0.9), 8)
                    y0, y1 = max(0, cy - r_w * 2), min(dem.shape[0], cy + r_w * 2)
                    x0, x1 = max(0, cx - r_w * 2), min(dem.shape[1], cx + r_w * 2)
                    win = lrm[y0:y1, x0:x1]
                    plaza = 0.0
                    if np.isfinite(win).sum() > 50:
                        low = (win < 0.25) & np.isfinite(win)
                        yy, xx = np.ogrid[y0:y1, x0:x1]
                        outside = (yy - cy) ** 2 + (xx - cx) ** 2 > (diam_px * 0.35) ** 2
                        plaza = float(min(1.0, (low & outside).sum() / max(200, r_w * r_w)))
                    score = (
                        2.2 * min(height, 8.0)
                        + 3.0 * rect
                        + 2.5 * ra
                        + 1.5 * flat
                        + 1.2 * plaza
                        + 0.002 * min(area, 4000)
                    ) * scale
                    if score >= 4.5:
                        ux, uy = rio_xy(ds.transform, cy, cx)
                        lons, lats = rio_transform(ds.crs, "EPSG:4326", [ux], [uy])
                        orig_hits.append(
                            {
                                "score": float(score),
                                "path": str(path),
                                "source": source,
                                "lat": float(lats[0]),
                                "lon": float(lons[0]),
                                "height_m": height,
                                "diameter_m": diam_m,
                                "area_m2": area,
                                "rectangularity": rect,
                                "right_angle_score": ra,
                            }
                        )

    crs_str = str(ds.crs)
    ds.close()
    return {
        "source": source,
        "path": str(path),
        "name": path.name,
        "crs": crs_str,
        "valid_px": valid_px,
        "area_km2": area_km2,
        "bbox_wgs84": bbox_wgs84,
        "poly_wgs84": poly_wgs84_geojson,
        "orig_hits": orig_hits,
    }


def nms_hits(hits: list[dict], dist_m: float = 90.0) -> list[dict]:
    hits = sorted(hits, key=lambda h: h["score"], reverse=True)
    kept: list[dict] = []
    for h in hits:
        if not (math.isfinite(h["lat"]) and math.isfinite(h["lon"])):
            continue
        ok = True
        for k in kept:
            if haversine_km(h["lat"], h["lon"], k["lat"], k["lon"]) * 1000.0 < dist_m:
                ok = False
                break
        if ok:
            kept.append(h)
    return kept


def inspect_known_site_in_dem(dem_path: str, site_lat: float, site_lon: float, radius_m: float = 400.0) -> dict:
    """Measure local relief & prominent mounds inside a DEM around a known site coordinate."""
    ds = rasterio.open(dem_path)
    xs, ys = rio_transform("EPSG:4326", ds.crs, [site_lon], [site_lat])
    sx, sy = xs[0], ys[0]
    r_px = int(radius_m / abs(ds.res[0]))
    row, col = rowcol(ds.transform, sx, sy)
    direct_inside = False
    if 0 <= row < ds.height and 0 <= col < ds.width:
        val = ds.read(1, window=((row, row + 1), (col, col + 1)))[0, 0]
        if ds.nodata is None or (val != ds.nodata and abs(val) < 1e6 and val > -100):
            direct_inside = True

    r0, r1 = max(0, row - r_px), min(ds.height, row + r_px)
    c0, c1 = max(0, col - r_px), min(ds.width, col + r_px)
    if r1 <= r0 or c1 <= c0:
        ds.close()
        return {"direct_inside_1m": direct_inside, "valid_window_px": 0}
    sub = ds.read(1, window=((r0, r1), (c0, c1))).astype(np.float32)
    if ds.nodata is not None:
        sub = np.where(sub == ds.nodata, np.nan, sub)
    sub = np.where((np.abs(sub) > 1e6) | (sub < -100), np.nan, sub)
    valid = np.isfinite(sub)
    vpx = int(valid.sum())
    if vpx < 200:
        ds.close()
        return {"direct_inside_1m": direct_inside, "valid_window_px": vpx}

    ind = ndi.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    fill = sub[tuple(ind)]
    lrm = sub - ndi.gaussian_filter(fill, sigma=25.0, mode="nearest")
    lrm[~valid] = np.nan
    max_lrm = float(np.nanmax(lrm))
    p99_lrm = float(np.nanpercentile(lrm[valid], 99))
    mounds_gt_2m = int(ndi.label((lrm >= 2.0) & valid)[1])
    mounds_gt_4m = int(ndi.label((lrm >= 4.0) & valid)[1])

    # Peak LRM location inside window
    mr, mc = np.unravel_index(np.nanargmax(lrm), lrm.shape)
    px, py = rio_xy(ds.transform, r0 + mr, c0 + mc)
    plon, plat = rio_transform(ds.crs, "EPSG:4326", [px], [py])
    ds.close()
    return {
        "direct_inside_1m": direct_inside,
        "valid_window_px": vpx,
        "max_lrm_m": round(max_lrm, 2),
        "p99_lrm_m": round(p99_lrm, 2),
        "mounds_gt_2m": mounds_gt_2m,
        "mounds_gt_4m": mounds_gt_4m,
        "peak_lat": round(float(plat[0]), 6),
        "peak_lon": round(float(plon[0]), 6),
    }


def main() -> None:
    eaams = load_eaams()
    sic = load_sic()
    print(f"Loaded {len(eaams)} EAAMS sites and {len(sic)} INAH SIC sites (Campeche + Q.Roo)")

    tiles: list[tuple[str, str, bool]] = []
    for label, root, fp in [
        ("ot_east", ROOT / "ot/east/dtm", True),
        ("ot_west", ROOT / "ot/west/dtm", True),
        ("nasa", ROOT / "nasa/extracted", False),
    ]:
        for p in sorted(root.rglob("*.tif")):
            tiles.append((label, str(p), fp))
    tiles.sort(key=lambda x: 0 if x[0].startswith("ot_") else 1)
    print(f"Scanning {len(tiles)} raster files (456 OT 1m + 66 NASA)...")

    ot_footprints = []
    raw_orig_hits = []
    workers = min(os.cpu_count() or 4, 14)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        # Preserve deterministic tile order for NMS tie-breaking if needed
        results = list(ex.map(extract_tile_footprint_and_original_hits, tiles))

    for r in results:
        raw_orig_hits.extend(r["orig_hits"])
        if r["source"].startswith("ot_"):
            ot_footprints.append(
                {
                    "source": r["source"],
                    "path": r["path"],
                    "name": r["name"],
                    "crs": r["crs"],
                    "valid_px": r["valid_px"],
                    "area_km2": r["area_km2"],
                    "bbox_wgs84": r["bbox_wgs84"],
                    "poly_wgs84": r["poly_wgs84"],
                }
            )

    total_ot_area = sum(f["area_km2"] for f in ot_footprints)
    print(f"456 OT transects valid bare-earth area sum: {total_ot_area:.2f} km^2")

    # Build Shapely geometries for all 456 OT transects
    strip_shapes = []
    for f in ot_footprints:
        if f["poly_wgs84"] is not None:
            sh = shape(f["poly_wgs84"])
            strip_shapes.append((f, sh))

    # Intersect EAAMS sites against exact valid G-LiHT swath polygons
    # Approximate degree prefilter: 3 km ~ 0.03 deg, so prefilter bbox padded by 0.04 deg
    eaams_intersections = []
    for s in eaams:
        slat, slon = s["lat"], s["lon"]
        pt = Point(slon, slat)
        best_dist_m = 1e9
        best_strips = []
        bbox_hits = []
        for f, sh in strip_shapes:
            minx, miny, maxx, maxy = f["bbox_wgs84"]
            if minx <= slon <= maxx and miny <= slat <= maxy:
                bbox_hits.append(f["name"])
            if not (minx - 0.04 <= slon <= maxx + 0.04 and miny - 0.04 <= slat <= maxy + 0.04):
                continue
            if sh.contains(pt):
                dist_m = 0.0
            else:
                # Exact nearest point on polygon boundary -> haversine meters
                np_pt = sh.exterior.interpolate(sh.exterior.project(pt)) if sh.geom_type == "Polygon" else None
                if np_pt is None:
                    # MultiPolygon: check all parts
                    d_min = 1e9
                    for part in sh.geoms:
                        cand_pt = part.exterior.interpolate(part.exterior.project(pt))
                        d_m = haversine_km(slat, slon, cand_pt.y, cand_pt.x) * 1000.0
                        if part.contains(pt):
                            d_m = 0.0
                        if d_m < d_min:
                            d_min = d_m
                    dist_m = d_min
                else:
                    dist_m = haversine_km(slat, slon, np_pt.y, np_pt.x) * 1000.0
            if dist_m <= 3000.0:
                best_strips.append({"strip": f["name"], "path": f["path"], "dist_m": round(dist_m, 1)})
            if dist_m < best_dist_m:
                best_dist_m = dist_m

        if best_dist_m <= 3000.0 or bbox_hits:
            best_strips.sort(key=lambda x: x["dist_m"])
            rec = {
                **s,
                "min_swath_dist_m": round(best_dist_m, 1),
                "in_raster_bbox": len(bbox_hits) > 0,
                "bbox_strips": bbox_hits[:5],
                "nearby_swaths": best_strips[:5],
            }
            if best_dist_m <= 400.0 and best_strips:
                rec["dem_check"] = inspect_known_site_in_dem(
                    best_strips[0]["path"], slat, slon, radius_m=450.0
                )
            eaams_intersections.append(rec)

    eaams_intersections.sort(key=lambda x: x["min_swath_dist_m"])

    # Same for INAH SIC
    sic_intersections = []
    for s in sic:
        slat, slon = s["lat"], s["lon"]
        pt = Point(slon, slat)
        best_dist_m = 1e9
        best_strip = None
        for f, sh in strip_shapes:
            minx, miny, maxx, maxy = f["bbox_wgs84"]
            if not (minx - 0.06 <= slon <= maxx + 0.06 and miny - 0.06 <= slat <= maxy + 0.06):
                continue
            if sh.contains(pt):
                dist_m = 0.0
            else:
                parts = [sh] if sh.geom_type == "Polygon" else list(sh.geoms)
                dist_m = min(
                    0.0
                    if part.contains(pt)
                    else haversine_km(
                        slat,
                        slon,
                        part.exterior.interpolate(part.exterior.project(pt)).y,
                        part.exterior.interpolate(part.exterior.project(pt)).x,
                    )
                    * 1000.0
                    for part in parts
                )
            if dist_m < best_dist_m:
                best_dist_m = dist_m
                best_strip = f
        if best_dist_m <= 5000.0 and best_strip is not None:
            rec = {
                **s,
                "min_swath_dist_m": round(best_dist_m, 1),
                "nearest_strip": best_strip["name"],
                "nearest_strip_path": best_strip["path"],
            }
            if best_dist_m <= 500.0:
                rec["dem_check"] = inspect_known_site_in_dem(
                    best_strip["path"], slat, slon, radius_m=500.0
                )
            sic_intersections.append(rec)
    sic_intersections.sort(key=lambda x: x["min_swath_dist_m"])

    # Cross-check ALL 518 unique DEM hits + saved hits.json (300) + Candidates A-E against 5,223 EAAMS & SIC
    eaams_lats = np.array([s["lat"] for s in eaams])
    eaams_lons = np.array([s["lon"] for s in eaams])
    eaams_tree = cKDTree(latlon_to_xyz_km(eaams_lats, eaams_lons))

    sic_lats = np.array([s["lat"] for s in sic])
    sic_lons = np.array([s["lon"] for s in sic])
    sic_tree = cKDTree(latlon_to_xyz_km(sic_lats, sic_lons))

    def query_catalogs(lat: float, lon: float) -> dict:
        q = latlon_to_xyz_km(np.array([lat]), np.array([lon]))[0]
        dists_e, idxs_e = eaams_tree.query(q, k=5)
        # Convert chord distance to exact arc distance (virtually identical at <50km)
        e_list = []
        for idx in idxs_e:
            s = eaams[int(idx)]
            dkm = haversine_km(lat, lon, s["lat"], s["lon"])
            e_list.append(
                {
                    "name": s["name"],
                    "rank": s["rank"],
                    "state": s["state"],
                    "dist_km": round(dkm, 3),
                    "lat": s["lat"],
                    "lon": s["lon"],
                    "reference": s["reference"],
                }
            )
        e_list.sort(key=lambda x: x["dist_km"])

        dists_s, idxs_s = sic_tree.query(q, k=3)
        s_list = []
        for idx in idxs_s:
            s = sic[int(idx)]
            dkm = haversine_km(lat, lon, s["lat"], s["lon"])
            s_list.append({"name": s["name"], "state": s["state"], "dist_km": round(dkm, 3)})
        s_list.sort(key=lambda x: x["dist_km"])
        return {
            "nearest_eaams": e_list[0],
            "top3_eaams": e_list[:3],
            "eaams_within_1km": sum(1 for x in e_list if x["dist_km"] <= 1.0),
            "eaams_within_3km": sum(1 for x in e_list if x["dist_km"] <= 3.0),
            "nearest_sic": s_list[0],
        }

    reproduced_unique_hits = nms_hits(raw_orig_hits, dist_m=90.0)
    saved_hits_300 = json.loads((ROOT / "work/hits.json").read_text())
    g1_hits_200 = json.loads((ROOT / "work/g1_hits.json").read_text())

    def summarize_hit_set(hit_list: list[dict], label: str) -> dict:
        d_eaams = []
        d_sic = []
        enriched = []
        for h in hit_list:
            c = query_catalogs(h["lat"], h["lon"])
            de = c["nearest_eaams"]["dist_km"]
            ds_km = c["nearest_sic"]["dist_km"]
            d_eaams.append(de)
            d_sic.append(ds_km)
            enriched.append({**h, "catalog": c})
        arr = np.array(d_eaams)
        return {
            "label": label,
            "total": len(hit_list),
            "within_0p5km": int((arr <= 0.5).sum()),
            "within_1km": int((arr <= 1.0).sum()),
            "within_2km": int((arr <= 2.0).sum()),
            "within_3km": int((arr <= 3.0).sum()),
            "within_5km": int((arr <= 5.0).sum()),
            "within_10km": int((arr <= 10.0).sum()),
            "gt_5km": int((arr > 5.0).sum()),
            "gt_10km": int((arr > 10.0).sum()),
            "sic_within_3km": int((np.array(d_sic) <= 3.0).sum()),
            "sic_within_5km": int((np.array(d_sic) <= 5.0).sum()),
            "enriched": enriched,
        }

    sum_full = summarize_hit_set(reproduced_unique_hits, "full_reproduced_dem_hits")
    sum_300 = summarize_hit_set(saved_hits_300, "saved_hits_json_300")
    sum_g1 = summarize_hit_set(g1_hits_200, "saved_g1_hits_200")

    candidates_abcde = [
        {"id": "A", "tag": "auto_01", "lat": 19.555263, "lon": -89.125789, "strip": "Centro_NFI_l10s481"},
        {"id": "B", "tag": "auto_02", "lat": 19.702893, "lon": -89.341191, "strip": "Centro_NFI_l10s485"},
        {"id": "C", "tag": "auto_05", "lat": 18.971267, "lon": -89.138857, "strip": "Centro_GLAS_l1s463"},
        {"id": "D", "tag": "visual_calakmul", "lat": 18.215700, "lon": -89.888500, "strip": "South_NFI_l2s452"},
        {"id": "E", "tag": "auto_09", "lat": 18.417362, "lon": -89.061804, "strip": "Centro_GLAS_l1s455"},
    ]
    for cand in candidates_abcde:
        cand["catalog"] = query_catalogs(cand["lat"], cand["lon"])

    # Save footprints separately so Part 2 can also reuse strip metadata
    (WORK / "ot_456_footprints.json").write_text(json.dumps(ot_footprints))

    inside_swath = [x for x in eaams_intersections if x["min_swath_dist_m"] == 0.0]
    within_250m = [x for x in eaams_intersections if x["min_swath_dist_m"] <= 250.0]
    within_1km = [x for x in eaams_intersections if x["min_swath_dist_m"] <= 1000.0]
    within_3km = [x for x in eaams_intersections if x["min_swath_dist_m"] <= 3000.0]
    in_bbox_only = [x for x in eaams_intersections if x["in_raster_bbox"]]

    summary = {
        "ot_transects_count": len(ot_footprints),
        "ot_valid_area_km2": round(total_ot_area, 2),
        "eaams_total_sites": len(eaams),
        "eaams_swath_intersection_counts": {
            "directly_inside_valid_swath_0m": len(inside_swath),
            "within_250m_of_valid_swath": len(within_250m),
            "within_1km_of_valid_swath": len(within_1km),
            "within_3km_of_valid_swath": len(within_3km),
            "inside_axis_aligned_raster_bbox": len(in_bbox_only),
        },
        "eaams_inside_or_within_3km": within_3km,
        "sic_within_5km": sic_intersections,
        "dem_hits_full_summary": {k: v for k, v in sum_full.items() if k != "enriched"},
        "dem_hits_top300_summary": {k: v for k, v in sum_300.items() if k != "enriched"},
        "g1_hits_top200_summary": {k: v for k, v in sum_g1.items() if k != "enriched"},
        "candidates_A_to_E_crosscheck": candidates_abcde,
    }
    (WORK / "part1_intersection_summary.json").write_text(json.dumps(summary, indent=2))
    (WORK / "part1_full_dem_hits_enriched.json").write_text(
        json.dumps(sum_full["enriched"], indent=2)
    )

    print("\n=== PART 1 HEADLINE STATS ===")
    print("EAAMS Swath Intersections:", json.dumps(summary["eaams_swath_intersection_counts"], indent=2))
    print("Full DEM Hits vs EAAMS:", json.dumps(summary["dem_hits_full_summary"], indent=2))
    print("Saved hits.json (300) vs EAAMS:", json.dumps(summary["dem_hits_top300_summary"], indent=2))
    print("\n=== KNOWN SITES DIRECTLY INSIDE OR <= 300m OF VALID LIDAR SWATH ===")
    for x in eaams_intersections:
        if x["min_swath_dist_m"] <= 300.0:
            dc = x.get("dem_check", {})
            print(
                f"  dist={x['min_swath_dist_m']:5.1f}m | rank={x['rank']} | {x['name']} ({x['state']}) "
                f"| strip={x['nearby_swaths'][0]['strip']} | maxLRM={dc.get('max_lrm_m')}m "
                f"mounds>2m={dc.get('mounds_gt_2m')} mounds>4m={dc.get('mounds_gt_4m')}"
            )
    print("\n=== INAH SIC SITES WITHIN 5 KM OF VALID LIDAR SWATH ===")
    for x in sic_intersections:
        print(
            f"  dist={x['min_swath_dist_m']:6.1f}m | {x['name']} ({x['state']}) | strip={x['nearest_strip']}"
        )
    print("\n=== CANDIDATES A-E VS 5,223 EAAMS & INAH SIC ===")
    for c in candidates_abcde:
        ne = c["catalog"]["nearest_eaams"]
        ns = c["catalog"]["nearest_sic"]
        print(
            f"  Cand {c['id']} ({c['strip']}): nearest EAAMS = {ne['name']} (rank {ne['rank']}, {ne['dist_km']} km); "
            f"nearest SIC = {ns['name']} ({ns['dist_km']} km)"
        )


if __name__ == "__main__":
    main()
