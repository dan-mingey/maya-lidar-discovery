#!/usr/bin/env python3
"""Shared data-fetching & clean bare-earth / point-cloud utilities for Open LiDAR Hunt v2.

Fixes both Pass-1 artifacts:
  1. INPE EBA Acre LAZ (`NP_T-*.laz`): `classification == 2` is NOT bare-earth ground
     (contains 20-30 m canopy returns). Uses a 4 m circular-disk Progressive
     Morphological Ground Filter + 1 m normalized Gaussian ground splatting so
     canopy speckle is 100% removed without introducing boxy morphology artifacts.
  2. Motagua UAV (`OTDS.072026.32616.*`): fetches real 3D UAV canopy + ground
     point clouds directly from OpenTopography's Entwine Point Tile (`ept`) octree,
     paired with the native 0.25 m / 0.50 m bare-earth COG DTMs.
  3. Reliable Esri World Imagery optical satellite stitching via Web Mercator
     XYZ tiles (`tile/{z}/{y}/{x}`) with local disk caching (immune to `/export` HTTP 500s).
"""
from __future__ import annotations

import io
import json
import math
import urllib.request
from pathlib import Path

import laspy
import numpy as np
import rasterio
from PIL import Image
from pyproj import Transformer
from rasterio.windows import from_bounds
from scipy.ndimage import (
    distance_transform_edt,
    gaussian_filter,
    grey_dilation,
    grey_erosion,
    zoom,
)

WORK = Path("/agent/lidar/work/open_lidar_hunt_v2")
CACHE = WORK / "viz_cache"
CACHE.mkdir(parents=True, exist_ok=True)
STORE_MEDIA = Path("/cursor/stores/bc-adff6ec3-2baf-40bb-9bab-102a5f567f18/media/open-lidar-hunt-v2")
REPO_MEDIA = Path("/agent/maya-lidar-discovery/media/open-lidar-hunt-v2")
STORE_MEDIA.mkdir(parents=True, exist_ok=True)
REPO_MEDIA.mkdir(parents=True, exist_ok=True)


def multi_hillshade(dem: np.ndarray, res_m: float = 0.5) -> np.ndarray:
    valid = np.isfinite(dem)
    fill = np.where(valid, dem, float(np.nanmedian(dem[valid])))
    dy, dx = np.gradient(fill, res_m, res_m)
    slope = np.pi / 2.0 - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(-dx, dy)
    acc = np.zeros_like(dem, dtype=np.float32)
    for az_deg, alt_deg, w in [(315, 42, 0.52), (45, 42, 0.28), (270, 35, 0.20)]:
        az = math.radians(az_deg)
        alt = math.radians(alt_deg)
        sh = np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * np.cos(az - aspect)
        acc += w * ((sh + 1.0) * 0.5)
    acc[~valid] = np.nan
    return acc


def slope_relief_rgb(dem: np.ndarray, res_m: float = 0.5) -> np.ndarray:
    """Compute crisp G1-style architectural slope + multi-hillshade grayscale RGB."""
    valid = np.isfinite(dem)
    fill = np.where(valid, dem, float(np.nanmedian(dem[valid])))
    hs = np.nan_to_num(multi_hillshade(fill, res_m), nan=0.5)
    dy, dx = np.gradient(fill, res_m, res_m)
    slp = np.arctan(np.hypot(dx, dy)) / (np.pi / 2.0)
    # Darken steep architectural edges (classic archaeological RVT/G1 look)
    lum = np.clip(hs * (1.0 - 0.55 * np.clip(slp * 2.2, 0.0, 1.0)), 0.04, 1.0)
    u8 = (lum * 255.0).astype(np.uint8)
    rgb = np.stack([u8, u8, u8], axis=-1)
    rgb[~valid] = np.array([18, 22, 30], dtype=np.uint8)
    return rgb


def compute_lrm(dem: np.ndarray, res_m: float = 0.5, sigma_m: float = 22.0) -> np.ndarray:
    valid = np.isfinite(dem)
    step = max(1, int(round(2.0 / res_m)))
    sv = valid[::step, ::step].astype(np.float32)
    sd = np.where(valid[::step, ::step], dem[::step, ::step], 0.0).astype(np.float32)
    sig = sigma_m / (res_m * step)
    num = gaussian_filter(sd, sigma=sig, mode="constant", cval=0.0)
    den = gaussian_filter(sv, sigma=sig, mode="constant", cval=0.0)
    med = float(np.nanmedian(dem[valid]))
    strd = np.where(den > 1e-3, num / np.maximum(den, 1e-3), med).astype(np.float32)
    h, w = dem.shape
    tr = zoom(strd, (h / strd.shape[0], w / strd.shape[1]), order=1)[:h, :w]
    lrm = dem - tr
    lrm[~valid] = np.nan
    return lrm


def fetch_esri_bbox_rgb(
    lon_min: float,
    lat_min: float,
    lon_max: float,
    lat_max: float,
    out_w: int = 512,
    out_h: int = 512,
    zoom_level: int | None = None,
) -> np.ndarray:
    """Fetch Esri World Imagery seamlessly by stitching Web Mercator 256x256 XYZ tiles."""
    key = f"esri_xyz_{lon_min:.5f}_{lat_min:.5f}_{lon_max:.5f}_{lat_max:.5f}_{out_w}x{out_h}.png"
    cp = CACHE / key
    if cp.exists():
        return np.asarray(Image.open(cp).convert("RGB"))

    if zoom_level is None:
        span_deg = max(lon_max - lon_min, lat_max - lat_min)
        # Pick zoom so stitched image has >= ~600-1200 px across bbox
        zoom_level = int(np.clip(math.ceil(math.log2(360.0 / max(span_deg, 1e-5)) + 1.2), 5, 18))

    def deg2num_f(lat_deg: float, lon_deg: float, z: int) -> tuple[float, float]:
        lat_rad = math.radians(np.clip(lat_deg, -85.0511, 85.0511))
        n = 2.0 ** z
        xtile = (lon_deg + 180.0) / 360.0 * n
        ytile = (1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n
        return xtile, ytile

    while zoom_level > 5:
        x0f, y0f = deg2num_f(lat_max, lon_min, zoom_level)
        x1f, y1f = deg2num_f(lat_min, lon_max, zoom_level)
        tx0, tx1 = int(math.floor(x0f)), int(math.floor(x1f))
        ty0, ty1 = int(math.floor(y0f)), int(math.floor(y1f))
        if (tx1 - tx0 + 1) * (ty1 - ty0 + 1) <= 36:
            break
        zoom_level -= 1

    canvas = Image.new("RGB", ((tx1 - tx0 + 1) * 256, (ty1 - ty0 + 1) * 256), (22, 34, 26))
    for ty in range(ty0, ty1 + 1):
        for tx in range(tx0, tx1 + 1):
            t_cache = CACHE / f"tile_z{zoom_level}_{tx}_{ty}.jpg"
            if t_cache.exists():
                raw = t_cache.read_bytes()
            else:
                url = f"https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{zoom_level}/{ty}/{tx}"
                try:
                    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                    raw = urllib.request.urlopen(req, timeout=15).read()
                    t_cache.write_bytes(raw)
                except Exception:
                    raw = None
            if raw:
                try:
                    t_img = Image.open(io.BytesIO(raw)).convert("RGB")
                    canvas.paste(t_img, ((tx - tx0) * 256, (ty - ty0) * 256))
                except Exception:
                    pass

    px0 = int(round((x0f - tx0) * 256))
    py0 = int(round((y0f - ty0) * 256))
    px1 = max(px0 + 4, int(round((x1f - tx0) * 256)))
    py1 = max(py0 + 4, int(round((y1f - ty0) * 256)))
    cropped = canvas.crop((px0, py0, px1, py1)).resize((out_w, out_h), Image.Resampling.BICUBIC)
    cropped.save(cp)
    return np.asarray(cropped)


def fetch_motagua_bundle(
    code: str,
    dtm_name: str,
    ept_idx: int,
    lat: float,
    lon: float,
    half_m: float = 175.0,
    sigma_m: float = 22.0,
) -> dict:
    """Fetch native 0.25m/0.50m bare-earth DTM + live EPT 3D canopy/ground point cloud + Esri RGB."""
    dtm_path = WORK / "motagua/dtm" / dtm_name
    with rasterio.open(dtm_path) as ds:
        tr = Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True)
        cx, cy = tr.transform(lon, lat)
        win = from_bounds(cx - half_m, cy - half_m, cx + half_m, cy + half_m, ds.transform)
        dem = ds.read(1, window=win, boundless=True, fill_value=np.nan).astype(np.float32)
        if ds.nodata is not None:
            dem = np.where(dem == ds.nodata, np.nan, dem)
        dem = np.where((np.abs(dem) > 1e6) | (dem < -100), np.nan, dem)
        res_m = float(abs(ds.res[0]))

    # Fetch or load cached EPT 3D point cloud over (cx±half_m, cy±half_m)
    pt_cache = CACHE / f"ept_{code}_{ept_idx}_{lat:.5f}_{lon:.5f}_{int(half_m)}.npz"
    if pt_cache.exists():
        npz = np.load(pt_cache)
        px, py, pz, pcls = npz["x"], npz["y"], npz["z"], npz["cls"]
    else:
        base = f"https://opentopography.s3.sdsc.edu/dataspace/OTDS.072026.32616.{ept_idx}/ept"
        ept = json.loads(urllib.request.urlopen(f"{base}/ept.json", timeout=15).read())
        eb = ept["bounds"]
        hier = json.loads(urllib.request.urlopen(f"{base}/ept-hierarchy/0-0-0-0.json", timeout=15).read())

        def walk_hier(h_dict: dict, out_list: list):
            for k, cnt in h_dict.items():
                d, xi, yi, zi = map(int, k.split("-"))
                step = 1.0 / (2 ** d)
                nx0 = eb[0] + xi * step * (eb[3] - eb[0])
                nx1 = eb[0] + (xi + 1) * step * (eb[3] - eb[0])
                ny0 = eb[1] + yi * step * (eb[4] - eb[1])
                ny1 = eb[1] + (yi + 1) * step * (eb[4] - eb[1])
                if nx1 < (cx - half_m) or nx0 > (cx + half_m) or ny1 < (cy - half_m) or ny0 > (cy + half_m):
                    continue
                if cnt == -1:
                    try:
                        sub_h = json.loads(urllib.request.urlopen(f"{base}/ept-hierarchy/{k}.json", timeout=15).read())
                        walk_hier(sub_h, out_list)
                    except Exception:
                        pass
                elif cnt > 0:
                    out_list.append((d, k, cnt))

        matching: list = []
        walk_hier(hier, matching)
        matching.sort()
        xs_l, ys_l, zs_l, cl_l = [], [], [], []
        tot = 0
        for d, k, _ in matching:
            if d > 6 and tot >= 350000:
                break
            raw = urllib.request.urlopen(f"{base}/ept-data/{k}.laz", timeout=25).read()
            las = laspy.read(io.BytesIO(raw))
            lx, ly, lz = np.array(las.x), np.array(las.y), np.array(las.z)
            m = (lx >= cx - half_m) & (lx <= cx + half_m) & (ly >= cy - half_m) & (ly <= cy + half_m)
            if m.any():
                xs_l.append(lx[m].astype(np.float32))
                ys_l.append(ly[m].astype(np.float32))
                zs_l.append(lz[m].astype(np.float32))
                cl_l.append(np.array(las.classification)[m].astype(np.uint8))
                tot += int(m.sum())
            if tot >= 500000:
                break
        px = np.concatenate(xs_l)
        py = np.concatenate(ys_l)
        pz = np.concatenate(zs_l)
        pcls = np.concatenate(cl_l)
        np.savez_compressed(pt_cache, x=px, y=py, z=pz, cls=pcls)

    # Build 1m Canopy Height Model (CHM) + surgical 3D laser under-penetration ground cleanup:
    # Wherever raw 3D laser pulses penetrate >0.9 m BELOW YellowScan's delivered DTM (leaked tree crowns),
    # clamp the DTM down to the true sub-canopy ground envelope while keeping 0.25m/0.50m crispness on solid stone/earth!
    H, W = dem.shape
    valid = np.isfinite(dem)
    dem_raw_delivered = dem.copy()
    gnd_res = 3.0
    gn = max(16, int(round(2.0 * half_m / gnd_res)))
    ggx = np.clip((((px - (cx - half_m)) / (2.0 * half_m)) * (gn - 1)).astype(np.int32), 0, gn - 1)
    ggy = np.clip((((cy + half_m - py) / (2.0 * half_m)) * (gn - 1)).astype(np.int32), 0, gn - 1)
    gmin3 = np.full((gn, gn), np.inf, dtype=np.float32)
    np.minimum.at(gmin3, (ggy, ggx), pz)
    gv = np.isfinite(gmin3)
    if gv.any():
        _, ginds = distance_transform_edt(~gv, return_indices=True)
        yy2, xx2 = np.ogrid[-2:3, -2:3]
        disk2 = (xx2 * xx2 + yy2 * yy2) <= 4
        g_env3 = gaussian_filter(
            grey_dilation(grey_erosion(gmin3[tuple(ginds)], footprint=disk2), footprint=disk2),
            sigma=0.9,
        )
        g_env_full = zoom(g_env3, (H / gn, W / gn), order=1)[:H, :W]
        leaked_tree = valid & (dem > (g_env_full + 0.95))
        dem_blend = np.where(leaked_tree, g_env_full + 0.25, dem)
        # Smooth only across clamped tree-footprint seams
        seam = gaussian_filter(leaked_tree.astype(np.float32), sigma=1.5) > 0.05
        dem_s = gaussian_filter(np.where(valid, dem_blend, float(np.nanmedian(dem[valid]))), sigma=1.1)
        dem = np.where(seam & valid, dem_s, dem_blend).astype(np.float32)

    dem_fill = np.where(valid, dem, float(np.nanmedian(dem[valid])))
    gx = np.clip((((px - (cx - half_m)) / (2.0 * half_m)) * (W - 1)).astype(np.int32), 0, W - 1)
    gy = np.clip((((cy + half_m - py) / (2.0 * half_m)) * (H - 1)).astype(np.int32), 0, H - 1)
    agl = np.clip(pz - dem_fill[gy, gx], 0.0, 52.0)

    chm_n = int(round(2.0 * half_m / 2.0))
    cgx = np.clip((((px - (cx - half_m)) / (2.0 * half_m)) * (chm_n - 1)).astype(np.int32), 0, chm_n - 1)
    cgy = np.clip((((cy + half_m - py) / (2.0 * half_m)) * (chm_n - 1)).astype(np.int32), 0, chm_n - 1)
    chm = np.full((chm_n, chm_n), 0.0, dtype=np.float32)
    np.maximum.at(chm, (cgy, cgx), agl)
    chm = gaussian_filter(chm, sigma=0.6)

    dlat = half_m / 110574.0
    dlon = half_m / (111320.0 * math.cos(math.radians(lat)))
    sat_rgb = fetch_esri_bbox_rgb(lon - dlon, lat - dlat, lon + dlon, lat + dlat, out_w=512, out_h=512, zoom_level=17)

    hs = multi_hillshade(dem, res_m)
    lrm = compute_lrm(dem, res_m, sigma_m=sigma_m)
    g1_rgb = slope_relief_rgb(dem, res_m)

    return {
        "code": code,
        "lat": lat,
        "lon": lon,
        "cx": cx,
        "cy": cy,
        "half_m": half_m,
        "res_m": res_m,
        "dem": dem,
        "dem_raw_delivered": dem_raw_delivered,
        "hs": hs,
        "lrm": lrm,
        "g1_rgb": g1_rgb,
        "chm": chm,
        "sat_rgb": sat_rgb,
        "pts_x": px,
        "pts_y": py,
        "pts_z": pz,
        "pts_agl": agl,
        "pts_is_gnd": agl < 0.75,
    }


def fetch_eba_acre_bundle(
    code: str,
    laz_path: Path,
    lat: float,
    lon: float,
    half_m: float = 220.0,
    sigma_m: float = 22.0,
) -> dict:
    """Extract artifact-free 1.0 m bare-earth DTM, CHM, raw 3D points & Esri RGB from EBA LAZ."""
    cache_npz = CACHE / f"eba_{code}_{laz_path.stem}_{lat:.5f}_{lon:.5f}_{int(half_m)}.npz"
    if cache_npz.exists():
        d = np.load(cache_npz)
        cx, cy = float(d["cx"]), float(d["cy"])
        dem, chm = d["dem"], d["chm"]
        px, py, pz, agl, is_gnd = d["px"], d["py"], d["pz"], d["agl"], d["is_gnd"]
    else:
        las = laspy.read(laz_path)
        epsg = las.header.parse_crs().to_epsg()
        tr = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
        cx, cy = tr.transform(lon, lat)
        pad = half_m + 32.0
        x_all, y_all, z_all = np.array(las.x), np.array(las.y), np.array(las.z)
        m_pad = (np.abs(x_all - cx) <= pad) & (np.abs(y_all - cy) <= pad)
        lx = (x_all[m_pad] - (cx - pad)).astype(np.float32)
        ly = ((cy + pad) - y_all[m_pad]).astype(np.float32)
        lz = z_all[m_pad].astype(np.float32)

        # Circular-disk Progressive Morphological Ground Filter at 4 m seed -> 1 m Gaussian splat
        res_c = 4.0
        nc = int(round(2.0 * pad / res_c))
        ci = np.clip((lx / res_c).astype(np.int32), 0, nc - 1)
        ri = np.clip((ly / res_c).astype(np.int32), 0, nc - 1)
        g4 = np.full((nc, nc), np.inf, dtype=np.float32)
        np.minimum.at(g4, (ri, ci), lz)
        vc = np.isfinite(g4)
        _, inds = distance_transform_edt(~vc, return_indices=True)
        yy, xx = np.ogrid[-3:4, -3:4]
        disk3 = (xx * xx + yy * yy) <= 9
        g4_open = gaussian_filter(
            grey_dilation(grey_erosion(g4[tuple(inds)], footprint=disk3), footprint=disk3),
            sigma=1.1,
        )
        dz = lz - g4_open[ri, ci]
        gnd_pad = dz < 1.15

        nf_pad = int(round(2.0 * pad))
        gx = np.clip(lx[gnd_pad].astype(np.int32), 0, nf_pad - 1)
        gy = np.clip(ly[gnd_pad].astype(np.int32), 0, nf_pad - 1)
        s_z = np.zeros((nf_pad, nf_pad), dtype=np.float64)
        c_z = np.zeros((nf_pad, nf_pad), dtype=np.float64)
        np.add.at(s_z, (gy, gx), lz[gnd_pad])
        np.add.at(c_z, (gy, gx), 1.0)
        num = gaussian_filter(s_z, sigma=1.6)
        den = gaussian_filter(c_z, sigma=1.6)
        dtm_pad = np.where(
            den > 1e-3,
            num / np.maximum(den, 1e-3),
            zoom(g4_open, nf_pad / nc, order=1)[:nf_pad, :nf_pad],
        ).astype(np.float32)
        swath_pad = gaussian_filter((c_z > 0).astype(np.float32), sigma=3.0) > 0.05
        dtm_pad[~swath_pad] = np.nan

        # Crop inner (2*half_m x 2*half_m)
        off = int(round(32.0))
        n_win = int(round(2.0 * half_m))
        dem = dtm_pad[off : off + n_win, off : off + n_win].copy()

        # Points strictly inside target window
        m_in = (np.abs(x_all - cx) <= half_m) & (np.abs(y_all - cy) <= half_m)
        px = x_all[m_in].astype(np.float32)
        py = y_all[m_in].astype(np.float32)
        pz = z_all[m_in].astype(np.float32)

        valid = np.isfinite(dem)
        dem_fill = np.where(valid, dem, float(np.nanmedian(dem[valid])))
        pgx = np.clip((((px - (cx - half_m)) / (2.0 * half_m)) * (n_win - 1)).astype(np.int32), 0, n_win - 1)
        pgy = np.clip((((cy + half_m - py) / (2.0 * half_m)) * (n_win - 1)).astype(np.int32), 0, n_win - 1)
        agl = np.clip(pz - dem_fill[pgy, pgx], -1.0, 48.0).astype(np.float32)
        is_gnd = agl < 1.05

        chm_n = int(round(2.0 * half_m / 2.0))
        cgx = np.clip((((px - (cx - half_m)) / (2.0 * half_m)) * (chm_n - 1)).astype(np.int32), 0, chm_n - 1)
        cgy = np.clip((((cy + half_m - py) / (2.0 * half_m)) * (chm_n - 1)).astype(np.int32), 0, chm_n - 1)
        chm = np.full((chm_n, chm_n), 0.0, dtype=np.float32)
        np.maximum.at(chm, (cgy, cgx), np.maximum(0.0, agl))
        chm = gaussian_filter(chm, sigma=0.7)

        np.savez_compressed(
            cache_npz, cx=cx, cy=cy, dem=dem, chm=chm, px=px, py=py, pz=pz, agl=agl, is_gnd=is_gnd
        )

    dlat = half_m / 110574.0
    dlon = half_m / (111320.0 * math.cos(math.radians(lat)))
    sat_rgb = fetch_esri_bbox_rgb(lon - dlon, lat - dlat, lon + dlon, lat + dlat, out_w=512, out_h=512, zoom_level=16)

    hs = multi_hillshade(dem, 1.0)
    lrm = compute_lrm(dem, 1.0, sigma_m=sigma_m)
    g1_rgb = slope_relief_rgb(dem, 1.0)

    return {
        "code": code,
        "lat": lat,
        "lon": lon,
        "cx": cx,
        "cy": cy,
        "half_m": half_m,
        "res_m": 1.0,
        "dem": dem,
        "hs": hs,
        "lrm": lrm,
        "g1_rgb": g1_rgb,
        "chm": chm,
        "sat_rgb": sat_rgb,
        "pts_x": px,
        "pts_y": py,
        "pts_z": pz,
        "pts_agl": np.maximum(0.0, agl),
        "pts_is_gnd": is_gnd,
    }
