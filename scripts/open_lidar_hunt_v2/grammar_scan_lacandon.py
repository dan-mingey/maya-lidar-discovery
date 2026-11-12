#!/usr/bin/env python3
"""Run the exact multi-structure Maya Architectural Grammar Hunter (part2_grammar_hunter)
across ALL surrounding Lacandon G-LiHT strips (l0s472..l0s480 + nearby Chiaps1 strips),
then verify the top ~30 complexes on 0.33m G1 + 0.33m DEM + raw .LAZ ground points + Esri
satellite, and render 3 large high-clarity screening sheets (10 cards per sheet) for visual curation.
"""
from __future__ import annotations

import json
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import laspy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/agent/maya-lidar-discovery/scripts")
sys.path.insert(0, "/agent/maya-lidar-discovery/scripts/open_lidar_hunt_v2")
import part2_grammar_hunter as p2
import raw_pointcloud_tools as rpt
from scan_lacandon_corridor import (
    LAT0,
    LON0,
    OUT_DIR,
    build_literature_and_registry_catalog,
    ensure_laz_only,
    get_chm_path_for_strip,
    haversine_km,
    verify_candidate_0p33m_and_laz,
)


def main() -> None:
    inv = json.loads((OUT_DIR / "lacandon_inventory.json").read_text())
    catalog = build_literature_and_registry_catalog()

    # Include l0s472..l0s480 + all Chiaps1 strips within 40 km of #08!
    pairs = []
    for s in inv["nearby_strips"]:
        sid = s["strip_id"]
        if (
            re.match(r"ChiapsCampeche_NFI_l0s4(7[2-9]|80)$", sid)
            or (s["approx_dist_km"] <= 40.0 and sid.startswith("Chiaps"))
        ):
            chm_p = get_chm_path_for_strip(sid) or ""
            pairs.append((s["path"], chm_p))

    print(f"Running part2_grammar_hunter.scan_single_transect on {len(pairs)} Lacandon strips...", flush=True)
    raw_hits = []
    with ProcessPoolExecutor(max_workers=min(12, os.cpu_count() or 4)) as ex:
        for res in ex.map(p2.scan_single_transect, pairs):
            raw_hits.extend(res)

    print(f"Total grammar complexes detected: {len(raw_hits)}", flush=True)

    # Filter to intact forest canopy (chm_med >= 11m and chm_cleared_frac <= 0.14)
    forest_hits = []
    for h in raw_hits:
        sid = rpt.strip_id_from_name(h["strip"])
        dkm = haversine_km(LAT0, LON0, h["lat"], h["lon"])
        # Exclude #08 itself (within 220m of 16.719521, -91.168757) since #08 has its own dedicated zoom sheet
        if dkm < 0.22:
            continue
        if h["chm_med_m"] >= 11.0 and h["chm_cleared_frac"] <= 0.14:
            h["strip_id"] = sid
            h["dist_from_card08_km"] = round(dkm, 2)
            # Composite rank favoring crisp rectangular mounds, orthogonal courtyard, and low karst hill inflation
            karst_pen = max(0.0, h["basal_acropolis_m"] - 4.5) * 1.6
            h["arch_rank_score"] = round(
                h["grammar_raw"]
                + 3.0 * min(h["n_rect_mounds"], 7)
                + 8.0 * h["mean_ra"]
                + 4.0 * min(1.0, h["plaza_area_m2"] / 600.0)
                - karst_pen,
                3,
            )
            forest_hits.append(h)

    # NMS at 260m
    forest_hits.sort(key=lambda x: x["arch_rank_score"], reverse=True)
    dedup = []
    for h in forest_hits:
        if any(haversine_km(h["lat"], h["lon"], k["lat"], k["lon"]) * 1000.0 < 260.0 for k in dedup):
            continue
        dedup.append(h)

    print(f"Unique closed-forest grammar complexes (>=260m spacing): {len(dedup)}", flush=True)

    # Guarantee ALL detections on l0s473 & l0s472 + top detections on l0s474..l0s477 + top Chiaps1 detections!
    on_473_472 = [h for h in dedup if h["strip_id"] in ("ChiapsCampeche_NFI_l0s473", "ChiapsCampeche_NFI_l0s472")][:8]
    on_474_477 = [h for h in dedup if re.match(r"ChiapsCampeche_NFI_l0s47[4-8]$", h["strip_id"])][:16]
    on_chiaps1 = [h for h in dedup if h["strip_id"].startswith("Chiaps1_")][:10]

    # Combine preserving uniqueness
    seen = set()
    shortlist_raw = []
    for h in on_473_472 + on_474_477 + on_chiaps1:
        key = (round(h["lat"], 5), round(h["lon"], 5))
        if key not in seen:
            seen.add(key)
            shortlist_raw.append(
                {
                    "strip_id": h["strip_id"],
                    "dem_path": h["dem_path"],
                    "chm_path": h["chm_path"],
                    "lat": round(h["lat"], 6),
                    "lon": round(h["lon"], 6),
                    "utm_e": h["utm_e"],
                    "utm_n": h["utm_n"],
                    "score_1m": h["arch_rank_score"],
                    "grammar_raw": h["grammar_raw"],
                    "n_mounds": h["n_mounds"],
                    "n_rect": h["n_rect_mounds"],
                    "n_quads": h["n_quadrants"],
                    "mean_rect": h["mean_rect"],
                    "mean_ra": h["mean_ra"],
                    "ortho_coh": h["ortho_coherence"],
                    "h1_m": h["max_height_m"],
                    "h2_m": h["second_height_m"],
                    "plaza_m2": h["plaza_area_m2"],
                    "basal_acropolis_m": h["basal_acropolis_m"],
                    "has_sacbe_hint": h["has_sacbe_hint"],
                    "has_aguada": h["has_aguada"],
                    "chm_med_m": h["chm_med_m"],
                    "chm_p90_m": h["chm_p90_m"],
                    "cleared_frac": h["chm_cleared_frac"],
                    "dist_from_card08_km": h["dist_from_card08_km"],
                }
            )

    needed_strips = sorted({c["strip_id"] for c in shortlist_raw})
    print(f"Shortlisted {len(shortlist_raw)} grammar complexes across {len(needed_strips)} strips: {needed_strips}", flush=True)

    laz_pts_by_strip = {}
    for sid in needed_strips:
        lp = ensure_laz_only(sid)
        las = laspy.read(lp)
        laz_pts_by_strip[sid] = (
            np.array(las.x, dtype=np.float64),
            np.array(las.y, dtype=np.float64),
            np.array(las.z, dtype=np.float32),
            np.array(las.classification, dtype=np.uint8),
        )

    verified = []
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(verify_candidate_0p33m_and_laz, c, laz_pts_by_strip, catalog) for c in shortlist_raw]
        for f in futs:
            verified.append(f.result())

    verified.sort(key=lambda x: (x["passes_laz_ground_gate"], x["score_1m"]), reverse=True)
    (OUT_DIR / "grammar_lacandon_verified.json").write_text(json.dumps(verified, indent=2))

    passing = [v for v in verified if v["passes_laz_ground_gate"]]
    failing = [v for v in verified if not v["passes_laz_ground_gate"]]
    print(f"\n=== GRAMMAR + .LAZ AUDIT: {len(passing)} PASSED | {len(failing)} VOID-REJECTED ===", flush=True)
    for idx, v in enumerate(verified, 1):
        st = "PASS" if v["passes_laz_ground_gate"] else "VOID-FAIL"
        nk = v["nearest_any_known"]
        ea = v["nearest_by_source"].get("EAAMS", {"name": "-", "dist_km": 99})
        print(
            f"G{idx:02d} [{st:9s}] {v['strip_id']:26s} ({v['lat']:.5f}, {v['lon']:.5f}) | "
            f"d#08={v['dist_from_card08_km']:5.2f}km | rk={v['score_1m']:5.1f} | "
            f"M={v['n_mounds']}(R:{v['n_rect']},Q:{v['n_quads']}) H={v['h1_m']:.1f}m Plz={v['plaza_m2']:.0f}m2 | "
            f"LAZ={v['laz_pts_270m']:6d} ({v['laz_pts_per_m2']:.2f}/m2, mnd={v['mound_laz_pts']}p, void={v['pct_mound_void_gt6m']:.0f}%) | "
            f"Sat={v['sat_1800m_forest_pct']:.0f}% CHM={v['chm_med_m']:.0f}m | "
            f"EAAMS={ea['name']}({ea['dist_km']}km) | Lit/OSM={nk['name'][:22]}({nk['dist_km']}km)"
        )

    # Render LARGE 2-column x 5-row (10 per page) crisp inspection sheets for ALL passing candidates!
    ext = [-135, 135, -135, 135]
    chunk_sz = 10
    for page_idx in range(math.ceil(len(passing) / chunk_sz)):
        sub = passing[page_idx * chunk_sz : (page_idx + 1) * chunk_sz]
        nrows = math.ceil(len(sub) / 2)
        fig, axes = plt.subplots(nrows, 8, figsize=(24, 3.6 * nrows), dpi=145, facecolor="#0b0f17")
        if nrows == 1:
            axes = axes[None, :]
        for i, v in enumerate(sub):
            gid = page_idx * chunk_sz + i + 1
            r_i, c_i = divmod(i, 2)
            ax_g1 = axes[r_i, c_i * 4 + 0]
            ax_lrm = axes[r_i, c_i * 4 + 1]
            ax_laz = axes[r_i, c_i * 4 + 2]
            ax_sat = axes[r_i, c_i * 4 + 3]

            npz = np.load(v["cache_npz"])
            ax_g1.imshow(npz["g1_33"], extent=ext, origin="upper")
            ax_lrm.imshow(npz["lrm33"], extent=ext, origin="upper", cmap="inferno", vmin=-1.2, vmax=4.8)
            ax_laz.imshow(npz["g1_33"], extent=ext, origin="upper", alpha=0.75)
            ax_laz.scatter(npz["px"], npz["py"], s=1.2, c="#00ffcc", alpha=0.75, linewidths=0)
            ax_sat.imshow(npz["sat_1800"])
            ax_sat.axhline(160, color="#ffd166", lw=0.8, ls="--")
            ax_sat.axvline(160, color="#ffd166", lw=0.8, ls="--")

            ea = v["nearest_by_source"].get("EAAMS", {"name": "-", "dist_km": 99})
            nk = v["nearest_any_known"]
            ax_g1.set_title(f"G{gid:02d} {v['strip_id'].replace('ChiapsCampeche_NFI_','CC_').replace('Chiaps1_','C1_')}\n({v['lat']:.5f}, {v['lon']:.5f})", color="#ffd166", fontsize=8.5, fontweight="bold")
            ax_lrm.set_title(f"M:{v['n_mounds']}(R:{v['n_rect']}) {v['n_quads']}/4Q\nH:{v['h1_m']:.1f}m Plz:{v['plaza_m2']:.0f}m²", color="#f8fafc", fontsize=8.5, fontweight="bold")
            ax_laz.set_title(f"LAZ:{v['laz_pts_270m']:,} ({v['laz_pts_per_m2']:.2f}/m²)\nMnd:{v['mound_laz_pts']:,}p Void:{v['pct_mound_void_gt6m']:.0f}%", color="#4ade80", fontsize=8.5, fontweight="bold")
            ax_sat.set_title(f"Sat:{v['sat_1800m_forest_pct']:.0f}% CHM:{v['chm_med_m']:.0f}m d#08:{v['dist_from_card08_km']}km\n{ea['name'][:14]}:{ea['dist_km']}km | {nk['name'][:14]}:{nk['dist_km']}km", color="#93c5fd", fontsize=7.8, fontweight="bold")

            for ax in (ax_g1, ax_lrm, ax_laz, ax_sat):
                ax.set_xticks([])
                ax.set_yticks([])

        for i in range(len(sub), nrows * 2):
            r_i, c_i = divmod(i, 2)
            for k in range(4):
                axes[r_i, c_i * 4 + k].axis("off")

        fig.tight_layout(pad=0.7)
        p_out = OUT_DIR / f"grammar_screen_page{page_idx+1}.png"
        fig.savefig(p_out, facecolor="#0b0f17")
        plt.close(fig)
        print("Saved:", p_out, flush=True)


if __name__ == "__main__":
    main()
