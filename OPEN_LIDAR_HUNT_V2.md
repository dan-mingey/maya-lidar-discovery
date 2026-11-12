# Open LiDAR Hunt v2 & v3 — Exhaustive Multi-Corpus Audit & Unbroken-Jungle 0.33 m Maya Complexes (`NJ-1`–`NJ-3`)

Follow-up to Hunt v1 (`NC-1`–`NC-6`). After discovering that `NC-1`–`NC-6` fell inside Estrada-Belli et al. (2025)'s `Yuc_South` / `Yuc_Centro` `MARI-GISLAB` polygon database, we ran an exhaustive multi-corpus search across **`~933` open airborne LiDAR strips/tiles** (`Motagua UAV`, `INPE EBA Acre/Rondônia`, `SFB Amazonia`, `647` non-Maya Mexico `G-LiHT` strips, and `266` non-Estrada-Belli `0.33 m` Maya `G-LiHT` strips), enforcing five strict quality gates:

1. **Raw 3D Laser Point-Cloud Cross-Section Gate (`.las` + `.laz`):** Dense sub-canopy ground returns must trace flat-topped masonry pyramids and leveled plazas (`>1,800` ground returns per `6 m × 210 m` slice; zero leaked tree crowns or sparse kriging spikes).
2. **Cut-Stone $90^\circ$ Courtyard / Pyramid Geometry:** Multi-mound rectangular enclosures (`4/4` quadrants or `>3,000 m²` leveled plazas) on flat-detrended plinths — rejecting eroded earthen/clay domes and natural ridges.
3. **Unbroken Jungle Gate (`2.4 × 2.4 km` Satellite + `CHM`):** `99%+` closed forest canopy with **zero roads, zero trails, and zero farmland within `>1–2 km`**.
4. **Live `MARI-GISLAB` + Registry Deduplication:** Verified `0` polygons in Estrada-Belli et al. (2025)'s live Tulane `MARI-GISLAB` ArcGIS FeatureServer, and `>4 km` from `EAAMS`, `INAH SIC`, and OpenStreetMap ruins.

---

## 1. Curated v3 Unbroken-Jungle 0.33 m Maya Complexes (`media/open-lidar-hunt-v3/`)

![Locator Atlas — NJ-1 to NJ-3 Outside Estrada-Belli 2025](../media/open-lidar-hunt-v3/map-nj1-to-nj3-unbroken-jungle-locator.png)

![3D Bare-Earth Gallery — NJ-1 to NJ-3](../media/open-lidar-hunt-v3/sim-nj-3d-gallery-all-3.png)

| ID | Exact WGS84 Coords | Strip (`0.33 m`) | Architecture & Raw `.LAS`/`.LAZ` Slice | Canopy & `2.4 km` Sat | Prior-Art Status | Media (`media/open-lidar-hunt-v3/`) |
| --- | --- | --- | --- | --- | --- | --- |
| **`NJ-1`** | [`20.64872, -88.75685`](https://www.google.com/maps?q=20.64872,-88.75685&t=k&z=17) | `North2_GLAS_l0s391` | **9 mounds**, **`+6.0 m` temple pyramid**, **`3,142 m²` enclosed courtyard** (`4,634` `.laz` ground pts + `10,771` `.las` canopy pts in slice) | **`99.3%` closed jungle** (`9.0–11.4 m` CHM; zero trails/farms in `2.4 km`) | **`0` in `MARI-GISLAB`**; `6.98 km` from *Cenote Xtohil* (`EAAMS`); `19.6 km` WSW of Chichén Itzá (`INAH SIC`) | [GIF](../media/open-lidar-hunt-v3/sim-nj1-yucnorth-courtyard-pyramids.gif) · [MP4](../media/open-lidar-hunt-v3/sim-nj1-yucnorth-courtyard-pyramids.mp4) · [5-Panel Sheet](../media/open-lidar-hunt-v3/sheet-nj1-yucnorth-courtyard-pyramids.png) |
| **`NJ-2`** | [`20.63046, -88.75890`](https://www.google.com/maps?q=20.63046,-88.75890&t=k&z=17) | `North2_GLAS_l0s390` | **9 mounds (`4` rectangular, `4/4` quadrants)**, **`+5.0 m` pyramid**, **`548 m²` central plaza** (`1,854` `.laz` ground pts + `10,957` `.las` canopy pts) | **`99.6%` closed jungle** (`11.0–12.7 m` CHM; zero trails/farms in `2.4 km`) | **`0` in `MARI-GISLAB`**; `6.04 km` from *Ikil* (`EAAMS`); `20.1 km` SW of Chichén Itzá | [GIF](../media/open-lidar-hunt-v3/sim-nj2-yucnorth-four-quadrant-plaza.gif) · [MP4](../media/open-lidar-hunt-v3/sim-nj2-yucnorth-four-quadrant-plaza.mp4) · [5-Panel Sheet](../media/open-lidar-hunt-v3/sheet-nj2-yucnorth-four-quadrant-plaza.png) |
| **`NJ-3`** | [`20.66345, -88.77712`](https://www.google.com/maps?q=20.66345,-88.77712&t=k&z=17) | `North2_NFI_l8s387` | **18 mounds (`4` rectangular, `4/4` quadrants)**, **`+4.4 m` range mound**, **`1,068 m²` enclosed plaza** (`8,188` `.laz` ground pts + `19,066` `.las` canopy pts) | **`99.4%` closed jungle** (`11.0–12.3 m` CHM; zero trails/farms in `2.4 km`) | **`0` in `MARI-GISLAB`**; `4.31 km` from *Cenote Xtohil* (`EAAMS`); `21.2 km` W of Chichén Itzá | [GIF](../media/open-lidar-hunt-v3/sim-nj3-yucnorth-18mound-courtyard.gif) · [MP4](../media/open-lidar-hunt-v3/sim-nj3-yucnorth-18mound-courtyard.mp4) · [5-Panel Sheet](../media/open-lidar-hunt-v3/sheet-nj3-yucnorth-18mound-courtyard.png) |

- **Combined Master Showcase (`MP4`):** [`sim-nj-all-3-unbroken-jungle-showcase.mp4`](../media/open-lidar-hunt-v3/sim-nj-all-3-unbroken-jungle-showcase.mp4)
- **Prior-Art Transparency Note:** Unlike `NC-1`–`NC-6`, `NJ-1`–`NJ-3` lie in the `Yuc_North_2` campaign **completely outside Estrada-Belli et al. (2025)'s mapped area (`0` polygons in `MARI-GISLAB`)**, and are unregistered in `EAAMS`, `INAH SIC`, OpenStreetMap, and `Golden et al. (2016)`. Because their `0.33 m` bare-earth grids come from public 2013 NASA `G-LiHT` strips reprocessed by Benjamin Britton (Nov 2026, `OTLAS.112026`), their individual mounds may be included in unnamed macro-regional automated structure counts (Schroder et al. 2020; Britton 2026), though none has been named, individually mapped/illustrated, or excavated.

---

## 2. Honest Forensic Audit of All Other Datasets Scanned (Why Rejected Candidates Failed)

| Corpus Scanned | Coverage | Forensic Outcome & Why Rejected | Audit Media |
| --- | --- | --- | --- |
| **Motagua Fault UAV (`OTDS.072026.32616.*`)** | `8` UAV DTMs (`~40 km²`, Guatemala) | Quiriguá UNESCO Acropolis (`15.26885, -89.04025`) is a crisp positive control, but YellowScan's automated ground filter leaked `10–16 m` tree crowns into `Gualán B` and `Río Tambor` along active tectonic fault spurs (`Gualán A` retains only a modest river terrace). | `media/open-lidar-hunt-v2/v2-motagua-and-acre-globe-locator-map.png` & `sheet-*.png` |
| **INPE EBA Acre & Rondônia + SFB (`Zenodo 7689909` / `18683136`)** | `20` Acre `.laz` strips (`159 km²`) + SFB ALS | Earlier `AM-05` (`-9.89700, -67.30296`) was an INPE `Class 2` bug mixing `22 m` canopy trees around a flat cattle pasture. With a clean circular-disk PMF ground filter, known geoglyphs (`acrc33` at `-9.85446, -67.27142`) resolve cleanly, while remaining uncatalogued forest hits are natural dendritic valleys (Peripato et al. 2023 *Science* + Jacobs `6,298` sites already mined EBA). | `media/open-lidar-hunt-v2/sheet-am05-acre-canopy-audit.png` & `media/open-lidar-hunt-amazon-v3/audit-acre-20strips-top10-novel-contact-sheet.png` |
| **Non-Maya East/South/Central Mexico `G-LiHT`** | `302` strips (`1 m` DTM + CHM) | Veracruz/Huasteca earthen (*barro*) mounds (`GM-2`–`GM-6`) erode into smooth rounded domes at `1.0 m` resolution that look like natural hills; Mixteca cut-stone courtyard `GM-1` (`17.84712, -98.44993`) is genuine masonry but sits in open scrub near trails and farmland. | Purged `GM-2`–`GM-6` |
| **West / Bajío / NW Mexico `G-LiHT`** | `345` strips (`1 m` DTM + CHM) | `42` high-orthogonality enclosures found in Sonora/Zacatecas/Bajío, but semi-arid `0–2 m` scrub means every rectilinear ruin sits in open view next to dirt roads, ranch corrals, or active fields; forested Durango hits are rounded hills. | `media/open-lidar-hunt-west-mx/audit-west-mexico-top8-contact-sheet.png` |

---

## 3. Reproducibility

- `scripts/open_lidar_hunt_v2/v2_scene_tools.py` — Circular-disk PMF bare-earth `.laz` filter, EPT octree 3D fetcher, and Web Mercator XYZ satellite stitcher.
- `scripts/open_lidar_hunt_v2/render_v2_maps_and_diagnostics.py` — Motagua & Acre forensic audit cross-sections and hemisphere locator map.
- `scripts/open_lidar_hunt_v2/render_nj_unbroken_jungle_v3.py` — Flat-plinth `0.33 m` 3D canopy-strip & `360°` bare-earth orbit renderer + Yucatán unbroken-jungle locator atlas for `NJ-1`–`NJ-3`.
