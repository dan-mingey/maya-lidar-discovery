# Open LiDAR Hunt v2 — Motagua UAV + Amazon SFB/EBA

Follow-up to the Yucatán G-LiHT / NC-1–NC-6 hunt. Same LRM + grammar / ditch-detection style pipeline, applied to **open forestry/geology LiDAR that archaeologists have not wall-to-wall digitized** the way Estrada-Belli 2025 did for G-LiHT.

## Data ingested

| Corpus | Source | Products used | Area |
| --- | --- | --- | --- |
| Motagua–Polochic UAV (Jan 2024 / pub. Jul 2026) | OpenTopography Dataspace OTDS.072026.32616.* | 8× 0.5 m (Quirigua 0.25 m) DTMs | ~40 km² |
| SFB Jacundá & Jamari (pre-2017 + 2011) | Zenodo 10.5281/zenodo.18683136 | Bare-earth `dtm.tif` tiles | multi-plot |
| EBA Acre LAZ (INPE biomass) | Zenodo 10.5281/zenodo.7689909 | Ground-classified LAZ → 1 m DTM | 4 strips |
| Prior-art catalogs | Jacobs geoglyphs KML/XLS; Quirigua; Taube/Hruby Río El Tambor jade sites | — | — |

## Pipeline (same spirit as Yucatán)

1. Bare-earth DTM / LAZ→DTM  
2. NaN-safe Gaussian **Local Relief Model**  
3. Motagua: multi-mound orthogonal / plaza grammar (architectural height cap ~1–10 m)  
4. Amazon: negative-relief ditch / ring detector (depth 0.7–5 m, diam 60–350 m)  
5. NMS + distance to known sites  
6. Esri World Imagery modernity / canopy screen  
7. Manual visual crops of shortlisted hits  

## Motagua results

**Positive control:** Quirigua Museo DTM lights up (e.g. QC-01 ~`15.2667, -89.0392`) with multi-mound orthogonal architecture — detector is working.

**Hard lesson:** first-pass scores were dominated by **fault scarps / hill ridges** (20–30 m “mounds”) and **modern towns** (La Tinta, parts of Carillabas). After architectural height capping + satellite roof/veg screening:

### Priority Motagua candidates (need field / literature GPS match)

| ID | Lat, Lon | Dataset | Notes | Prior-art |
| --- | --- | --- | --- | --- |
| **MG-Gualan A** | **`15.13427, -89.31641`** | Gualan | River-terrace mound / platform under canopy; high veg; orthogonal hints | ~33 km from Quirigua; Lower Motagua valley surveyed historically (Schortman etc.) but this exact LiDAR footprint is new (fault UAV, 2026) |
| **MG-Gualan B** | **`15.12571, -89.32880`** | Gualan | Raised rectilinear platform S of modern road | Same |
| **MG-RioTambor cluster** | **`14.93932, -89.81672`** (and neighbors) | Rio Tambor | Dense hillside platform pack under canopy | ~21 km from Carrizal Grande / Taube–Hruby–Romero 2004 **jade workshop** zone (Aguilucho, Carrizal Grande, La Ceiba). Likely same cultural landscape; exact GPS match to their 7 mapped sites not published in open text — treat as **probable related / possibly unmapped terrace groups**, not “lost city” |
| Finca del Eden | `15.54993, -88.89825` | Finca del Eden | Mid-score platform cluster | Needs more screening |

**Reject / downgrade:** La Tinta top hits sit on modern grid settlement; Carillabas high scores overlap modern village roofs in Esri imagery.

Maps: [Gualan A](https://www.google.com/maps?q=15.13427,-89.31641&t=k&z=18) · [Gualan B](https://www.google.com/maps?q=15.12571,-89.32880&t=k&z=18) · [Rio Tambor](https://www.google.com/maps?q=14.93932,-89.81672&t=k&z=17)

## Amazon results

### Strongest Acre hit (EBA)

| ID | Lat, Lon | Strip | Notes | Prior-art |
| --- | --- | --- | --- | --- |
| **AM-05** | **`-9.89700, -67.30296`** | `NP_T-0414` | ~84 m-scale rectangular / ditch feature in LRM; under/near pasture–forest edge | **2.4 km** from a Jacobs-catalogued “geoglyph”; ~25 km E of Fazenda Colorada-type Acre core. **Not a named published site at this exact coord** in open catalogs checked; sits inside the well-studied Acre earthwork province. Peripato 2023 scanned EBA basin-wide — this may be one of their 24 or a near neighbor; treat as **high-priority verify**, not guaranteed first discovery |

Map: [AM-05](https://www.google.com/maps?q=-9.89700,-67.30296&t=k&z=16)

### SFB Jacundá / Jamari (Rondônia)

Strict filter yielded many candidates **~40–60 km** from nearest Jacobs geoglyph placemarks (Jamari/Jacundá forestry plots). Visual crops of top scores (AM-01, AM-08, etc.) are mostly **natural drainage confluences**, not classic circular geoglyphs. Honest verdict: **no clean, publishable first-discovery geoglyph locked yet** in SFB tiles with the current ditch detector — region remains promising (sparse published geoglyph inventory vs Acre) but needs a stronger ring/Hough detector + human visual pass on full-tile hillshades.

## What is *not* a free lunch

- Motagua UAV = **fault** survey → dense tectonic false positives  
- Amazon ditch detector without depth/circ caps → **ravines** score as “17 m ditches”  
- Acre EBA was already mined by Peripato et al. 2023; SFB forestry ALS is the fresher corpus  
- “Undiscovered” here means **absent from open named catalogs / >2 km from Jacobs points**, not proof of zero archaeological knowledge  

## Repro paths (workspace)

- Work root: `/agent/lidar/work/open_lidar_hunt_v2/`  
- Scripts: `motagua_hunter.py`, `motagua_refine_and_crop.py`, `amazon_earthwork_hunter.py`  
- Outputs: `out/motagua/`, `out/amazon/`, `previews/`  
- Media pack: `media/open-lidar-hunt-v2/` in this store  

## Recommended next levers

1. Human visual sweep of full Motagua Gualan + Rio Tambor LRM mosaics (not just auto peaks)  
2. Ring/Hough geoglyph detector on SFB Jamari/Jacundá + more EBA Acre LAZ  
3. Pull Taube 2004 GPS appendix / IGN San Diego map to hard-exclude known jade workshop mounds  
4. Optional: download additional SFB Altamira / Acre-adjacent raw MDT years if disk allows  
