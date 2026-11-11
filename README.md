# Pan-Yucatán Airborne LiDAR Discovery Pipeline: Six Unmapped Maya Architectural Complexes Beneath the Jungle Canopy

*End-to-end computational and visual discovery pass across all **456** NASA G-LiHT / OpenTopography airborne LiDAR transects spanning Yucatán, Campeche, Quintana Roo, and Chiapas (`1,128.92 km²` net bare-earth coverage), paired 1-to-1 with **456** NASA `1 m` Canopy Height Models (`CHM`), **5,223** Electronic Atlas of Ancient Maya Sites (`EAAMS`) records, **INAH `SIC`** directories, **Esri World Imagery** optical satellite tiles, **`0.33 m` (`33 cm/px`)** bare-earth DEMs/G1s, and **raw 3D laser point clouds** (`.las` full canopy returns + `.laz` classified ground returns).*

> **Archaeological & Legal Framing:** This repository presents an **independent remote-sensing analysis of open airborne LiDAR data**, not a field-excavated discovery or an official INAH registry clearance. Under Mexican federal law (*Ley Federal sobre Monumentos y Zonas Arqueológicos, Artísticos e Históricos*), official site registration and field verification belong exclusively to **INAH** (*Instituto Nacional de Antropología e Historia*). Following Society for American Archaeology (SAA) ethical principles and INAH norms, exact GPS coordinates for the six unmapped complexes (`NC-1` through `NC-6`) are withheld from public figures and tables.

---

## Contents
1. **[Global & National Context: Where `NC-1`–`NC-6` Sit in Mexico](#1-global--national-context-where-nc-1nc-6-sit-in-mexico)**
2. **[Summary Table of the 6 Unmapped Maya Complexes (`NC-1`–`NC-6`)](#2-summary-table-of-the-6-unmapped-maya-complexes-nc-1nc-6)**
3. **[Methodology: From 456 Raw LiDAR Strips to 6 Unmapped Complexes](#3-methodology-from-456-raw-lidar-strips-to-6-unmapped-complexes)**
4. **[3D Point-Cloud Simulations & Raw `.las`/`.laz` Verification (`NC-1`–`NC-6`)](#4-3d-point-cloud-simulations--raw-laslaz-verification-nc-1nc-6)**
5. **[Known-Site Benchmark & Elimination Filter](#5-known-site-benchmark--elimination-filter)**
6. **[Reproducibility, Code & Data Credits](#6-reproducibility-code--data-credits)**

---

## 1. Global & National Context: Where `NC-1`–`NC-6` Sit in Mexico

![Where the New Maya LiDAR Discoveries Sit in Mexico & on the Globe](media/maya-discoveries-mexico-and-globe-map.png)

All six unmapped architectural complexes (`NC-1` through `NC-6`) sit in the southern and central jungle interior of the **Yucatán Peninsula, Mexico** (`Campeche` and western `Quintana Roo` states), roughly `950–1,050 km` east-southeast of **Mexico City (CDMX)** and `200–300 km` south of **Mérida** and **Cancún**.

---

## 2. Summary Table of the 6 Unmapped Maya Complexes (`NC-1`–`NC-6`)

Every complex below was cross-checked against **5,223** named sites in the *Electronic Atlas of Ancient Maya Sites* (`EAAMS`, Witschey & Brown), **INAH’s official `SIC` archaeological zone registry**, Florencia Müller’s 1959/1960 state archaeological atlases, Karl Ruppert & John Denison (1943), Ivan Šprajc’s 1996–2026 SE/Central Campeche surveys (*Chactún*, *Tamchén*, *Lagunita*, *Ocomtún*, *Minanbé*), Luke Auld-Thomas et al. (2024, *Valeriana*), OpenStreetMap, Wikidata, and all prior academic studies of the 2013 NASA G-LiHT Yucatán dataset (`Golden et al. 2016`, `Ruhl et al. 2018`, `Schroder et al. 2020`, `Hutson et al. 2021`, `Britton et al. 2025/2026`). **None of these six complexes have ever been individually named, mapped, or illustrated in published literature.**

| Site ID | G-LiHT Strip ID | State & Municipality *(Approx. Region)* | Nearest Modern Cities & Federal Highways | Nearest Published Maya Ruins (`EAAMS` / `INAH`) | Architectural Layout (`0.33 m` DEM & Raw `.laz`) | Forest Canopy (`2013 CHM` & Satellite) |
| :---: | --- | --- | --- | --- | --- | --- |
| **`NC-1`** | `South_NFI_l2s450` | Calakmul Mun., Campeche (`~18.3°N, 90.0°W`) | `28 km` S of Conhuás (`Mex 186`); `42 km` WSW of **Xpujil**; `72 km` SE of **Escárcega** | `12.1 km` NNE of *Esperanza* (Rank 3); `29 km` NW of **Calakmul** | **Twin-Pyramid Rectangular Acropolis & Courtyard:** `17` mounds (`8` rectangular, `4/4` quadrants), twin **`+4.9 m` & `+4.8 m` pyramids** + **`979 m²` leveled plaza** (`+8.55 m` total profile relief) | `16.2 m` med CHM (`p90 = 18.9 m`, **`99.5%` closed**); `100%` Esri jungle |
| **`NC-2`** | `South_NFI_l2s446` | Calakmul / Escárcega, Campeche (`~18.4°N, 90.2°W`) | `19 km` S of Constitución / Silvituc (`Mex 186`); `57 km` ESE of **Escárcega** | `8.0 km` WNW of *Tuxpeña* (Rank 4); `30 km` WSW of **Balamkú**; `55 km` SW of **Ocomtún** | **`8.5 m` Monumental Flat-Topped Pyramid & Terraced Plaza:** `#1` grammar score (`68.01`), `9` mounds (`4/4` quadrants), square-based **`+8.5 m` pyramid** (`+24.67 m` hilltop prominence) + `+4.0 m` secondary mound | `14.2 m` med CHM (`p90 = 19.2 m`, **`99.4%` closed**); `100%` Esri jungle |
| **`NC-3`** | `Centro_GLAS_l1s465` | Hopelchén / Bacalar Borderlands (`~19.2°N, 89.2°W`) | `26 km` NNW of Ejido Miguel Alemán; `71 km` NNE of **Xpujil**; `64 km` SE of **Hopelchén** | `12.1 km` NNE of *Xpaicheil* (Rank 4); `31 km` NE of **Valeriana**; `>60 km` from INAH zones | **Hilltop Square Temple-Platform & Courtyard:** `10` mounds (`4/4` quadrants), rectangular **`+5.6 m` summit temple** + `+3.8 m` flanking mound + **`476 m²` leveled hilltop plaza** (`+11.64 m` profile relief) | `12.2 m` med CHM (`94.3%` closed); `100%` Esri jungle |
| **`NC-4`** | `South_GLAS_l2s437` | Escárcega / Champotón, Campeche (`~18.8°N, 90.7°W`) | `22 km` NNE of **Escárcega**; `15 km` E of `Mex 261`; `42 km` S of **Champotón** | `19.6 km` WNW of *Lake Pan Caliente* (Rank 4); `96 km` SSW of **Edzná** *(most isolated)* | **Ultra-Remote Orthogonal Platform & Courtyard Complex:** `11` mounds (`5` crisp rectangles, `4/4` quadrants), `+3.4 m` & `+2.7 m` platforms around a **`302 m²` leveled courtyard** (`+8.20 m` relief) | `15.4 m` med CHM (`p90 = 18.5 m`); `100%` Esri jungle |
| **`NC-5`** | `South_NFI_l2s445` | Escárcega / Calakmul, Campeche (`~18.5°N, 90.2°W`) | `9 km` S of Constitución / Silvituc (`Mex 186`); `51 km` E of **Escárcega** | `10.7 km` WNW of *Chen Laguna* (Rank 3); `31 km` W of **Balamkú** | **Rectangular Platform & L-Shaped Courtyard Complex:** `12` mounds (`4/4` quadrants), L-shaped range buildings + **`+4.3 m` & `+3.7 m` platforms** around a **`454 m²` leveled plaza** (`+11.31 m` relief) | `12.1 m` med CHM (`p90 = 15.7 m`, `94.8%` closed); `100%` Esri jungle |
| **`NC-6`** | `Centro_GLAS_l1s463` | Bacalar Mun., Q. Roo (`~830 m` E of Campeche; `~19.0°N, 89.1°W`) | `8 km` NW of Ejido Miguel Alemán; `58 km` NNE of **Xpujil**; `72 km` NW of **Bacalar** | `8.6 km` NW of *Sacbecan* (Rank 4); **`20.3 km` ENE of Valeriana**; `52 km` NW of **Kinichná** | **`650×300 m` Borderland Acropolis & Platform Complex:** **`40×40 m` rimmed square platform** (`+2.2 m`) + **`140×110 m` Western Acropolis** (`+3.6 m` mounds) + `2` northern patio groups | `6–19 m` forest over ruins (`30 m` N of unpaved rural road) |

---

## 3. Methodology: From 456 Raw LiDAR Strips to 6 Unmapped Complexes

![Pan-Yucatán Airborne LiDAR Discovery Pipeline & Filter Funnel](media/raw-hunt-peninsula-overview.png)

Our pipeline processes three distinct layers of open remote-sensing data across all **456** G-LiHT Yucatán transects (`1,128.92 km²` net valid pixels):

1. **Layer 1 — Multi-Scale Bare-Earth Relief Extraction (`1 m` & `0.33 m` DEMs + `G1` Composites):**
   - For every transect in Benjamin Britton’s (2026) OpenTopography bare-earth reprocessing (`DOIs: 10.5069/G9TH8JX9` & `10.5069/G9Z899N1`), we compute normalized-convolution **Local Relief Models** at fine ($\sigma = 18\text{ m}$) and meso ($\sigma = 45\text{ m}$) scales to separate human-built platforms, pyramids, and basal acropolises from natural hillside slopes.
2. **Layer 2 — Multi-Structure Maya Architectural Grammar Scoring (`5,693` Clusters):**
   - Rather than flagging isolated bumps, our detector requires **$\ge 3$ distinct steep-flanked structures within a `115 m` radius** and scores:
     - **Orthogonal Grid Coherence ($R_{90}$):** Circular resultant length of PCA major-axis azimuths folded modulo $90^\circ$, rewarding co-aligned rectilinear range structures.
     - **Four-Quadrant Courtyard Enclosure:** Identifies leveled interior plazas ($|\text{LRM}_{18}| \le 0.42\text{ m}$, slope $\le 4.2^\circ$, area $\ge 180\text{ m}^2$) bounded on `3–4` orthogonal sides by mounds.
     - **Anti-Quarry Gate:** Rejects modern roadside gravel borrow pits (*sascaberas* with sharp negative cuts $< -1.15\text{ m}$).
3. **Layer 3 — Dual-Epoch Rainforest Canopy Gate (`456` NASA `CHM`s + Esri World Imagery):**
   - Pairing all `456` DEMs 1-to-1 with NASA’s `1 m` Canopy Height Models (`CHM`) eliminates `80.4%` (`4,576 / 5,693`) of clusters located in modern pastures, ranches, towns, or road clearings, leaving **`1,117` complexes under strict intact rainforest canopy (`CHM median ≥ 11.0 m`, `≤ 10%` open)**, verified against contemporary **Esri World Imagery** optical satellite tiles.
4. **Layer 4 — Known-Site Elimination Filter (`5,223` EAAMS Sites + `INAH SIC`):**
   - Cross-referencing all `1,117` under-canopy complexes against `5,223` published EAAMS sites and INAH zones removes **`318` already-known complexes within `< 5 km`**, leaving **`799` remote complexes at `≥ 5 km`** and **`506` at `≥ 10 km`**.
5. **Layer 5 — `0.33 m` Visual Screening & Raw 3D Laser Point-Cloud Slicing (`.las` + `.laz`):**
   - Visual inspection at `0.33 m` (`33 cm/px`) removes natural tropical cone-karst hills (*mogotes*). Finally, we pull the **raw 3D NASA `.las` (all-return canopy + ground)** and **OpenTopography `.laz` (`Class 2` bare-earth + `Class 8` near-ground)** point clouds to slice `6 m`-wide vertical cross-sections directly through each complex and render 3D orbital simulations.

---

## 4. 3D Point-Cloud Simulations & Raw `.las`/`.laz` Verification (`NC-1`–`NC-6`)

Each simulation below is generated directly from the site’s **raw 3D NASA `.las` canopy returns (`12–20 m` trees)**, **Esri optical satellite imagery**, and **`0.33 m` bare-earth `.laz` ground surface**:
- **Phase 1 (`0.0–1.1 s`):** Unbroken 3D rainforest canopy concealing the ruins.
- **Phase 2 (`1.1–3.3 s`):** Airborne LiDAR laser sweep peeling away the jungle canopy in real time.
- **Phase 3 (`3.3–7.5 s`):** Full `360°` right-handed 3D orbital flyaround (`+Z UP`) of the bare-earth pyramids, plazas, and acropolises with camera-locked sun shading and a synchronized 2D HUD.

- **Master 45-Second Compilation Video (`.mp4`, All 6 Unmapped Sites):** [media/sim-all-new-maya-complexes-showcase.mp4](media/sim-all-new-maya-complexes-showcase.mp4)

![3D LiDAR Perspective Gallery — All 6 Unmapped Maya Complexes (NC-1 Through NC-6)](media/sim-3d-perspective-gallery-all-6.png)

---

### 4.1 `NC-1`: Twin-Pyramid Rectangular Acropolis & Enclosed Courtyard (`South_NFI_l2s450`)
- **Location:** Central-Southern Campeche (`12.1 km` NNE of *Esperanza*, `29 km` NW of *Calakmul*, `28 km` S of Conhuás on `Mex 186`).
- **Architectural & Point-Cloud Proof:** `17` mounds (`8` crisp rectangles, `4/4` enclosed plaza quadrants, `ortho = 0.51`), twin **`+4.9 m` & `+4.8 m` pyramids**, and a **`979 m²` leveled plaza** (`+8.55 m` total profile relief across `1,435` `.laz` ground returns beneath `14,038` `.las` canopy returns at `16.2 m` median CHM, `99.5%` closed canopy). A companion `10`-mound platform complex (`+6.1 m` relief) sits `2.1 km` NW on the same flight strip.
- **Full `1280×720` `.mp4` Video:** [media/sim-nc1-twin-pyramid-acropolis.mp4](media/sim-nc1-twin-pyramid-acropolis.mp4)

![3D Animated Simulation of NC-1](media/sim-nc1-twin-pyramid-acropolis.gif)

![New Candidate NC-1: Twin-Pyramid Rectangular Acropolis & Courtyard](media/raw-hunt-new-candidate-1-pointcloud.png)

---

### 4.2 `NC-2`: `8.5 m` Monumental Flat-Topped Pyramid & Terraced Plaza (`South_NFI_l2s446`)
- **Location:** Balamkú–Tuxpeña Forest Interior, Campeche (`8.0 km` WNW of *Tuxpeña*, `30 km` WSW of *Balamkú*, `19 km` S of Constitución / Silvituc on `Mex 186`).
- **Architectural & Point-Cloud Proof:** Highest architectural grammar score in the 456-strip corpus (`68.01`); `9` mounds (`4/4` quadrants) anchored by a massive square-based **`+8.5 m` local-relief flat-topped pyramid** (`+24.67 m` hilltop-pyramid prominence above the `180.3 m` basin datum in `4,166` `.laz` ground returns) and a `+4.0 m` secondary structure under `19,945` `.las` canopy returns (`14.2 m` median CHM, `99.4%` closed).
- **Full `1280×720` `.mp4` Video:** [media/sim-nc2-monumental-pyramid.mp4](media/sim-nc2-monumental-pyramid.mp4)

![3D Animated Simulation of NC-2](media/sim-nc2-monumental-pyramid.gif)

![New Candidate NC-2: 8.5 m Monumental Flat-Topped Pyramid & Terraced Plaza](media/raw-hunt-new-candidate-2-pointcloud.png)

---

### 4.3 `NC-3`: Hilltop Square Temple-Platform & Leveled Courtyard (`Centro_GLAS_l1s465`)
- **Location:** Campeche–Quintana Roo Borderlands (`12.1 km` NNE of *Xpaicheil*, `31 km` NE of *Valeriana*, `>60 km` from any INAH zone).
- **Architectural & Point-Cloud Proof:** Modified hilltop civic-ceremonial center (`10` mounds, `4/4` quadrants, `ortho = 0.48`) featuring a sharp rectangular **`+5.6 m` summit temple**, a `+3.8 m` flanking mound, and a **`476 m²` leveled hilltop plaza** (`+11.64 m` stepped profile relief across `12,491` `.laz` ground returns under `21,678` `.las` canopy returns, `12.2 m` median CHM).
- **Full `1280×720` `.mp4` Video:** [media/sim-nc3-hilltop-temple-plaza.mp4](media/sim-nc3-hilltop-temple-plaza.mp4)

![3D Animated Simulation of NC-3](media/sim-nc3-hilltop-temple-plaza.gif)

![New Candidate NC-3: Hilltop Square Temple-Platform & Courtyard](media/raw-hunt-new-candidate-3-pointcloud.png)

---

### 4.4 `NC-4`: Ultra-Remote Orthogonal Platform & Courtyard Group (`South_GLAS_l2s437`)
- **Location:** West-Central Campeche Interior (`19.6 km` WNW of *Lake Pan Caliente*, `22 km` NNE of **Escárcega**).
- **Architectural & Point-Cloud Proof:** Most geographically isolated complex (`11` mounds, `5` crisp rectangles, `4/4` quadrants), with `+3.4 m` and `+2.7 m` focal platforms framing a **`302 m²` leveled courtyard** (`+8.20 m` relief across `5,727` `.laz` ground returns beneath `12,556` `.las` canopy returns at `15.4 m` median CHM).
- **Full `1280×720` `.mp4` Video:** [media/sim-nc4-remote-courtyard.mp4](media/sim-nc4-remote-courtyard.mp4)

![3D Animated Simulation of NC-4](media/sim-nc4-remote-courtyard.gif)

![New Candidate NC-4: Ultra-Remote Orthogonal Platform & Courtyard Group](media/raw-hunt-new-candidate-4-pointcloud.png)

---

### 4.5 `NC-5`: Rectangular Platform & L-Shaped Courtyard Complex (`South_NFI_l2s445`)
- **Location:** Southern Campeche Forest Interior (`10.7 km` WNW of *Chen Laguna*, `9 km` S of Constitución / Silvituc on `Mex 186`).
- **Architectural & Point-Cloud Proof:** `12` mounds (`4/4` quadrants) combining crisp L-shaped courtyard range structures with **`+4.3 m` and `+3.7 m` platforms** around a **`454 m²` leveled plaza** (`+11.31 m` profile relief across `4,975` `.laz` ground returns under `12,449` `.las` canopy returns at `12.1 m` median CHM).
- **Full `1280×720` `.mp4` Video:** [media/sim-nc5-lshaped-courtyard.mp4](media/sim-nc5-lshaped-courtyard.mp4)

![3D Animated Simulation of NC-5](media/sim-nc5-lshaped-courtyard.gif)

![New Candidate NC-5: Rectangular Platform & L-Shaped Courtyard Complex](media/raw-hunt-new-candidate-5-pointcloud.png)

---

### 4.6 `NC-6`: Western Bacalar Borderland Acropolis & `40×40 m` Rimmed Platform Complex (`Centro_GLAS_l1s463`)
- **Location:** Western Municipio de Bacalar, Quintana Roo (`~830 m` E of the Campeche state line, `8.6 km` NW of *Sacbecan*, `9.5 km` SE of *Xpaicheil*, and **`20.3 km` ENE of Valeriana**, bridging the unmapped borderland corridor between eastern Campeche and western Bacalar).
- **Architectural & Point-Cloud Proof:** A continuous **`~650 × 300 m` Maya settlement zone** comprising a **`40 × 40 m` raised square platform with perimeter rims** (`+2.2 m` local relief), a **`140 × 110 m` Western Acropolis / L-shaped plaza complex** (`+3.6 m` mounds), and **two northern rectangular patio groups** (`1,804,451` total `.las` returns = `15.0 pts/m²`; `1,020,714` `.laz` ground returns = `8.5 pts/m²` under `6–19 m` forest canopy immediately north of an unpaved rural road).
- **Full `1280×720` `.mp4` Video:** [media/sim-nc6-borderland-acropolis.mp4](media/sim-nc6-borderland-acropolis.mp4)

![3D Animated Simulation of NC-6](media/sim-nc6-borderland-acropolis.gif)

![New Candidate NC-6: 40×40 m Rimmed Platform & Western Acropolis](media/raw-hunt-new-candidate-6-pointcloud.png)

![NC-6 Full 500×650 m Settlement Neighborhood Architecture Map](media/candidate-c-complex-architecture-map.png)

---

## 5. Known-Site Benchmark & Elimination Filter

To validate our detection pipeline and ensure that our final discoveries (`NC-1`–`NC-6`) were genuinely unmapped, we intersected the exact valid-pixel swath polygons of all `456` transects against all **`5,223` published EAAMS sites** and **INAH `SIC` zones**:
- **`17` known EAAMS sites** fall at **`0.0 m` (directly inside valid LiDAR pixels)**, **`33` within `≤ 250 m`**, **`102` within `≤ 1.0 km`**, and **`259` within `≤ 3.0 km`**.
- In our multi-structure under-canopy grammar hunt (`1,117` complexes), **`318` complexes (`28.5%`)** were **automatically identified and ruled out because they lay within `< 5.0 km` of an already-discovered site** (`132` within `≤ 2.5 km`).

![Ground-Truth Benchmark on Confirmed Maya Sites: Cobá South and Nadzca'an](media/raw-hunt-known-sites-benchmark.png)

### Top Known Maya Cities & Towns Rediscovered by Our Pipeline (and Filtered Out)

| Known Site Name | State | Distance to Catalog Pin | Strip ID | Pipeline Detection Signature (Before Filtering Out as Known) | Published Reference |
| --- | --- | ---: | --- | --- | --- |
| **Healy** *(Montes Azules)* | Chiapas | `1.78 km` | `Chiaps1_GLAS_l0s446` | **Grammar Score `70.40` (Rank #1 overall):** `10` mounds (`4/4` quadrants), `5.6 m` & `5.2 m` pyramids + `3.9 m` basal acropolis under `22.3 m` rainforest | Blom & Healey / Ekholm 1992; EAAMS Rank 4 |
| **Nadzca'an** *(INAH Zone)* | Campeche | **`0.39 km` / `2.35 km`** | `South_GLAS_l3s459` | **Grammar Score `66.12` (Rank #2 overall):** `12` mounds, `7.7 m` & `6.5 m` pyramids + `+12.2 m` rectangular rainforest acropolis under `13.8 m` canopy | Carrasco & Wolf 1996; INAH SIC Zone |
| **Santa Clara** *(Major Center)* | Campeche | **`0.00 km` (Direct Hit)** | `South_GLAS_l1s406` | **`74` mounds `> 2 m` (`52` mounds `> 4 m`):** Massive pyramids up to **`9.86 m` LRM** directly inside the laser swath | Müller 1960 (*Atlas de Campeche*); EAAMS Rank 3 |
| **Cobá South: Pech Mul, Dzib Mul & Kubulté** | Quintana Roo | **`0.00 km` – `0.16 km`** | `North2_GLAS_l3s454` | **`31` mounds `> 2 m`:** **`+9.7 m` pyramid**, orthogonal courtyard plazas, and diagonal raised stone **`sacbé` causeway** | Robles Castellanos 1990 (*Cobá*) |
| **El Resbalón** *(Hieroglyphic Stairway City)* | Quintana Roo | **`0.00 km` (Direct Hit)** | `Centro_NFI_l5s409` | **`146` mounds `> 2 m` (`13` mounds `> 4 m`):** Monumental core up to `5.74 m` LRM | Cortés de Brasdefer 1984; Esparza Olguín 2022 |
| **Montículo Destruido** *(Chakanbakán N)* | Quintana Roo | `2.03 km` | `Centro_GLAS_l1s458` | **Grammar Score `58.41`:** `19` mounds (`6` rectangular, `4/4` quadrants), `4.4 m` pyramid + `aguada` under `16.7 m` canopy | EAAMS Rank 4; Hutson et al. 2021 Fig. 2 |
| **Planchón de las Figuras** | Chiapas | `2.14 km` | `Chiaps1_NFI_l1s456` | **Grammar Score `56.51`:** `7` mounds (`5.1 m` pyramid + `4.1 m` basal acropolis) under `19.0 m` Lacandón canopy | EAAMS Rank 4; Golden & Scherer |
| **Dos Aguadas** *(Calakmul Biosphere)* | Campeche | `1.27 km` | `South_GLAS_l3s456` | **Grammar Score `56.25`:** `7` mounds (`ortho = 0.71`), twin `5.6 m` & `5.4 m` pyramids + `aguada` under `12.9 m` canopy | Ruppert & Denison 1943; EAAMS Rank 3 |
| **Balancanché & Kanluum** *(Chichén Itzá East)* | Yucatán | **`0.08 km` – `1.69 km`** | `North2_GLAS_l1s412` & `l1s413` | **Grammar Score `54.43`:** `22` mounds (`1,507 m²` leveled plaza) + `90` mounds `> 2 m` at Kanluum (`1.6 km` E of Chichén Itzá) | Garza & Kurjack 1980 (*Atlas de Yucatán*) |
| **Xkuncheil** | Yucatán | `1.22 km` | `North2_NFI_l3s372` | **Grammar Score `54.23`:** `11` mounds (`4` rectangular, `4/4` quadrants) around a massive leveled plaza | Garza & Kurjack 1980; EAAMS Rank 4 |
| **Benemérito Primera Sección** | Chiapas | **`0.35 km`** | `Chiaps1_GLAS_l1s457` | **Grammar Score `54.14`:** `11` mounds (`3.7 m` & `3.5 m` pyramids + `3.4 m` acropolis); ground-truthed by Golden & Scherer | Tovalín & Ortiz 2005; Golden et al. 2016 |
| **Plan de Ayutla** | Chiapas | `2.45 km` | `ChiapsCampeche_NFI_l0s475` | **Grammar Score `53.36`:** `11` mounds (`4.7 m` & `4.4 m` structures) near the famous fortified acropolis | Martos López 2009; EAAMS Rank 4 |
| **Tanchi & Kapelchultun** *(Puuc)* | Campeche | **`0.00 km` – `0.14 km`** | `Norte_GLAS_l1s426` | **`58` mounds `> 2 m` (`37` mounds `> 4 m`):** Puuc monumental platforms up to `6.37 m` LRM | Williams-Beck 1998 |
| **Tikincacab & San José Km. 29** | Yucatán | **`0.00 km` (Direct Hits)** | `North2_GLAS_l1s411` & `l2s433` | **`126` & `45` mounds `> 2 m`:** Leveled central plazas & range buildings (`3.9–4.6 m` LRM) | Ringle, Bey et al. 2004; Kepecs 1999 |
| **Calakmul Periphery (*Central Buenfil / Multun*)** | Campeche | `2.32 km` – `9.40 km` | `South_GLAS_l3s454` & `South_NFI_l2s452` | **Grammar Score `51.46`:** Calakmul urban hinterland palace & courtyard complexes (`197` mounds on `s452`) | Folan / Reese-Taylor; Golden et al. 2016 Table 2 |
| **Chakanbakán & Laguna Om Periphery** | Quintana Roo | `4.30 km` – `4.80 km` | `Centro_GLAS_l1s455` | **4-mound enclosed courtyard (`+6.4 m` relief) + causeways:** `220` mounds logged on strip | INAH SIC Zone; Golden et al. 2016 Table 2 |

---

## 6. Reproducibility, Code & Data Credits

### Pipeline Scripts (`scripts/`)
- `scripts/part1_intersect_known_sites.py` — Extracts exact valid-pixel swath polygons across all `456` G-LiHT transects and computes distances to all `5,223` EAAMS sites and INAH `SIC` zones.
- `scripts/fetch_nasa_chm.py` — Downloads and extracts all `456` NASA G-LiHT `1 m` Canopy Height Models (`CHM`).
- `scripts/part2_grammar_hunter.py` — Runs the multi-scale Local Relief Model ($\sigma = 18\text{ m}, 45\text{ m}$), multi-structure orthogonal courtyard grammar scorer, and `CHM` forest canopy gate across all `456` transects.
- `scripts/raw_pointcloud_tools.py` & `scripts/render_part3_all_showcases.py` — Streams `0.33 m` bare-earth `DEM`/`G1` rasters and raw 3D NASA `.las` + OpenTopography `.laz` point clouds to generate the 5-panel cross-section sheets (`raw-hunt-new-candidate-1..6-pointcloud.png`).
- `scripts/render_3d_video_simulations.py` — Renders the right-handed 3D point-cloud canopy-strip and `360°` orbital flyaround simulations (`sim-nc1..nc6-*.mp4`, `.gif`, and `sim-3d-perspective-gallery-all-6.png`).
- `scripts/render_mexico_and_globe_map.py` — Generates the 4-panel global, national, peninsula, and satellite heartland locator atlas (`maya-discoveries-mexico-and-globe-map.png`).

### Data Sources & Citations
- **Airborne LiDAR Acquisition:** Cook, B. D. et al. (2013). *NASA Goddard’s LiDAR, Hyperspectral and Thermal (G-LiHT) Airborne Imager*, AMIGA-Carb Yucatán Campaign ([glihtdata.gsfc.nasa.gov](https://glihtdata.gsfc.nasa.gov/files/G-LiHT/), U.S. Public Domain).
- **Bare-Earth Reprocessing & `G1` Relief Rasters:** Britton, B. (2026). *G-LiHT Yucatan East & West 2013: Reprocessed Ground Points, DEMs and Visualizations*. Distributed by OpenTopography ([10.5069/G9TH8JX9](https://doi.org/10.5069/G9TH8JX9), [10.5069/G9Z899N1](https://doi.org/10.5069/G9Z899N1), CC BY 4.0); Britton et al. (2025), *Journal of Archaeological Method and Theory*, [10.1007/s10816-025-09741-5](https://doi.org/10.1007/s10816-025-09741-5).
- **Archaeological Context & Prior G-LiHT Literature:** Auld-Thomas, L. et al. (2024, *Antiquity*); Golden, C. et al. (2016, *JAS: Reports*); Ruhl, T., Dunning, N. P. & Carr, C. (2018, *Mexicon*); Schroder, W. et al. (2020, *JAS: Reports*); Hutson, S. R. et al. (2021, *J. Anthropol. Res.*); Šprajc, I. (1996–2026); Witschey, W. R. T. & Brown, C. T. (*EAAMS*); Müller, F. (1959, 1960, *INAH*); Pietro Miotti ([pietromiotti.com/latest/maya-site](https://pietromiotti.com/latest/maya-site/)).
