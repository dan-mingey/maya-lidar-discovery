#!/usr/bin/env python3
"""Download and extract all 13 NASA G-LiHT Canopy Height Model (CHM) tarballs
for the Yucatan/Campeche/Chiapas 2013 campaigns so Part 2 & Part 3 can verify
exact 2013 forest canopy height over every candidate structure and plaza.
"""
from __future__ import annotations

import os
import tarfile
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

OUT_DIR = Path("/agent/lidar/work/carbon_raw_hunt/chm")
TAR_DIR = OUT_DIR / "tarballs"
EXTRACTED = OUT_DIR / "extracted"
TAR_DIR.mkdir(parents=True, exist_ok=True)
EXTRACTED.mkdir(parents=True, exist_ok=True)

CAMPAIGNS = [
    "AMIGACarb_Yuc_Centro_GLAS_Apr2013",
    "AMIGACarb_Yuc_Centro_NFI_Apr2013",
    "AMIGACarb_Yuc_Norte_GLAS_Apr2013",
    "AMIGACarb_Yuc_Norte_NFI_Apr2013",
    "AMIGACarb_Yuc_North_2_GLAS_Apr2013",
    "AMIGACarb_Yuc_North_2_NFI_Apr2013",
    "AMIGACarb_Yuc_South_GLAS_Apr2013",
    "AMIGACarb_Yuc_South_NFI_Apr2013",
    "AMIGACarb_Out_of_the_Yuc_GLAS_May2013",
    "AMIGACarb_Out_of_the_Yuc_NFI_May2013",
    "AMIGACarb_Chiap_Campeche_NFI_Apr2013",
    "AMIGACarb_Chiaps_1_GLAS_Apr2013",
    "AMIGACarb_Chiaps_1_NFI_Apr2013",
]


def fetch_one(campaign: str) -> str:
    url = f"https://glihtdata.gsfc.nasa.gov/files/G-LiHT/{campaign}/lidar/geotiff/{campaign}_CHM.tar.gz"
    tar_path = TAR_DIR / f"{campaign}_CHM.tar.gz"
    dest = EXTRACTED / campaign
    if dest.exists() and any(dest.rglob("*.tif")):
        n = len(list(dest.rglob("*.tif")))
        return f"SKIP {campaign} ({n} tifs)"
    dest.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, tar_path)
    with tarfile.open(tar_path, "r:gz") as tf:
        tf.extractall(dest)
    # Remove tarball after extraction to keep disk tidy
    tar_path.unlink(missing_ok=True)
    n = len(list(dest.rglob("*.tif")))
    return f"OK {campaign} ({n} CHM tifs)"


def main() -> None:
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(fetch_one, c): c for c in CAMPAIGNS}
        for fut in as_completed(futs):
            print(fut.result(), flush=True)
    total = len(list(EXTRACTED.rglob("*.tif")))
    print(f"TOTAL extracted CHM GeoTIFFs: {total}")


if __name__ == "__main__":
    main()
