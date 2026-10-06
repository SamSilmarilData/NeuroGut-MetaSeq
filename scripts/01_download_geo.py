#!/usr/bin/env python3
"""
scripts/01_download_geo.py
Automated acquisition of public RNA-seq count matrices and metadata from NCBI GEO.
Reads configurations from config/datasets.yaml and downloads to data/raw/<cohort>/.
"""

import os
import sys
import yaml
import gzip
import shutil
import urllib.request
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("GEO-Downloader")

def download_file(url: str, dest_path: str):
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 0:
        logger.info(f"[SKIP] File already exists: {dest_path}")
        return

    logger.info(f"[DOWNLOAD] Fetching {url} -> {dest_path}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Pipeline)"})
    try:
        with urllib.request.urlopen(req) as response, open(dest_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        logger.info(f"[OK] Downloaded {dest_path} ({os.path.getsize(dest_path):,} bytes)")
    except Exception as e:
        logger.error(f"[ERROR] Failed to download {url}: {e}")
        raise

def main():
    config_path = "config/datasets.yaml"
    if not os.path.exists(config_path):
        logger.error(f"Configuration file not found: {config_path}")
        sys.exit(1)

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    cohorts = config.get("cohorts", {})
    logger.info(f"Loaded {len(cohorts)} cohort configurations from {config_path}")

    for cohort_id, cohort_info in cohorts.items():
        logger.info(f"--- Processing {cohort_id}: {cohort_info.get('title')} ---")
        files = cohort_info.get("files", {})
        for file_key, file_meta in files.items():
            url = file_meta.get("url")
            local_path = file_meta.get("local_path")
            if url and local_path:
                try:
                    download_file(url, local_path)
                except Exception as e:
                    logger.warning(f"Could not download {file_key} for {cohort_id}: {e}")

    logger.info("[SUCCESS] All configured datasets processed.")

if __name__ == "__main__":
    main()
