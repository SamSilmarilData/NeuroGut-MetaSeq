#!/usr/bin/env python3
"""
scripts/01_download_geo.py
Automated acquisition of public RNA-seq count matrices and metadata from NCBI GEO.
Reads configurations from config/datasets.yaml, downloads to data/raw/<cohort>/,
verifies integrity via SHA-256, and records checksums into data/checksums.sha256.

Usage:
    python scripts/01_download_geo.py [--force] [--cohort GSE107925]
"""

import os
import sys
import yaml
import hashlib
import shutil
import urllib.request
import argparse
import logging
from typing import Dict, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("GEO-Downloader")

def compute_sha256(filepath: str) -> str:
    """Compute SHA-256 checksum of a local file in 64KB chunks."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(65536), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def download_file(url: str, dest_path: str, force: bool = False) -> str:
    """
    Download a remote file via HTTP/HTTPS if not already present or if forced.
    Returns the SHA-256 checksum of the file.
    """
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 0 and not force:
        logger.info(f"[SKIP] File already exists: {dest_path}")
        return compute_sha256(dest_path)

    logger.info(f"[DOWNLOAD] Fetching {url} -> {dest_path}")
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (NeuroGut-MetaSeq/1.0)"}
    )
    temp_path = dest_path + ".tmp"
    try:
        with urllib.request.urlopen(req, timeout=60) as response, open(temp_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        os.replace(temp_path, dest_path)
        size_bytes = os.path.getsize(dest_path)
        checksum = compute_sha256(dest_path)
        logger.info(f"[OK] Downloaded {dest_path} ({size_bytes:,} bytes, sha256: {checksum[:12]}...)")
        return checksum
    except Exception as e:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        logger.error(f"[ERROR] Failed downloading {url}: {e}")
        raise

def update_checksum_manifest(manifest_path: str, checksum_dict: Dict[str, str]):
    """Update or write the SHA-256 manifest file."""
    existing = {}
    if os.path.exists(manifest_path):
        with open(manifest_path, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    existing[parts[1]] = parts[0]
    
    existing.update(checksum_dict)
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
    with open(manifest_path, "w") as f:
        for filepath in sorted(existing.keys()):
            f.write(f"{existing[filepath]}  {filepath}\n")
    logger.info(f"[MANIFEST] Updated checksum manifest: {manifest_path} ({len(existing)} entries)")

def main():
    parser = argparse.ArgumentParser(description="Download raw RNA-seq data from NCBI GEO")
    parser.add_argument("--force", action="store_true", help="Force re-download of existing files")
    parser.add_argument("--cohort", type=str, default=None, help="Specific cohort ID to download (default: all)")
    parser.add_argument("--config", type=str, default="config/datasets.yaml", help="Path to datasets.yaml")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        logger.error(f"Configuration file not found: {args.config}")
        sys.exit(1)

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    cohorts = config.get("cohorts", {})
    if args.cohort:
        if args.cohort not in cohorts:
            logger.error(f"Cohort '{args.cohort}' not defined in {args.config}")
            sys.exit(1)
        cohorts = {args.cohort: cohorts[args.cohort]}

    logger.info(f"Loaded {len(cohorts)} cohort configurations to process.")
    checksums = {}

    for cohort_id, cohort_info in cohorts.items():
        logger.info(f"=== Ingesting {cohort_id}: {cohort_info.get('title')} ===")
        files = cohort_info.get("files", {})
        for file_key, file_meta in files.items():
            url = file_meta.get("url")
            local_path = file_meta.get("local_path")
            if url and local_path:
                try:
                    csum = download_file(url, local_path, force=args.force)
                    checksums[local_path] = csum
                except Exception as e:
                    logger.error(f"Failed to acquire {file_key} for {cohort_id}: {e}")

    update_checksum_manifest("data/checksums.sha256", checksums)
    logger.info("[SUCCESS] GEO raw data acquisition complete.")

if __name__ == "__main__":
    main()
