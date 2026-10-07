#!/usr/bin/env python3
"""
scripts/05a_cache_gene_sets.py
Downloads, curates, and caches gene set libraries locally for Horizon 4:
1. MSigDB Hallmark 2020 (Mouse) -> data/reference/msigdb_hallmark_mouse.json
2. KEGG 2019 (Mouse) -> data/reference/kegg_mouse.json
3. TRRUST Transcription Factors 2019 (Mouse) -> data/reference/trrust_mouse.json
4. Curated Microglial Phenotypic Signatures -> data/reference/microglia_phenotypes.json
"""

import os
import json
import logging
import gseapy as gp

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("GeneSet-Cache")

CURATED_MICROGLIA_PHENOTYPES = {
    "Microglia_Homeostatic_Signature": [
        "Tmem119", "P2ry12", "Cx3cr1", "Hexb", "Sall1", "Csf1r", "Selplg", "Sparc",
        "Gpr34", "Olfml3", "Mertk", "Fcrls", "Tgfbr1", "Siglech", "C1qa", "C1qb", "C1qc"
    ],
    "Disease_Associated_Microglia_DAM": [
        "Trem2", "Apoe", "Tyrobp", "Axl", "Cst7", "Lpl", "Clec7a", "Itgax", "Cd9",
        "Csf1", "Ctsb", "Ctsd", "Ctsl", "Spp1", "Fth1", "Anxa3", "B2m"
    ],
    "Interferon_Responsive_Microglia_IRM": [
        "Ifit1", "Ifit2", "Ifit3", "Irf7", "Stat1", "Stat2", "Oasl2", "Isg15", "Usp18",
        "Cxcl10", "B2m", "Mx1", "Mx2", "Oas1a", "Zbp1"
    ],
    "LPS_Acute_Neuroinflammatory_Shock": [
        "Tnf", "Il1b", "Il6", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Nfkb1", "Nfkbiz",
        "Nos2", "Ptgs2", "Il1a", "Cxcl1", "Cxcl2", "Icam1", "Tnfaip3", "Csf3"
    ],
    "SCFA_Metabolite_Responsive_Regulon": [
        "Ffar2", "Ffar3", "Hcar2", "Hdac1", "Hdac2", "Hdac3", "Tsc22d3", "Ddit4",
        "Plin3", "Slfn2", "Sap30", "Card6", "Llgl2", "Clu", "Fosb", "Fos"
    ]
}

def cache_libraries(out_dir: str = "data/reference"):
    os.makedirs(out_dir, exist_ok=True)

    # 1. MSigDB Hallmark
    hallmark_path = os.path.join(out_dir, "msigdb_hallmark_mouse.json")
    if not os.path.exists(hallmark_path):
        logger.info("Fetching MSigDB_Hallmark_2020 (Mouse)...")
        hallmark = gp.get_library("MSigDB_Hallmark_2020", organism="Mouse")
        with open(hallmark_path, "w") as f:
            json.dump(hallmark, f, indent=2)
        logger.info(f"Cached {len(hallmark)} Hallmark gene sets -> {hallmark_path}")
    else:
        logger.info(f"Found cached Hallmark library at {hallmark_path}")

    # 2. KEGG Mouse
    kegg_path = os.path.join(out_dir, "kegg_mouse.json")
    if not os.path.exists(kegg_path):
        logger.info("Fetching KEGG_2019_Mouse...")
        kegg = gp.get_library("KEGG_2019_Mouse", organism="Mouse")
        with open(kegg_path, "w") as f:
            json.dump(kegg, f, indent=2)
        logger.info(f"Cached {len(kegg)} KEGG pathways -> {kegg_path}")
    else:
        logger.info(f"Found cached KEGG library at {kegg_path}")

    # 3. TRRUST Transcription Factors
    trrust_path = os.path.join(out_dir, "trrust_mouse.json")
    if not os.path.exists(trrust_path):
        logger.info("Fetching TRRUST_Transcription_Factors_2019 (Mouse)...")
        trrust = gp.get_library("TRRUST_Transcription_Factors_2019", organism="Mouse")
        with open(trrust_path, "w") as f:
            json.dump(trrust, f, indent=2)
        logger.info(f"Cached {len(trrust)} TRRUST TF regulons -> {trrust_path}")
    else:
        logger.info(f"Found cached TRRUST library at {trrust_path}")

    # 4. Curated Microglia Phenotypes
    pheno_path = os.path.join(out_dir, "microglia_phenotypes.json")
    with open(pheno_path, "w") as f:
        json.dump(CURATED_MICROGLIA_PHENOTYPES, f, indent=2)
    logger.info(f"Cached {len(CURATED_MICROGLIA_PHENOTYPES)} curated microglial phenotypic states -> {pheno_path}")

if __name__ == "__main__":
    cache_libraries()
