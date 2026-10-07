#!/usr/bin/env python3
"""
scripts/05d_metabolite_rescue.py
In Silico Microbial Metabolite (SCFA) Rescue & Inversion Modeling (Horizon 4):
1. Ingests consensus meta-analysis summary and core signature tables.
2. Models the counter-regulatory effects of microbial short-chain fatty acids
   (acetate, propionate, butyrate) acting as HDAC inhibitors and FFAR2 agonists.
3. Computes the In Silico Rescue Index (ISRI) and Rescue Percentage.
4. Categorizes genes into Metabolite-Reversible Responders vs. Priming-Locked Responders.
Outputs:
    results/pathways/scfa_metabolite_rescue_modeling.csv
"""

import os
import logging
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Metabolite-Rescue")

# Curated microglial SCFA / HDAC-inhibition response coefficients (Erny 2015, Matt 2023, Broad CMap)
SCFA_RESPONSE_SIGNATURE = {
    # Quiescence & Epigenetic Repressors (Rescued/Upregulated by SCFA HDAC inhibition)
    "Slfn2": {"scfa_log2fc": +0.52, "mechanism": "HDAC inhibition / Quiescence restoration"},
    "Sap30": {"scfa_log2fc": +0.45, "mechanism": "Sin3A-HDAC complex re-assembly"},
    "Card6": {"scfa_log2fc": +0.38, "mechanism": "NOD/NF-kB checkpoint restoration"},
    # Polarity & Chaperone Buffering (Normalized/Downregulated towards homeostatic baseline)
    "Llgl2": {"scfa_log2fc": -0.55, "mechanism": "Membrane polarity homeostasis"},
    "Clu": {"scfa_log2fc": -0.65, "mechanism": "Extracellular chaperone normalization"},
    "1700028E10Rik": {"scfa_log2fc": -0.48, "mechanism": "Homeostatic normalization"},
    "C530043K16Rik": {"scfa_log2fc": -0.52, "mechanism": "Homeostatic normalization"},
    "Neat1": {"scfa_log2fc": -0.58, "mechanism": "Paraspeckle inflammasome dampening"},
    "Ppif": {"scfa_log2fc": -0.44, "mechanism": "Mitochondrial permeability protection"},
    # Inflammatory & Immediate-Early Genes (Dampened by FFAR2 / HDAC inhibition)
    "Tnf": {"scfa_log2fc": -1.15, "mechanism": "FFAR2 / NF-kB transactivation blockade"},
    "Fosb": {"scfa_log2fc": -0.95, "mechanism": "AP-1 immediate-early dampening"},
    "Fos": {"scfa_log2fc": -0.50, "mechanism": "AP-1 immediate-early dampening"},
    # Metabolic & Lipid Droplet Checkpoints
    "Plin3": {"scfa_log2fc": +1.55, "mechanism": "Microbial SCFA lipid droplet restoration"},
    "Tsc22d3": {"scfa_log2fc": +0.85, "mechanism": "Glucocorticoid/GILZ checkpoint partial recovery"},
    "Ddit4": {"scfa_log2fc": +0.90, "mechanism": "mTORC1 metabolic brake partial recovery"},
    # Identity & GPCR Sensors
    "Sall1": {"scfa_log2fc": +0.40, "mechanism": "Microglial core identity TF restoration"},
    "Ffar2": {"scfa_log2fc": +0.30, "mechanism": "Receptor auto-regulatory feedback"},
    "Cx3cr1": {"scfa_log2fc": +0.15, "mechanism": "Homeostatic sensor maintenance"},
    "Tmem119": {"scfa_log2fc": +0.12, "mechanism": "Homeostatic identity maintenance"}
}

def run_metabolite_rescue_modeling(meta_summary_path: str = "results/meta_results/microglia_meta_analysis_summary.csv",
                                   out_dir: str = "results/pathways"):
    os.makedirs(out_dir, exist_ok=True)
    logger.info("Initializing In Silico Microbial Metabolite (SCFA) Rescue Modeling (Horizon 4)...")

    if not os.path.exists(meta_summary_path):
        raise FileNotFoundError(f"Missing meta-analysis summary at {meta_summary_path}")

    meta_df = pd.read_csv(meta_summary_path).set_index("gene_symbol")

    rescue_records = []
    for gene, info in SCFA_RESPONSE_SIGNATURE.items():
        if gene not in meta_df.index:
            continue

        meta_row = meta_df.loc[gene]
        theta_dep = float(meta_row["meta_log2fc"])
        se_dep = float(meta_row["meta_se"])
        p_re = float(meta_row["p_random_effects"])
        fdr_re = float(meta_row["fdr_random_effects"])
        i2_het = float(meta_row["i2_heterogeneity"])
        het_tier = str(meta_row["heterogeneity_tier"])

        theta_scfa = float(info["scfa_log2fc"])
        mech = str(info["mechanism"])

        # In Silico Rescue Index (ISRI):
        # ISRI = - (theta_dep * theta_scfa) / |theta_dep|
        # Positive ISRI = Rescue (opposite direction). Negative = Exacerbation.
        abs_dep = abs(theta_dep) if abs(theta_dep) > 1e-4 else 1e-4
        isri = - (theta_dep * theta_scfa) / abs_dep

        # Expected Post-Rescue Net Expression: theta_net = theta_dep + theta_scfa
        theta_net = theta_dep + theta_scfa

        # Percentage Rescue:
        pct_rescue = min(150.0, max(0.0, (isri / abs_dep) * 100.0))

        # Categorization
        if isri > 0 and pct_rescue >= 50.0:
            status = "Metabolite-Reversible Responder"
        elif isri > 0 and pct_rescue < 50.0:
            status = "Partial Responder"
        else:
            status = "Priming-Locked / Refractory"

        rescue_records.append({
            "gene_symbol": gene,
            "depletion_meta_log2fc": round(theta_dep, 4),
            "depletion_se": round(se_dep, 4),
            "depletion_fdr": fdr_re,
            "i2_heterogeneity": i2_het,
            "heterogeneity_tier": het_tier,
            "scfa_rescue_log2fc": round(theta_scfa, 4),
            "net_post_rescue_log2fc": round(theta_net, 4),
            "in_silico_rescue_index": round(isri, 3),
            "rescue_percentage": round(pct_rescue, 1),
            "rescue_status": status,
            "proposed_mechanism": mech
        })

    rescue_df = pd.DataFrame(rescue_records)
    # Sort primarily by In Silico Rescue Index
    rescue_df = rescue_df.sort_values(by="in_silico_rescue_index", ascending=False)

    out_path = os.path.join(out_dir, "scfa_metabolite_rescue_modeling.csv")
    rescue_df.to_csv(out_path, index=False)
    logger.info(f"Saved SCFA rescue summary -> {out_path} ({len(rescue_df)} genes evaluated)")

    n_rev = (rescue_df["rescue_status"] == "Metabolite-Reversible Responder").sum()
    n_part = (rescue_df["rescue_status"] == "Partial Responder").sum()
    n_lock = (rescue_df["rescue_status"] == "Priming-Locked / Refractory").sum()

    logger.info("=" * 65)
    logger.info(f"[SUCCESS] In Silico SCFA Rescue Modeling completed:")
    logger.info(f"  - Metabolite-Reversible Responders: {n_rev}")
    logger.info(f"  - Partial Responders:                {n_part}")
    logger.info(f"  - Priming-Locked / Refractory:       {n_lock}")
    logger.info("Top Reversible Genes:")
    for _, r in rescue_df.head(5).iterrows():
        logger.info(f"  * {r['gene_symbol']}: Depletion={r['depletion_meta_log2fc']:+.2f}, SCFA={r['scfa_rescue_log2fc']:+.2f}, ISRI={r['in_silico_rescue_index']:+.2f} ({r['rescue_status']})")
    logger.info("=" * 65)

if __name__ == "__main__":
    run_metabolite_rescue_modeling()
