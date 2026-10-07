#!/usr/bin/env python3
"""
scripts/05c_coexpression_network.py
Weighted Co-Expression Network Analysis (WGCNA) & Hub Interactome (Horizon 4):
1. Ingests normalized log2(CPM+1) expression across all 60 biological samples.
2. Selects top ~3,500 most variable omnipresent genes (anchored by consensus & landmark genes).
3. Computes adjacency matrix (soft-threshold beta=6) and Topological Overlap Matrix (TOM).
4. Identifies co-expression modules (M1: Quiescence, M2: PolarityStress, M3: InflammatoryPriming, M4: MetabolicLipid).
5. Computes Module Eigengenes (MEs) and correlates with experimental traits.
6. Identifies intramodular hub genes and builds NetworkX interactome subgraph.
Outputs:
    results/networks/coexpression_module_assignments.csv
    results/networks/module_trait_correlations.csv
    results/networks/hub_genes_summary.csv
    results/networks/consensus_network_edges.csv
"""

import os
import glob
import logging
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.decomposition import PCA
import networkx as nx

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("WGCNA-Network")

def load_and_standardize_expression(counts_dir: str = "data/processed",
                                     meta_dir: str = "data/metadata"):
    """Loads, CPM-normalizes, and z-score standardizes expression across cohorts."""
    count_files = sorted(glob.glob(os.path.join(counts_dir, "*_counts.csv")))
    if not count_files:
        raise FileNotFoundError(f"No count files found in {counts_dir}")

    cohort_exprs = {}
    cohort_metas = {}
    gene_sets = []

    for cf in count_files:
        cname = os.path.basename(cf).replace("_counts.csv", "")
        df = pd.read_csv(cf).drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")
        sample_cols = [c for c in df.columns if c != "gene_id"]

        # log2(CPM + 1)
        lib_sizes = df[sample_cols].sum(axis=0)
        cpm = np.log2((df[sample_cols].div(lib_sizes, axis=1) * 1e6) + 1.0)

        # Z-score standardization within cohort
        z_cpm = cpm.sub(cpm.mean(axis=1), axis=0).div(cpm.std(axis=1).replace(0, 1), axis=0)

        meta_f = os.path.join(meta_dir, f"{cname}_metadata.csv")
        mdf = pd.read_csv(meta_f).set_index("sample_id")

        cohort_exprs[cname] = z_cpm
        cohort_metas[cname] = mdf
        gene_sets.append(set(z_cpm.index))

    # Common omnipresent genes across all 4 cohorts
    omnipresent = sorted(list(set.intersection(*gene_sets)))
    logger.info(f"Identified {len(omnipresent)} omnipresent expressed genes across all 4 cohorts.")

    # Combine matrices
    combined_dfs = [cohort_exprs[c].loc[omnipresent] for c in sorted(cohort_exprs.keys())]
    expr_matrix = pd.concat(combined_dfs, axis=1)

    # Combine metadata
    meta_matrix = pd.concat([cohort_metas[c] for c in sorted(cohort_metas.keys())], axis=0)
    meta_matrix = meta_matrix.loc[expr_matrix.columns]

    return expr_matrix, meta_matrix

def compute_tom(adj_matrix: np.ndarray) -> np.ndarray:
    """Computes Topological Overlap Matrix (TOM) from adjacency matrix."""
    L = np.dot(adj_matrix, adj_matrix)
    k = np.sum(adj_matrix, axis=1)
    k_min = np.minimum.outer(k, k)
    tom = (L + adj_matrix) / (k_min + 1.0 - adj_matrix)
    np.fill_diagonal(tom, 1.0)
    return np.clip(tom, 0.0, 1.0)

def run_coexpression_analysis(out_dir: str = "results/networks",
                              meta_summary_path: str = "results/meta_results/microglia_meta_analysis_summary.csv"):
    os.makedirs(out_dir, exist_ok=True)
    logger.info("Initializing Co-Expression Network Analysis (WGCNA) & Hub Deconvolution...")

    expr_matrix, meta_matrix = load_and_standardize_expression()
    logger.info(f"Total biological samples: {expr_matrix.shape[1]} across 4 cohorts.")

    # Prioritize consensus & landmark genes
    landmarks = [
        "Fosb", "Clu", "1700028E10Rik", "Llgl2", "C530043K16Rik", "Slfn2", "Neat1", "Ppif",
        "Sap30", "Card6", "Tsc22d3", "Ddit4", "Plin3", "Tnf", "Sall1", "Ffar2", "Cx3cr1", "Tmem119", "Fos"
    ]
    avail_landmarks = [g for g in landmarks if g in expr_matrix.index]

    # Select top variable genes
    variances = expr_matrix.var(axis=1)
    target_n = 3500
    top_var_genes = variances.nlargest(target_n).index.tolist()
    selected_genes = sorted(list(set(top_var_genes).union(set(avail_landmarks))))
    logger.info(f"Retained {len(selected_genes)} genes for network construction (including {len(avail_landmarks)} landmarks).")

    sub_expr = expr_matrix.loc[selected_genes]
    X = sub_expr.values.T  # Shape: (N_samples, N_genes)

    # 1. Pearson Correlation & Adjacency (beta = 6)
    logger.info("Computing Pearson correlation matrix and soft-thresholded adjacency (beta=6)...")
    corr = np.corrcoef(sub_expr.values)
    np.fill_diagonal(corr, 1.0)
    beta = 6
    adj = np.power(np.abs(corr), beta)

    # 2. Topological Overlap Matrix (TOM)
    logger.info("Computing Topological Overlap Matrix (TOM)...")
    tom = compute_tom(adj)
    dissim_tom = 1.0 - tom

    # 3. Hierarchical Clustering into 4 Primary Biological Modules
    logger.info("Performing hierarchical clustering of topological dissimilarity...")
    # Condense dissimilarity matrix for scipy linkage
    from scipy.spatial.distance import squareform
    condensed_dist = squareform(dissim_tom, checks=False)
    Z = linkage(condensed_dist, method="average")

    # Cut tree into 4 modules
    n_modules = 4
    cluster_ids = fcluster(Z, t=n_modules, criterion="maxclust")

    # Map cluster IDs to biological names based on marker membership
    module_names = {
        1: "M1_Quiescence",
        2: "M2_PolarityStress",
        3: "M3_InflammatoryPriming",
        4: "M4_MetabolicLipid"
    }
    # Dynamically assign meaningful names based on where key landmarks land
    gene_module_df = pd.DataFrame({
        "gene_symbol": selected_genes,
        "cluster_id": cluster_ids
    })

    # Intramodular connectivity k_in
    k_in_list = []
    for idx, g in enumerate(selected_genes):
        c_id = cluster_ids[idx]
        mod_mask = (cluster_ids == c_id)
        k_in = np.sum(tom[idx, mod_mask]) - 1.0  # exclude self
        k_in_list.append(round(float(k_in), 3))
    gene_module_df["k_in"] = k_in_list

    # Name modules by their dominant landmark or highest k_in hub
    named_modules = {}
    for cid in range(1, n_modules + 1):
        mod_genes = gene_module_df[gene_module_df["cluster_id"] == cid]
        top_hub = mod_genes.sort_values(by="k_in", ascending=False).iloc[0]["gene_symbol"]
        if "Slfn2" in mod_genes["gene_symbol"].values or "Sap30" in mod_genes["gene_symbol"].values:
            mname = "M_Quiescence"
        elif "Llgl2" in mod_genes["gene_symbol"].values or "Clu" in mod_genes["gene_symbol"].values:
            mname = "M_PolarityStress"
        elif "Tnf" in mod_genes["gene_symbol"].values or "Fosb" in mod_genes["gene_symbol"].values:
            mname = "M_InflammatoryPriming"
        elif "Plin3" in mod_genes["gene_symbol"].values or "Ddit4" in mod_genes["gene_symbol"].values:
            mname = "M_MetabolicLipid"
        else:
            mname = f"M_{top_hub}"
        named_modules[cid] = mname

    gene_module_df["module_name"] = gene_module_df["cluster_id"].map(named_modules)

    # 4. Compute Module Eigengenes (MEs)
    logger.info("Computing Module Eigengenes (MEs)...")
    me_dict = {}
    for mname in gene_module_df["module_name"].unique():
        m_genes = gene_module_df[gene_module_df["module_name"] == mname]["gene_symbol"].tolist()
        m_expr = sub_expr.loc[m_genes].values.T  # (N_samples, N_mgenes)
        pca = PCA(n_components=1, random_state=42)
        me = pca.fit_transform(m_expr).flatten()
        # Sign ME so it correlates positively with mean expression
        if np.corrcoef(me, m_expr.mean(axis=1))[0, 1] < 0:
            me = -me
        me_dict[mname] = me

    me_df = pd.DataFrame(me_dict, index=sub_expr.columns)

    # 5. Module-Trait Correlations across 60 samples
    logger.info("Computing Module-Trait Correlations...")
    traits_df = pd.DataFrame(index=sub_expr.columns)
    traits_df["All_Perturbed"] = (meta_matrix["condition"] == "perturbed").astype(int)
    traits_df["GF_Perturbed"] = ((meta_matrix["cohort"] == "GSE107925") & (meta_matrix["condition"] == "perturbed")).astype(int)
    traits_df["ABX_Perturbed"] = ((meta_matrix["cohort"] == "GSE108045") & (meta_matrix["condition"] == "perturbed")).astype(int)
    traits_df["Fiber_Depleted"] = ((meta_matrix["cohort"] == "GSE186210") & (meta_matrix["condition"] == "perturbed")).astype(int)
    traits_df["Sham_Perturbed"] = ((meta_matrix["cohort"] == "GSE266602") & (meta_matrix["condition"] == "perturbed")).astype(int)

    trait_corr_records = []
    for mname in me_df.columns:
        row = {"module": mname}
        for trait in traits_df.columns:
            r_val, p_val = stats.pearsonr(me_df[mname], traits_df[trait])
            row[f"r_{trait}"] = round(float(r_val), 3)
            row[f"p_{trait}"] = float(p_val)
        trait_corr_records.append(row)

    trait_corr_df = pd.DataFrame(trait_corr_records)

    # 6. Extract Top Hub Genes per Module
    hub_records = []
    for mname, group in gene_module_df.groupby("module_name"):
        top_hubs = group.sort_values(by="k_in", ascending=False).head(10)
        for rank, (_, row) in enumerate(top_hubs.iterrows(), 1):
            hub_records.append({
                "module": mname,
                "hub_rank": rank,
                "gene_symbol": row["gene_symbol"],
                "intramodular_k_in": row["k_in"]
            })
    hubs_df = pd.DataFrame(hub_records)

    # 7. Construct NetworkX Interactome Graph around Consensus Hits
    logger.info("Constructing NetworkX interactome subgraph...")
    G = nx.Graph()
    # Add nodes
    for g in avail_landmarks:
        g_mod = gene_module_df[gene_module_df["gene_symbol"] == g]["module_name"].iloc[0] if g in gene_module_df["gene_symbol"].values else "Unassigned"
        G.add_node(g, module=g_mod)

    # Add edges between landmarks based on TOM similarity (top k or threshold >= 0.005)
    landmark_indices = [selected_genes.index(g) for g in avail_landmarks]
    edge_records = []
    for i in range(len(avail_landmarks)):
        for j in range(i + 1, len(avail_landmarks)):
            g1, g2 = avail_landmarks[i], avail_landmarks[j]
            idx1, idx2 = landmark_indices[i], landmark_indices[j]
            weight = float(tom[idx1, idx2])
            if weight >= 0.005:  # Sensitive edge threshold
                G.add_edge(g1, g2, weight=round(weight, 4))
                edge_records.append({
                    "source": g1,
                    "target": g2,
                    "tom_weight": round(weight, 4)
                })

    edges_df = pd.DataFrame(edge_records)

    # Save CSV outputs
    gene_mod_path = os.path.join(out_dir, "coexpression_module_assignments.csv")
    gene_module_df.to_csv(gene_mod_path, index=False)

    trait_corr_path = os.path.join(out_dir, "module_trait_correlations.csv")
    trait_corr_df.to_csv(trait_corr_path, index=False)

    hubs_path = os.path.join(out_dir, "hub_genes_summary.csv")
    hubs_df.to_csv(hubs_path, index=False)

    edges_path = os.path.join(out_dir, "consensus_network_edges.csv")
    edges_df.to_csv(edges_path, index=False)

    logger.info("=" * 65)
    logger.info(f"[SUCCESS] WGCNA co-expression network completed:")
    logger.info(f"  - Module assignments -> {gene_mod_path} ({len(gene_module_df)} genes)")
    logger.info(f"  - Module-trait correlations -> {trait_corr_path}")
    logger.info(f"  - Top hub genes -> {hubs_path}")
    logger.info(f"  - Network edges -> {edges_path} ({len(edges_df)} edges)")
    logger.info("=" * 65)

if __name__ == "__main__":
    run_coexpression_analysis()
