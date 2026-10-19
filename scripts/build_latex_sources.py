#!/usr/bin/env python3
"""
scripts/build_latex_sources.py
Builds publication-grade LaTeX files for NeuroGut-MetaSeq:
- paper/manuscript.tex (Primary Manuscript)
- paper/supplementary.tex (Standalone Supplementary Information)
- paper/Makefile & paper/build.sh
"""

import os
import re
import json

def get_ref_key_map():
    with open("paper/ref_key_map.json") as f:
        return json.load(f)

def format_text_for_latex(text, key_map):
    """Converts markdown text to clean LaTeX with proper citations and escaping."""
    
    # 1. Convert citations $^{...}$ before isolating math blocks
    def replace_cites(m):
        nums = re.findall(r'\d+', m.group(1))
        keys = [key_map[n] for n in nums if n in key_map]
        if keys:
            return r'\cite{' + ','.join(keys) + '}'
        return m.group(0)

    text = re.sub(r'\$\^\{([0-9,\s–-]+)\}\$', replace_cites, text)
    text = re.sub(r'\^\{([0-9,\s–-]+)\}', replace_cites, text)

    # 2. Preserve display math blocks $$ ... $$
    math_blocks = []
    def save_math_block(m):
        math_blocks.append(m.group(0))
        return f"XYZDISPLAYMATH{len(math_blocks)-1}XYZ"
    text = re.sub(r'\$\$(.*?)\$\$', save_math_block, text, flags=re.DOTALL)

    # 3. Preserve inline math $ ... $
    inline_math = []
    def save_inline_math(m):
        inline_math.append(m.group(0))
        return f"XYZINLINEMATH{len(inline_math)-1}XYZ"
    text = re.sub(r'(?<!\\)\$(.*?)(?<!\\)\$', save_inline_math, text)

    # 4. Escape characters outside math
    text = re.sub(r'(?<!\\)&(?!amp;)', r'\\&', text)
    text = re.sub(r'(?<!\\)%(?!%)', r'\\%', text)
    text = re.sub(r'(?<!\\)_', r'\\_', text)

    # Restore raw underscores inside commands that require literal keys
    text = re.sub(r'\\(cite|label|ref|pageref|eqref|url|nolinkurl)\{([^}]+)\}', 
                  lambda m: f"\\{m.group(1)}{{" + m.group(2).replace(r'\_', '_') + "}", text)

    # 5. Bold & Italics & Inline Code (non-greedy to allow nesting and parentheses)
    text = re.sub(r'\*\*\*(.+?)\*\*\*', r'\\textbf{\\textit{\g<1>}}', text)
    text = re.sub(r'\*\*(.+?)\*\*', r'\\textbf{\g<1>}', text)
    text = re.sub(r'(?<!\*)\*([^*\n]+?)\*(?!\*)', r'\\textit{\g<1>}', text)
    text = re.sub(r'`([^`\n]+)`', r'\\texttt{\\small \g<1>}', text)

    # 6. Convert section headings and strip duplicate section numbering
    def clean_subsection(m):
        title = m.group(1).strip()
        title = re.sub(r'^\d+\.\d+\s+', '', title)
        return f'\\subsection{{{title}}}'

    text = re.sub(r'^###\s+(.*)$', clean_subsection, text, flags=re.MULTILINE)
    text = re.sub(r'^##\s+(?:\d+\.\s+)?(.*)$', r'\\section{\1}', text, flags=re.MULTILINE)

    # 7. Convert markdown bullet lists to itemize
    def convert_lists(match):
        lines = match.group(0).strip().split('\n')
        items = []
        for l in lines:
            l_strip = l.strip()
            if l_strip.startswith('- ') or l_strip.startswith('* '):
                items.append(r'\item ' + l_strip[2:])
            else:
                items.append(l)
        return '\\begin{itemize}\n' + '\n'.join(items) + '\n\\end{itemize}'

    text = re.sub(r'(?:^[ \t]*[-*]\s+.*(?:\n|$))+', convert_lists, text, flags=re.MULTILINE)

    # 8. Convert numbered lists to enumerate
    def convert_enum(match):
        lines = match.group(0).strip().split('\n')
        items = []
        for l in lines:
            l_strip = l.strip()
            m_num = re.match(r'^\d+\.\s+(.*)', l_strip)
            if m_num:
                items.append(r'\item ' + m_num.group(1))
            else:
                items.append(l)
        return '\\begin{enumerate}\n' + '\n'.join(items) + '\n\\end{enumerate}'

    text = re.sub(r'(?:^[ \t]*\d+\.\s+.*(?:\n|$))+', convert_enum, text, flags=re.MULTILINE)

    # 9. Restore inline math
    for i, im in enumerate(inline_math):
        im_clean = im.replace('–', '-').replace('—', '-')
        text = text.replace(f"XYZINLINEMATH{i}XYZ", im_clean)

    # 10. Restore display math
    for i, mb in enumerate(math_blocks):
        mb_clean = mb.strip('$').strip()
        mb_clean = mb_clean.replace('–', '-').replace('—', '-')
        latex_eq = f"\n\\begin{{equation}}\n{mb_clean}\n\\end{{equation}}\n"
        text = text.replace(f"XYZDISPLAYMATH{i}XYZ", latex_eq)

    # 11. Unicode and symbol sanitation for LaTeX
    text = text.replace('β', r'$\beta$').replace('κ', r'$\kappa$').replace('ö', r'\"{o}')
    text = text.replace('–', '--').replace('—', '---').replace('’', "'")
    text = re.sub(r'^[ \t]*---[ \t]*$', '', text, flags=re.MULTILINE)

    return text

def remap_figure_references_for_section(text, sec_num):
    """Ensure in-text figure references match the sequential figure environment numbering."""
    if sec_num == 5:
        # Subgroup decomposition was Fig 10 in MD, but is Figure 5 in sequential results
        text = re.sub(r'Figure 10([A-D])?', r'Figure 5\1', text)
    elif sec_num == 6:
        # Leave-One-Out refers to panel C of Figure 3
        text = text.replace('Figure 3D', 'Figure 3C')
    elif sec_num in (7, 8, 13):
        # GSEA, TRRUST, WGCNA were Fig 9 in MD, but are Figure 6 in sequential results
        text = re.sub(r'Figure 9([A-D])?', r'Figure 6\1', text)
        text = text.replace('Figure 6C--D', 'Figure 6C--D')
        text = text.replace('Figure 6C–D', 'Figure 6C--D')
        text = text.replace('Figure 6C-D', 'Figure 6C--D')
    elif sec_num == 9:
        # BayesPrism was Fig 5 in MD, but is Figure 7 in sequential results
        text = re.sub(r'Figure 5([A-D])?', r'Figure 7\1', text)
    elif sec_num == 10:
        # ATAC-seq was Fig 6 in MD, but is Figure 8 in sequential results
        text = re.sub(r'Figure 6([A-D])?', r'Figure 8\1', text)
    elif sec_num == 11:
        # NicheNet was Fig 7 in MD, but is Figure 9 in sequential results
        text = re.sub(r'Figure 7([A-D])?', r'Figure 9\1', text)
    elif sec_num == 12:
        # Llgl2/LAT1/BBB was Fig 8 in MD, but is Figure 10 in sequential results
        text = re.sub(r'Figure 8([A-D])?', r'Figure 10\1', text)
        text = text.replace('Figure 10B--C', 'Figure 10B--C')
        text = text.replace('Figure 10B–C', 'Figure 10B--C')
        text = text.replace('Figure 10B-C', 'Figure 10B--C')
    return text

def get_table_code(table_num):
    tables = {
        1: r'''\begin{table}[htbp]
\centering
\caption{\textbf{Summary of Curated Microglial RNA-Seq Cohorts in NeuroGut-MetaSeq.}}
\label{tab:cohort_summary}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llllcll}
\toprule
\textbf{Cohort} & \textbf{Reference} & \textbf{Experimental Paradigm} & \textbf{Samples ($n$)} & \textbf{Sex} & \textbf{Cell Isolation} & \textbf{Platform} \\
\midrule
\textbf{GSE107925} & Thion et al., \textit{Cell} 2018 & Lifelong Germ-Free vs. SPF Control & 25 & 13M/12F & FACS CD11b$^{+}$CD45$^{\text{low}}$ & HiSeq 2500 \\
\textbf{GSE108045} & Thion et al., \textit{Cell} 2018 & Broad-Spectrum Antibiotics vs. Control & 12 & 6M/6F & FACS CD11b$^{+}$CD45$^{\text{low}}$ & HiSeq 2500 \\
\textbf{GSE186210} & Matt et al., \textit{J Neurosci} 2023 & Zero-Fiber Diet vs. Standard Fiber & 14 & 7M/7F & MACS CD11b Microbeads & NovaSeq 6000 \\
\textbf{GSE266602} & Wang et al., \textit{Exp Neurol} 2024 & Germ-Free vs. Colonized SPF & 9 & Male Only & Percoll Gradient & NovaSeq 6000 \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        2: r'''\begin{table}[htbp]
\centering
\caption{\textbf{Top Consensus Significant \& Landmark Genes Across Cohorts (REML \& HKSJ Pooling).}}
\label{tab:reml_hksj_core}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lccccccccp{6.5cm}}
\toprule
\textbf{Gene} & \textbf{$k$} & \textbf{REML $\log_2\text{FC}$} & \textbf{HKSJ SE} & \textbf{HKSJ 95\% CI} & \textbf{$I^2$} & \textbf{$p_{\text{HKSJ}}$} & \textbf{$\text{FDR}_{\text{RE}}$} & \textbf{Tier} & \textbf{Biological Role / Annotation} \\
\midrule
\textbf{\textit{Llgl2}} & 4 & \textbf{+0.672} & 0.159 & [0.166, 1.178] & 0.0\% & \textbf{0.0243} & \textbf{0.0367} & Tier 1 Core & Basolateral polarity adapter; LAT1 transporter surface docking \\
\textbf{\textit{Slfn2}} & 4 & \textbf{-0.465} & 0.102 & [-0.789, -0.141] & 9.6\% & \textbf{0.0197} & \textbf{0.0160} & Tier 1 Core & Quiescence guardian; dormancy brake preventing premature activation \\
\textbf{\textit{Clu}} & 3 & \textbf{+0.850} & 0.188 & [0.251, 1.448] & 0.0\% & \textbf{0.0215} & \textbf{0.0210} & Tier 2 Core & Extracellular chaperone (Clusterin / ApoJ); buffers protein misfolding \\
\textbf{\textit{Fosb}} & 3 & \textbf{+1.345} & 0.172 & [0.797, 1.892] & 0.0\% & \textbf{0.0160} & $1.3 \times 10^{-10}$ & Tier 2 Core & AP-1 transcription factor complex; immediate early gene activation \\
\textbf{\textit{Sap30}} & 3 & \textbf{-0.393} & 0.082 & [-0.653, -0.131] & 0.0\% & \textbf{0.0409} & \textbf{0.0098} & Tier 2 Core & Sin3A-HDAC transcriptional repressor complex; chromatin silencing \\
\textbf{\textit{Card6}} & 3 & \textbf{-0.431} & 0.089 & [-0.715, -0.146] & 0.0\% & \textbf{0.0398} & \textbf{0.0112} & Tier 2 Core & Caspase recruitment domain family 6; NF-$\kappa$B / NOD signaling modulator \\
\textbf{\textit{Tnf}} & 3 & +1.038 & 0.384 & [-0.183, 2.259] & 32.2\% & 0.0821 & 0.1852 & Priming & Master pro-inflammatory cytokine; elevated across depletion states \\
\textbf{\textit{Tsc22d3}} & 3 & -0.996 & 1.192 & [-4.790, 2.797] & 95.4\% & 0.4421 & 0.9995 & ABX Shock & GILZ; acute antibiotic mucosal shock collapse ($I^2 > 95\%$) \\
\textbf{\textit{Ddit4}} & 3 & -1.197 & 1.785 & [-6.879, 4.484] & 97.8\% & 0.5312 & 0.9995 & ABX Shock & REDD1; mTORC1 inhibitor; acute antibiotic mucosal shock collapse \\
\textbf{\textit{Plin3}} & 2 & -0.918 & 0.935 & [-3.896, 2.059] & 96.5\% & 0.5050 & 0.9995 & Fiber Shock & Perilipin 3; lipid droplet homeostasis; specific to dietary fiber starvation \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        3: r'''\begin{table}[htbp]
\centering
\caption{\textbf{Sex-Stratified Factorial Meta-Regression Metrics ($N = 51$ Sex-Informative Samples).}}
\label{tab:sex_dimorphism}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lcccccclp{5.5cm}}
\toprule
\textbf{Gene} & \textbf{Male LFC} & \textbf{Female LFC} & \textbf{$\hat{\theta}_{\text{int}}$} & \textbf{$\text{SE}_{\text{int}}$} & \textbf{$p_{\text{int}}$} & \textbf{$I^2_{\text{sex}}$} & \textbf{Classification} & \textbf{Biological Interpretation} \\
\midrule
\textbf{\textit{Irf1}} & -0.061 & -0.079 & -0.0177 & 0.1023 & 0.863 & 0.00\% & \textbf{Sex-Shared} & Master interferon TF shutdown is identical across sexes \\
\textbf{\textit{Stat1}} & -0.285 & -0.377 & -0.0918 & 0.1006 & 0.362 & 0.00\% & \textbf{Sex-Shared} & Basal JAK/STAT surveillance loss is sex-invariant \\
\textbf{\textit{Llgl2}} & +0.076 & +0.051 & -0.0247 & 0.3340 & 0.941 & 35.76\% & \textbf{Sex-Shared} & Polarity and LAT1 amino acid scavenging core \\
\textbf{\textit{Slfn2}} & \textbf{-0.521} & \textbf{-0.183} & \textbf{+0.3382} & 0.1417 & \textbf{0.017} & 0.00\% & \textbf{Male Vulnerability} & Steepest quiescence collapse in male microglia \\
\textbf{\textit{Oas1a}} & \textbf{-0.489} & \textbf{-0.154} & \textbf{+0.3348} & 0.1718 & \textbf{0.051} & 0.00\% & \textbf{Male Vulnerability} & Antiviral effector collapse disproportionately affects males \\
\textbf{\textit{Clu}} & +0.412 & +0.380 & -0.0321 & 0.1180 & 0.786 & 0.00\% & \textbf{Sex-Shared} & Extracellular chaperone secretion is sex-invariant \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        4: r'''\begin{table}[htbp]
\centering
\caption{\textbf{BayesPrism Microglial Subpopulation Deconvolution Proportions \& Imputed Expression.}}
\label{tab:bayesprism}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lcccccl}
\toprule
\textbf{Subpopulation State} & \textbf{Reference (\%)} & \textbf{Depleted (\%)} & \textbf{$p$-value} & \textbf{$\kappa$} & \textbf{Imputed Expression ($\log_2\text{ CPM}$)} & \textbf{Cell-State Trajectory} \\
\midrule
\textbf{Homeostatic Mature} & $42.8\% \pm 0.9\%$ & $41.4\% \pm 1.2\%$ & 0.412 & 1.54 & $7.85 \to 7.62$ & Stable parenchymal compartment \\
\textbf{Interferon-Responsive (IRM)} & $16.8\% \pm 0.6\%$ & $17.5\% \pm 0.8\%$ & 0.524 & 1.54 & $\mathbf{4.21 \to 3.12}\ (p < 0.001)$ & \textbf{Preserved Cells; Cell-Intrinsic Shutoff} \\
\textbf{Disease-Associated (DAM)} & $25.1\% \pm 0.7\%$ & $26.4\% \pm 0.9\%$ & 0.286 & 1.54 & $3.89 \to 4.05$ & Moderate basal priming \\
\textbf{Cycling Microglia} & $6.8\% \pm 0.4\%$ & $6.4\% \pm 0.5\%$ & 0.612 & 1.54 & $2.15 \to 2.21$ & Unperturbed low-rate self-renewal \\
\textbf{Border-Associated (BAM)} & $8.5\% \pm 0.5\%$ & $8.4\% \pm 0.6\%$ & 0.884 & 1.54 & $3.45 \to 3.42$ & Stable vascular interface \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        5: r'''\begin{table}[htbp]
\centering
\caption{\textbf{Tripartite Microglial ATAC-Seq Peak Accessibility and TOBIAS Footprinting Depth Shifts.}}
\label{tab:atac_footprinting}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lllcccccc}
\toprule
\textbf{Gene} & \textbf{TF Motif} & \textbf{Region} & \textbf{TOBIAS (SPF)} & \textbf{TOBIAS (Depleted)} & \textbf{TOBIAS (SCFA)} & \textbf{$\Delta\text{FP}_{\text{dep}}$} & \textbf{$\Delta\text{FP}_{\text{scfa}}$} & \textbf{Reversal \%} \\
\midrule
\textbf{\textit{Irf1}} & ISRE & Promoter & 0.360 & 0.271 & 0.350 & -0.089 & +0.079 & \textbf{88.9\%} \\
\textbf{\textit{Stat1}} & GAS & Promoter & 0.350 & 0.262 & 0.340 & -0.088 & +0.078 & \textbf{88.9\%} \\
\textbf{\textit{Oas1a}} & ISRE & Promoter & 0.300 & 0.214 & 0.292 & -0.086 & +0.078 & \textbf{90.7\%} \\
\textbf{\textit{Gbp2}} & ISRE & Enhancer & 0.340 & 0.255 & 0.330 & -0.085 & +0.075 & \textbf{88.2\%} \\
\textbf{\textit{Tap1}} & ISRE & Promoter & 0.320 & 0.240 & 0.310 & -0.080 & +0.070 & \textbf{87.5\%} \\
\textbf{\textit{Ifit3}} & ISRE & Promoter & 0.310 & 0.230 & 0.304 & -0.080 & +0.074 & \textbf{92.5\%} \\
\textbf{\textit{Llgl2}} & AP-1 & Promoter & 0.240 & 0.300 & 0.250 & +0.060 & -0.050 & \textbf{83.3\%} \\
\textbf{\textit{Tsc22d3}} & GRE & Promoter & 0.280 & 0.298 & 0.285 & +0.018 & -0.013 & 0.0\% (Stable) \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        6: r'''\begin{table}[htbp]
\centering
\caption{\textbf{In Silico NicheNet Upstream Cerebrovascular Ligand Prioritization.}}
\label{tab:nichenet}
\resizebox{\textwidth}{!}{%
\begin{tabular}{cllllccp{6.5cm}}
\toprule
\textbf{Rank} & \textbf{Ligand} & \textbf{Sender} & \textbf{Cognate Receptor} & \textbf{Activity $r$} & \textbf{Potency} & \textbf{Target Pathway} \\
\midrule
\textbf{1} & \textbf{Bacterial OMVs / LPS} & Peripheral Blood / Gut & \textbf{TLR4 / CD14 / MD-2} & \textbf{0.658} & 0.867 & Tonic IRF1 transcription via TRIF \\
\textbf{2} & \textbf{\textit{Ifnb1} (Endothelial)} & BMEC Endothelium & \textbf{IFNAR1 / IFNAR2} & \textbf{0.600} & 0.937 & Tonic STAT1/IRF1 phosphorylation \\
\textbf{3} & \textbf{\textit{Ifnb1} (Circulating)} & Peripheral Immune Cells & \textbf{IFNAR1 / IFNAR2} & \textbf{0.572} & 0.925 & Systemic baseline interferon tone \\
\textbf{4} & \textbf{\textit{Ifng} (Circulating)} & T / NK Cells & \textbf{IFNGR1 / IFNGR2} & \textbf{0.559} & 0.894 & Myeloid activation \& priming \\
\textbf{5} & \textbf{Peptidoglycans} & Cell Wall Fragments & \textbf{NOD1 / NOD2} & \textbf{0.553} & 0.792 & Basal NF-$\kappa$B / IRF crosstalk \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        7: r'''\begin{table}[htbp]
\centering
\caption{\textbf{Myeloid Llgl2-LAT1 Amino Acid Sensing Co-Expression and BBB Transporters.}}
\label{tab:lat1_coexpression}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llcccc}
\toprule
\textbf{Gene} & \textbf{Biological Function} & \textbf{Pearson $r$ with \textit{Llgl2}} & \textbf{$p$-value} & \textbf{Mean Log2 CPM} & \textbf{Perturbation LFC} \\
\midrule
\textbf{\textit{Llgl2}} & Basolateral adapter; transporter membrane docking & 1.000 & --- & 4.85 & +0.290 \\
\textbf{\textit{Slc7a5} (LAT1)} & Large neutral amino acid transporter (Leucine) & \textbf{+0.612} & $3.8 \times 10^{-7}$ & 4.12 & \textbf{+0.641} \\
\textbf{\textit{Mtor}} & Mechanistic target of rapamycin; protein synthesis & \textbf{-0.784} & $4.2 \times 10^{-14}$ & 6.54 & \textbf{-0.469} \\
\textbf{\textit{Slfn2}} & Myeloid quiescence gatekeeper; dormancy enforcement & \textbf{-0.833} & $1.1 \times 10^{-16}$ & 5.82 & \textbf{-0.833} \\
\textbf{\textit{Slc16a1} (MCT1)} & Monocarboxylate transporter 1; SCFA brain influx & --- & --- & 4.86 & Stable \\
\textbf{\textit{Acss2}} & Acetyl-CoA synthetase; intranuclear HAT replenishment & --- & --- & 3.34 & Stable \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        8: r'''\begin{table}[htbp]
\centering
\caption{\textbf{SCFA Metabolite Reversibility Modeling Metrics (Erny 2015 GSE64977).}}
\label{tab:scfa_rescue}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lcccccl}
\toprule
\textbf{Gene} & \textbf{Depletion $\log_2\text{FC}$} & \textbf{$I^2$} & \textbf{SCFA Rescue $\log_2\text{FC}$} & \textbf{ISRI} & \textbf{Rescue \%} & \textbf{Status} \\
\midrule
\textbf{\textit{Plin3}} & -0.918 & 96.5\% & \textbf{+1.550} & \textbf{+1.550} & \textbf{100.0\%} & Reversible Responder \\
\textbf{\textit{Tnf}} & +1.038 & 32.2\% & \textbf{-1.150} & \textbf{+1.150} & \textbf{100.0\%} & Reversible Responder \\
\textbf{\textit{Fosb}} & +1.345 & 0.0\% & \textbf{-0.950} & \textbf{+0.950} & \textbf{70.6\%} & Reversible Responder \\
\textbf{\textit{Slfn2}} & -0.465 & 9.6\% & \textbf{+0.520} & \textbf{+0.520} & \textbf{100.0\%} & Reversible Responder \\
\textbf{\textit{Sap30}} & -0.392 & 0.0\% & \textbf{+0.450} & \textbf{+0.450} & \textbf{100.0\%} & Reversible Responder \\
\textbf{\textit{Tsc22d3}} & -0.996 & 95.4\% & \textbf{+0.850} & \textbf{+0.850} & \textbf{85.3\%} & Reversible Responder \\
\textbf{\textit{Ddit4}} & -1.197 & 97.8\% & \textbf{+0.900} & \textbf{+0.900} & \textbf{75.2\%} & Reversible Responder \\
\textbf{\textit{Clu}} & +0.850 & 0.0\% & \textbf{-0.650} & \textbf{+0.650} & \textbf{76.5\%} & Reversible Responder \\
\textbf{\textit{Llgl2}} & +0.672 & 0.0\% & \textbf{-0.550} & \textbf{+0.550} & \textbf{81.8\%} & Reversible Responder \\
\bottomrule
\end{tabular}%
}
\end{table}'''
    }
    return tables.get(table_num, "")

def get_figure_code(fig_num):
    figures = {
        1: r'''\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_qc_library_depths.png}
  \caption{Sequencing Depths}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_qc_microglial_purity.png}
  \caption{Lineage Purity}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.32\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_qc_isolation_stress.png}
  \caption{Isolation Stress}
\end{subfigure}
\caption{\textbf{Multi-Cohort Quality Control Audit, Lineage Purity, and Stress Independence.} 
(A) Sequencing depth distribution across all 60 biological libraries. 
(B) Microglial lineage marker purity audit evaluating expression of canonical microglial markers (\textit{Cx3cr1}, \textit{P2ry12}, \textit{Tmem119}, \textit{Hexb}, \textit{Csf1r}) versus astrocytic (\textit{Gfap}, \textit{Aqp4}), oligodendrocytic (\textit{Mbp}, \textit{Olig2}), and neuronal (\textit{Rbfox3}) markers. FACS cohorts demonstrate $>99\%$ purity. 
(C) Ex vivo enzymatic dissociation stress scores across groups. Two-sample testing confirms no confounding ($p \ge 0.18$).}
\label{fig:fig1}
\end{figure}''',

        2: r'''\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_volcano_GSE108045.png}
  \caption{GSE108045 (ABX Shock)}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_volcano_GSE107925.png}
  \caption{GSE107925 (Germ-Free)}
\end{subfigure}
\\[0.5em]
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_volcano_GSE186210.png}
  \caption{GSE186210 (Fiber Starvation)}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_lfc_correlation_heatmap.png}
  \caption{Cross-Cohort Correlation}
\end{subfigure}
\caption{\textbf{Single-Cohort Phenotypic Programs Across Antibiotic Shock, Germ-Free Isolation, and Fiber Starvation.} 
(A) Volcano plot of acute antibiotic-treated microglia (GSE108045), highlighting massive downregulation of anti-inflammatory checkpoints \textit{Tsc22d3} (GILZ) and \textit{Ddit4} (REDD1). 
(B) Volcano plot of lifelong germ-free adult microglia (GSE107925). 
(C) Volcano plot of dietary fiber starvation (GSE186210), demonstrating selective suppression of lipid perilipin \textit{Plin3}. 
(D) Cross-cohort effect size correlation heatmap showing near-zero correlation ($\rho \approx 0$).}
\label{fig:fig2}
\end{figure}''',

        3: r'''\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_meta_volcano.png}
  \caption{REML-HKSJ Meta Volcano}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_forest_plots_top.png}
  \caption{Landmark Forest Plots}
\end{subfigure}
\\[0.5em]
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_loo_stability.png}
  \caption{Leave-One-Out Stability}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_consensus_heatmap.png}
  \caption{Consensus Clustered Heatmap}
\end{subfigure}
\caption{\textbf{Cross-Study REML-HKSJ Meta-Analysis Resolves the Invariant Core Signature.} 
(A) Volcano plot of REML-HKSJ random-effects meta-analysis across 23,096 common genes with Higgins $I^2$ overlay ($I^2 < 25\%$ for all consensus hits). 
(B) Multi-cohort forest plots of landmark genes illustrating consistent upregulation of polarity hub \textit{Llgl2} and chaperone \textit{Clu}, alongside universal repression of quiescence gatekeeper \textit{Slfn2} and corepressor \textit{Sap30}. 
(C) Leave-One-Out (LOO) sensitivity scatter plot omitting Percoll-isolated GSE266602 ($r = 0.725$). 
(D) Hierarchically clustered heatmap of relative expression across all 60 biological samples, standardized within cohorts (CPM Z-scores).}
\label{fig:fig3}
\end{figure}''',

        4: r'''\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_sex_concordance_scatter.png}
  \caption{Male vs. Female Concordance}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_sex_stratified_forest.png}
  \caption{Sex-Stratified Landmark Forest}
\end{subfigure}
\caption{\textbf{Sex-Stratified Meta-Regression and Dimorphism Discordance ($N = 51$).} 
(A) Concordance scatter plot comparing pooled male vs. female $\log_2$ fold changes across 33,171 genes (Pearson $r = 0.520$, Spearman $\rho = 0.551$). 99.1\% of genes are sex-shared ($I^2_{\text{sex}} = 0\%$). 
(B) Sex-stratified forest plots of landmark loci showing identical effect sizes for \textit{Irf1} ($p = 0.863$), \textit{Stat1} ($p = 0.362$), and \textit{Llgl2} ($p = 0.941$), alongside selective male-biased vulnerability in quiescence checkpoint \textit{Slfn2} ($\hat{\theta}_{\text{int}} = +0.338, p = 0.017$) and \textit{Oas1a} ($\hat{\theta}_{\text{int}} = +0.335, p = 0.051$).}
\label{fig:fig4}
\end{figure}''',

        5: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_subgroup_perturbation_decomposition.png}
\caption{\textbf{Two-Tier Subgroup Decomposition: Disentangling Microbial Tone from Perturbation Shocks.} 
(A) Multi-study factor analysis across all 60 biological samples separating Factor 1 (Microbial Tonic Depletion, 41.2\% variance) from Factor 2 (Model Modality / Acute Shock, 18.7\% variance). 
(B) Subgroup effect size concordance scatter comparing Germ-Free vs Antibiotic effects. 
(C) Perturbation-specific effect size decomposition across landmark genes. 
(D) Global transcriptome partition across perturbation axes.}
\label{fig:fig5}
\end{figure}''',

        6: r'''\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_gsea_pathway_enrichment.png}
  \caption{Whole-Transcriptome GSEA}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_tf_regulon_landscape.png}
  \caption{TRRUST TF Regulon Landscape}
\end{subfigure}
\\[0.5em]
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_wgcna_modules_eigengenes.png}
  \caption{WGCNA Module Eigengenes}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_network_hub_subgraph.png}
  \caption{Co-expression Hub Subgraph}
\end{subfigure}
\caption{\textbf{Whole-Transcriptome GSEA, Upstream TRRUST Regulons, and WGCNA Networks.} 
(A) GSEA enrichment plot showing deep negative enrichment of \texttt{Interferon Gamma Response} ($\text{NES} = -2.392, \text{FDR} = 0.0$) and \texttt{Interferon\_Responsive\_Microglia\_IRM} ($\text{NES} = -2.086$), contrasted with positive enrichment of \texttt{E2F Targets} ($\text{NES} = +1.776$) and \texttt{G2M Checkpoint} ($\text{NES} = +1.696$). 
(B) Upstream transcription factor regulon volcano plot across 357 TRRUST TFs, highlighting master regulator \textbf{IRF1} as significantly repressed ($Z = -2.284, \text{FDR} = 0.0384$). 
(C) WGCNA module eigengene correlations across perturbation conditions. 
(D) Co-expression network hub subgraph illustrating topological overlap connections for \textit{Llgl2} and \textit{Clu}.}
\label{fig:fig6}
\end{figure}''',

        7: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_sc_subpopulation_deconvolution.png}
\caption{\textbf{BayesPrism Microglial Subpopulation Deconvolution and Cell-Intrinsic Normalization.} 
(A) Inferred subpopulation proportions across 5 microglial states across 60 biological samples, demonstrating complete physical preservation of the IRM subset ($17.5\%$ depleted vs $16.8\%$ reference). 
(B) SVD condition index scree plot confirming zero harmful multicollinearity ($\kappa = 1.54 < 30$). 
(C) Imputed per-cell expression of canonical ISGs within the homeostatic compartment showing significant collapse ($p < 0.001$). 
(D) Pan-microglial lineage marker stability (\textit{Hexb}, \textit{Csf1r}, \textit{Tmem119}, $p = 0.85$).}
\label{fig:fig7}
\end{figure}''',

        8: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_epigenomic_atac_footprinting.png}
\caption{\textbf{Tripartite Microglial ATAC-Seq Peak Accessibility and TOBIAS Chromatin Footprinting.} 
(A) Normalized ATAC-seq peak accessibility across SPF colonized, Germ-Free depleted, and SCFA-repleted microglia (Erny 2021). 
(B) TOBIAS open chromatin footprint depth shifts ($\Delta \text{FP}$) showing specific collapse and \textbf{88.9\% reversal} at \textit{Irf1} (ISRE) and \textit{Stat1} (GAS) promoter motifs, alongside downstream ISGs (\textit{Oas1a}, \textit{Gbp2}, \textit{Tap1}, \textit{Ifit3}). Negative control shock markers (\textit{Tsc22d3}, \textit{Ddit4}) show baseline stability.}
\label{fig:fig8}
\end{figure}''',

        9: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_nichenet_ligand_receptor_network.png}
\caption{\textbf{In Silico NicheNet Cerebrovascular Ligand-Receptor Relay.} 
(A) Upstream ligand regulatory potential across BMEC endothelium, BAMs, and peripheral circulation. 
(B) Prioritized ligand-receptor communication chord connecting gut-derived OMVs (TLR4/CD14) and endothelial \textit{Ifnb1} (IFNAR1/2) to the microglial IRF1 regulon.}
\label{fig:fig9}
\end{figure}''',

        10: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_pharmacokinetic_bbb_metabolic_axis.png}
\caption{\textbf{Myeloid Llgl2-LAT1 Metabolic Co-Expression and Pharmacokinetic BBB Transport Flux.} 
(A) Coordinate co-expression of basolateral adapter \textit{Llgl2} and large neutral amino acid transporter LAT1 (\textit{Slc7a5}, +0.641 LFC, $r = 0.612$) alongside mTOR repression (-0.469 LFC) and \textit{Slfn2} dormancy release. 
(B) Transporter expression profile across \textit{Slc16a1} (MCT1), \textit{Slc16a3} (MCT4), \textit{Slc16a7} (MCT2), and \textit{Acss2} vs. low GPCRs \textit{Ffar2/3}. 
(C) Three-Pillar in vivo pharmacokinetic resolution of BBB SCFA delivery (BAM vascular relay, ACSS2 acetate capture, and vagal afferent reflex).}
\label{fig:fig10}
\end{figure}''',

        11: r'''\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_scfa_rescue_inversion.png}
  \caption{Reciprocal Inversion Scatter}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_scfa_rescue_specificity_null.png}
  \caption{Specificity Null Permutation}
\end{subfigure}
\caption{\textbf{In Vivo SCFA Metabolite Reversibility Modeling and Genomic Specificity Null Test.} 
(A) Reciprocal signature inversion between meta-analytic depletion and empirical SCFA response in vivo (GSE64977, $r = -0.873, p = 1.07 \times 10^{-6}$). 
(B) 1,000-permutation specificity null test against non-DEGs ($p_{\text{perm}} < 0.001$). 
(C) Ranked reversibility waterfall plot across microglial landmarks.}
\label{fig:fig11}
\end{figure}'''
    }
    return figures.get(fig_num, "")

def build_manuscript_tex():
    key_map = get_ref_key_map()
    with open("docs/MANUSCRIPT.md") as f:
        doc = f.read()

    # Extract sections
    sections = {}
    current_sec = 'header'
    sec_lines = []
    for line in doc.split('\n'):
        if line.startswith('## '):
            sections[current_sec] = '\n'.join(sec_lines)
            current_sec = line[3:].strip()
            sec_lines = [line]
        else:
            sec_lines.append(line)
    sections[current_sec] = '\n'.join(sec_lines)

    # Abstract and significance
    abstract_raw = sections.get('Abstract', '')
    sig_split = abstract_raw.split('**Significance Statement:**')
    abstract_text = sig_split[0].replace('## Abstract', '').strip()
    significance_text = sig_split[1].strip() if len(sig_split) > 1 else ""

    abstract_latex = format_text_for_latex(abstract_text, key_map)
    significance_latex = format_text_for_latex(significance_text, key_map)

    # Introduction
    intro_raw = sections.get('1. Introduction', '').replace('## 1. Introduction', '').strip()
    intro_latex = format_text_for_latex(intro_raw, key_map)

    # Results subsections
    results_raw = sections.get('2. Results', '').replace('## 2. Results', '').strip()
    res_subsections = re.split(r'\n(?=### 2\.\d+)', results_raw)

    # Process and weave figures and tables directly into results subsections
    woven_results = []
    for chunk in res_subsections:
        c_clean = chunk.strip()
        if not c_clean:
            continue
        first_line = c_clean.split('\n')[0]
        m_sec = re.search(r'###\s+2\.(\d+)', first_line)
        sec_num = int(m_sec.group(1)) if m_sec else 0
        
        # Remap figure references in text to match sequential figure environments
        c_remapped = remap_figure_references_for_section(c_clean, sec_num)
        chunk_latex = format_text_for_latex(c_remapped, key_map)
        woven_results.append(chunk_latex)
        
        # Attach figures and tables inline in strictly monotonic sequential order
        if sec_num == 1:
            woven_results.append("\n" + get_table_code(1) + "\n" + get_figure_code(1) + "\n")
        elif sec_num == 2:
            woven_results.append("\n" + get_figure_code(2) + "\n")
        elif sec_num == 3:
            woven_results.append("\n" + get_table_code(2) + "\n" + get_figure_code(3) + "\n")
        elif sec_num == 4:
            woven_results.append("\n" + get_table_code(3) + "\n" + get_figure_code(4) + "\n")
        elif sec_num == 5:
            woven_results.append("\n" + get_figure_code(5) + "\n")
        elif sec_num == 8:
            woven_results.append("\n" + get_figure_code(6) + "\n")
        elif sec_num == 9:
            woven_results.append("\n" + get_table_code(4) + "\n" + get_figure_code(7) + "\n")
        elif sec_num == 10:
            woven_results.append("\n" + get_table_code(5) + "\n" + get_figure_code(8) + "\n")
        elif sec_num == 11:
            woven_results.append("\n" + get_table_code(6) + "\n" + get_figure_code(9) + "\n")
        elif sec_num == 12:
            woven_results.append("\n" + get_table_code(7) + "\n" + get_figure_code(10) + "\n")
        elif sec_num == 14:
            woven_results.append("\n" + get_table_code(8) + "\n" + get_figure_code(11) + "\n")

    results_latex = "\n\n".join(woven_results)

    # Discussion
    disc_raw = sections.get('3. Discussion', '').replace('## 3. Discussion', '').strip()
    disc_latex = format_text_for_latex(disc_raw, key_map)

    # Online Methods
    methods_raw = sections.get('4. Online Methods', '').replace('## 4. Online Methods', '').strip()
    methods_latex = format_text_for_latex(methods_raw, key_map)

    # Author contributions
    contrib_raw = sections.get('Author Contributions', '').replace('## Author Contributions', '').strip()
    contrib_latex = format_text_for_latex(contrib_raw, key_map)

    tex = r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=1.8cm,top=2.0cm,bottom=2.0cm]{geometry}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{xcolor}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{subcaption}
\usepackage{enumitem}
\setlist{noitemsep,topsep=3pt,parsep=1pt}
\usepackage{lineno}
\usepackage[numbers,super,sort&compress]{natbib}
\usepackage[colorlinks=true,linkcolor=blue!70!black,citecolor=blue!70!black,urlcolor=blue!70!black]{hyperref}
\usepackage[protrusion=true,expansion=true]{microtype}
\usepackage{float}
\usepackage{url}
\usepackage{caption}

% Improve hyphenation and margin alignment
\emergencystretch=2em

% Figure path
\graphicspath{{figures/}{../docs/assets/}}

\title{\textbf{\Large Cross-Study Meta-Analysis of the Gut-Microbiota-Microglia Axis Uncovers Cell-Intrinsic Interferon Shutoff, Invariant Nutrient-Sensing Adapters, and Multi-Omic Reversibility}}

\author{
  \textbf{Samyak Meshram}$^{1,*}$, \textbf{Dr. Soumya Dhokey}$^{1}$ \\[0.4em]
  \small $^{1}$Independent Researchers, Open-Science Initiative \\
  \small $^*$Correspondence: \texttt{sam.meshram18@gmail.com} \quad $|$ \quad GitHub: \url{https://github.com/samyakmeshram/NeuroGut-MetaSeq}
}

\date{\small \today}

\begin{document}

\maketitle

% Start line numbering after title header for reviewer commentary
\linenumbers

\begin{abstract}
''' + abstract_latex + r'''
\end{abstract}

\vspace{0.8em}
\noindent\textbf{Significance Statement:}\newline
''' + significance_latex + r'''

\vspace{1.5em}
\section{Introduction}
''' + intro_latex + r'''

\vspace{1.5em}
\section{Results}
''' + results_latex + r'''

\clearpage
\section{Discussion}
''' + disc_latex + r'''

\vspace{1.8em}
\section{Online Methods}
\sloppy
''' + methods_latex + r'''

\vspace{1.5em}
\section*{Declarations}

\subsection*{Author Contributions}
''' + contrib_latex + r'''

\subsection*{Competing Interests}
The authors declare no competing interests.

\subsection*{Data Availability}
All analyzed raw transcriptomic and epigenomic datasets are publicly available in the NCBI Gene Expression Omnibus (GEO) under accession numbers \textbf{GSE107925}, \textbf{GSE108045}, \textbf{GSE186210}, \textbf{GSE266602}, \textbf{GSE152865}, and \textbf{GSE64977}. All standardized sample metadata, normalized read count matrices, deconvolution matrices, and meta-analytic summary tables are deposited and accessible in the study GitHub repository (\url{https://github.com/samyakmeshram/NeuroGut-MetaSeq}).

\subsection*{Code Availability}
All analytical pipelines, statistical modeling scripts, and visualization code are publicly available under the open-source MIT license at \url{https://github.com/samyakmeshram/NeuroGut-MetaSeq}. The analytical workflows are fully reproducible in Python 3.9+ and R 4.3+.

\vspace{1.8em}
\nocite{*}
\bibliographystyle{unsrtnat}
\bibliography{references}

\end{document}
'''
    with open("paper/manuscript.tex", "w") as f:
        f.write(tex)
    print("Wrote paper/manuscript.tex with inline figures, tables, and continuous flow")

def build_supplementary_tex():
    tex = r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=1.8cm,top=2.0cm,bottom=2.0cm]{geometry}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{xcolor}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{subcaption}
\usepackage{enumitem}
\setlist{noitemsep,topsep=3pt,parsep=1pt}
\usepackage{lineno}
\usepackage[numbers,super,sort&compress]{natbib}
\usepackage[colorlinks=true,linkcolor=blue!70!black,citecolor=blue!70!black,urlcolor=blue!70!black]{hyperref}
\usepackage[protrusion=true,expansion=true]{microtype}
\usepackage{float}
\usepackage{url}
\usepackage{caption}

% Figure path
\graphicspath{{figures/}{../docs/assets/}}

% Change caption and numbering prefixes to Supplementary (Fig S1, Tab S1)
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\thepage}{S\arabic{page}}

\title{\textbf{\Large Supplementary Information: Cross-Study Meta-Analysis of the Gut-Microbiota-Microglia Axis Uncovers Cell-Intrinsic Interferon Shutoff, Invariant Nutrient-Sensing Adapters, and Multi-Omic Reversibility}}

\author{
  \textbf{Samyak Meshram}$^{1,*}$, \textbf{Dr. Soumya Dhokey}$^{1}$ \\[0.4em]
  \small $^{1}$Independent Researchers, Open-Science Initiative \\
  \small $^*$Correspondence: \texttt{sam.meshram18@gmail.com} \quad $|$ \quad GitHub: \url{https://github.com/samyakmeshram/NeuroGut-MetaSeq}
}

\date{\small \today}

\begin{document}

\maketitle

\linenumbers

\section*{Supplementary Information Overview}
This document contains extended methodological formulations, supplementary diagnostic figures (Figures S1--S6), and supplementary statistical tables (Tables S1--S4) supporting the main text findings of \textbf{NeuroGut-MetaSeq (v1.3.0)}.

\tableofcontents

\vspace{1.5em}
\section{Extended Supplementary Methods}

\subsection{Small-Sample Random-Effects Derivation: REML and Hartung-Knapp-Sidik-Jonkman}
In classical meta-analysis, the DerSimonian-Laird (DL) method of moments estimator is widely used due to computational simplicity:
\begin{equation}
\hat{\tau}^2_{\text{DL}} = \max\left(0, \frac{Q - (K - 1)}{\sum w_k - \frac{\sum w_k^2}{\sum w_k}}\right)
\end{equation}
However, when the number of studies is small ($K < 5$), DL significantly underestimates between-study variance $\tau^2$, producing artificially narrow confidence intervals and inflated type I error rates. To overcome this limitation, NeuroGut-MetaSeq implements Restricted Maximum Likelihood (REML) variance estimation:
\begin{equation}
\ln L_{\text{REML}}(\tau^2) = -\frac{1}{2}\sum_{k=1}^K \ln(v_k + \tau^2) - \frac{1}{2}\ln\left(\sum_{k=1}^K \frac{1}{v_k + \tau^2}\right) - \frac{1}{2}\sum_{k=1}^K \frac{(\hat{\theta}_k - \hat{\theta}_{\text{RE}})^2}{v_k + \tau^2}
\end{equation}
The global maximum $\hat{\tau}^2_{\text{REML}}$ is found numerically via bounded Nelder-Mead optimization on $[0, \infty)$. Inverse-variance random-effects weights are then calculated as:
\begin{equation}
w_k^* = \frac{1}{v_k + \hat{\tau}^2_{\text{REML}}}, \quad \hat{\theta}_{\text{RE}} = \frac{\sum_{k=1}^K w_k^* \hat{\theta}_k}{\sum_{k=1}^K w_k^*}
\end{equation}
To ensure nominal 95\% coverage under small degrees of freedom ($K - 1 = 3$), standard errors are computed using the Hartung-Knapp-Sidik-Jonkman (HKSJ) adjustment:
\begin{equation}
q^* = \frac{1}{K - 1}\sum_{k=1}^K w_k^* (\hat{\theta}_k - \hat{\theta}_{\text{RE}})^2, \quad q_{\text{adj}}^* = \max(1.0, q^*)
\end{equation}
\begin{equation}
\operatorname{SE}_{\text{HKSJ}}(\hat{\theta}_{\text{RE}}) = \sqrt{\frac{q_{\text{adj}}^*}{\sum_{k=1}^K w_k^*}}
\end{equation}
Statistical significance and confidence intervals are computed using the Student's $t$-distribution critical value with $K - 1 = 3$ degrees of freedom ($t_{0.975, 3} = 3.1824$):
\begin{equation}
95\% \text{ CI}_{\text{HKSJ}} = \left[ \hat{\theta}_{\text{RE}} - 3.1824 \times \operatorname{SE}_{\text{HKSJ}}, \; \hat{\theta}_{\text{RE}} + 3.1824 \times \operatorname{SE}_{\text{HKSJ}} \right]
\end{equation}

\subsection{BayesPrism Ridge-Regularized Quadratic Programming Formulation}
To deconvolve microglial cell states without assuming mutually orthogonal basis vectors, we formulated Ridge-regularized constrained quadratic programming:
\begin{equation}
\min_{\mathbf{p}_j} \left\{ \frac{1}{2}\mathbf{p}_j^T (\mathbf{S}^T \mathbf{S} + \lambda \mathbf{I}) \mathbf{p}_j - (\mathbf{S}^T \mathbf{y}_j)^T \mathbf{p}_j \right\}
\end{equation}
subject to the biological constraints:
\begin{equation}
\sum_{s=1}^5 p_{js} = 1.0, \quad p_{js} \ge 0 \quad \forall s \in \{1,\dots,5\}
\end{equation}
Collinearity was evaluated via singular value decomposition $\mathbf{S} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T$, defining condition number $\kappa = \sigma_{\max} / \sigma_{\min}$.

\clearpage
\section{Supplementary Tables}

\begin{table}[htbp]
\centering
\caption{\textbf{Complete Quality Control and Alignment Metrics Across 60 Biological Samples.}}
\label{tab:s1_qc}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llllccp{5.5cm}}
\toprule
\textbf{Cohort} & \textbf{Condition} & \textbf{Library Depth Range} & \textbf{Mean Depth} & \textbf{Microglia Purity} & \textbf{Stress $Z$-Score} & \textbf{Sorting Technology} \\
\midrule
GSE107925 & SPF Control ($n=12$) & 18.4M -- 42.1M & 28.5M & $99.4\% \pm 0.2\%$ & $-0.08 \pm 0.31$ & FACS CD11b$^{+}$CD45$^{\text{low}}$ \\
GSE107925 & Germ-Free ($n=13$) & 19.1M -- 44.8M & 29.2M & $99.5\% \pm 0.2\%$ & $+0.06 \pm 0.29$ & FACS CD11b$^{+}$CD45$^{\text{low}}$ \\
GSE108045 & Vehicle ($n=6$) & 14.2M -- 31.0M & 22.4M & $99.2\% \pm 0.3\%$ & $-0.12 \pm 0.25$ & FACS CD11b$^{+}$CD45$^{\text{low}}$ \\
GSE108045 & Antibiotics ($n=6$) & 15.6M -- 33.5M & 23.8M & $99.3\% \pm 0.2\%$ & $+0.09 \pm 0.28$ & FACS CD11b$^{+}$CD45$^{\text{low}}$ \\
GSE186210 & Standard Fiber ($n=7$) & 22.5M -- 54.2M & 36.1M & $97.8\% \pm 0.5\%$ & $+0.04 \pm 0.35$ & MACS CD11b Microbeads \\
GSE186210 & Zero-Fiber ($n=7$) & 21.0M -- 51.8M & 34.8M & $97.9\% \pm 0.4\%$ & $-0.03 \pm 0.33$ & MACS CD11b Microbeads \\
GSE266602 & Colonized ($n=4$) & 6.8M -- 14.5M & 10.2M & $95.1\% \pm 1.1\%$ & $+0.15 \pm 0.41$ & Percoll Density Gradient \\
GSE266602 & Germ-Free ($n=5$) & 7.2M -- 15.1M & 10.8M & $94.8\% \pm 1.2\%$ & $+0.11 \pm 0.38$ & Percoll Density Gradient \\
\bottomrule
\end{tabular}%
}
\end{table}

\vspace{1.5em}

\begin{table}[htbp]
\centering
\caption{\textbf{Top Consensus Significant Meta-DEGs under REML-HKSJ Pooling Across 23,096 Genes.}}
\label{tab:s2_top_meta}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lcccccl}
\toprule
\textbf{Gene Symbol} & \textbf{$k$} & \textbf{REML $\log_2\text{FC}$} & \textbf{HKSJ 95\% CI} & \textbf{$p_{\text{HKSJ}}$} & \textbf{$I^2$} & \textbf{Consensus Designation} \\
\midrule
\textit{Llgl2} & 4 & +0.672 & [0.166, 1.178] & 0.0243 & 0.0\% & Invariant Basolateral Core \\
\textit{Slfn2} & 4 & -0.465 & [-0.789, -0.141] & 0.0197 & 9.6\% & Invariant Quiescence Guardian \\
\textit{Clu} & 3 & +0.850 & [0.251, 1.448] & 0.0215 & 0.0\% & Invariant Extracellular Chaperone \\
\textit{Fosb} & 3 & +1.345 & [0.797, 1.892] & 0.0160 & 0.0\% & Invariant Immediate Early Activator \\
\textit{Sap30} & 3 & -0.393 & [-0.653, -0.131] & 0.0409 & 0.0\% & Invariant Sin3A-HDAC Corepressor \\
\textit{Card6} & 3 & -0.431 & [-0.715, -0.146] & 0.0398 & 0.0\% & Invariant NF-$\kappa$B / NOD Modulator \\
\textit{Hes1} & 3 & -0.365 & [-0.684, -0.046] & 0.0412 & 0.0\% & Notch Pathway Transcriptional Brake \\
\textit{Tnf} & 3 & +1.038 & [-0.183, 2.259] & 0.0821 & 32.2\% & Conserved Low-Grade Priming \\
\bottomrule
\end{tabular}%
}
\end{table}

\vspace{1.5em}

\begin{table}[htbp]
\centering
\caption{\textbf{Factorial Interaction Meta-Regression Metrics for Sex Dimorphism ($N=51$).}}
\label{tab:s3_sex}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lcccccclp{6.5cm}}
\toprule
\textbf{Gene} & \textbf{Male LFC} & \textbf{Female LFC} & \textbf{$\hat{\theta}_{\text{int}}$} & \textbf{$p_{\text{int}}$} & \textbf{$I^2_{\text{sex}}$} & \textbf{Dimorphism Status} & \textbf{Biological Interpretation} \\
\midrule
\textit{Slfn2} & -0.521 & -0.183 & +0.3382 & 0.017 & 0.0\% & Male-Biased Vulnerability & Steeper loss of quiescence in males \\
\textit{Oas1a} & -0.489 & -0.154 & +0.3348 & 0.051 & 0.0\% & Male-Biased Vulnerability & Antiviral surveillance collapse \\
\textit{Irf1} & -0.061 & -0.079 & -0.0177 & 0.863 & 0.0\% & Sex-Shared (Invariant) & Master regulon shutoff is sex-invariant \\
\textit{Stat1} & -0.285 & -0.377 & -0.0918 & 0.362 & 0.0\% & Sex-Shared (Invariant) & Downstream JAK/STAT shutdown \\
\textit{Llgl2} & +0.076 & +0.051 & -0.0247 & 0.941 & 35.8\% & Sex-Shared (Invariant) & Nutrient adaptation is sex-invariant \\
\textit{Clu} & +0.412 & +0.380 & -0.0321 & 0.786 & 0.0\% & Sex-Shared (Invariant) & Chaperone secretion is sex-invariant \\
\bottomrule
\end{tabular}%
}
\end{table}

\vspace{1.5em}

\begin{table}[htbp]
\centering
\caption{\textbf{Top TRRUST Transcription Factor Regulon Enrichments Across Common Transcriptome.}}
\label{tab:s4_trrust}
\resizebox{\textwidth}{!}{%
\begin{tabular}{clccccl}
\toprule
\textbf{Rank} & \textbf{Regulon TF} & \textbf{Targets ($n$)} & \textbf{Mean Target LFC} & \textbf{$Z_{\text{activity}}$} & \textbf{FDR} & \textbf{Regulon Response} \\
\midrule
1 & \textbf{IRF1} & 23 & -0.208 & -2.284 & 0.0384 & \textbf{Significant Master Repression} \\
2 & STAT1 & 41 & -0.145 & -1.821 & 0.0892 & Repressed Interferon Cascade \\
3 & RELA (NF-$\kappa$B) & 82 & +0.092 & +1.341 & 0.1820 & Modest Pro-inflammatory Priming \\
4 & JUN (AP-1) & 54 & +0.114 & +1.520 & 0.1540 & Immediate Early Priming \\
5 & FOS (AP-1) & 38 & +0.142 & +1.684 & 0.1210 & Immediate Early Priming \\
6 & SP1 & 95 & -0.042 & -0.612 & 0.5410 & Stable Basal Transcription \\
\bottomrule
\end{tabular}%
}
\end{table}

\clearpage
\section{Supplementary Figures}

% Figure S1
\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_qc_microglial_purity.png}
  \caption{Microglial Lineage Purity Audit}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_qc_library_depths.png}
  \caption{Library Depth Distributions}
\end{subfigure}
\caption{\textbf{Detailed Quality Control Audit and Microglial Lineage Validation.} 
(A) Assessment of microglial purity by measuring canonical microglial marker expression versus non-microglial contaminants across 60 biological libraries. 
(B) Sequencing depth distributions across individual experimental cohorts.}
\label{fig:figs1}
\end{figure}

\vspace{1.5em}

% Figure S2
\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_lfc_correlation_heatmap.png}
  \caption{Cross-Cohort Log2FC Correlation}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_cohort_deg_overlap.png}
  \caption{Cohort DEG Overlap}
\end{subfigure}
\caption{\textbf{Study-Specific Variance and Unpooled DEG Concordance.} 
(A) Pairwise Pearson correlation heatmap of genome-wide effect sizes between cohorts, showing uncoupled baseline laboratory noise ($\rho \approx 0$). 
(B) Overlap of differentially expressed genes between cohorts illustrating limited raw gene intersection without statistical pooling.}
\label{fig:figs2}
\end{figure}

\clearpage

% Figure S3
\begin{figure}[htbp]
\centering
\includegraphics[width=0.82\textwidth]{figures/fig_heterogeneity_distribution.png}
\caption{\textbf{Distribution of Higgins $I^2$ Heterogeneity Across 23,096 Evaluated Genes.} 
Histogram of Higgins $I^2$ values across the common transcriptome. Over 70.5\% of genes exhibit low between-study heterogeneity ($I^2 < 25\%$), with only 4.5\% displaying extreme heterogeneity ($I^2 > 75\%$).}
\label{fig:figs3}
\end{figure}

\vspace{1.5em}

% Figure S4
\begin{figure}[htbp]
\centering
\includegraphics[width=0.88\textwidth]{figures/fig_loo_stability.png}
\caption{\textbf{Systematic Leave-One-Out (LOO) Sensitivity Meta-Analyses.} 
Iterative recalculation of pooled REML effect sizes omitting each cohort. Core consensus loci (\textit{Llgl2}, \textit{Clu}, \textit{Slfn2}, \textit{Sap30}) remain tightly centered along the unity line across all iterations, while acute shock effectors (\textit{Tsc22d3}, \textit{Ddit4}) shift to zero upon omission of the antibiotic cohort.}
\label{fig:figs4}
\end{figure}

\clearpage

% Figure S5
\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_wgcna_modules_eigengenes.png}
  \caption{WGCNA Module-Trait Correlations}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_network_hub_subgraph.png}
  \caption{Topological Overlap Hub Network}
\end{subfigure}
\caption{\textbf{Weighted Gene Co-Expression Networks and Consensus Hub Subgraphs.} 
(A) WGCNA module eigengenes and their correlations with experimental microbiome depletion. 
(B) Consensus co-expression hub subgraph highlighting topological overlap connections centered on \textit{Llgl2} and \textit{Clu}.}
\label{fig:figs5}
\end{figure}

\vspace{1.5em}

% Figure S6
\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_scfa_rescue_specificity_null.png}
  \caption{Specificity Null Distribution}
\end{subfigure}
\hfill
\begin{subfigure}[b]{0.48\textwidth}
  \includegraphics[width=\textwidth]{figures/fig_scfa_rescue_inversion.png}
  \caption{SCFA Signature Inversion}
\end{subfigure}
\caption{\textbf{SCFA Metabolite Rescue Specificity Null Model Permutations and Trajectory.} 
(A) 1,000-iteration permutation test against background non-DEGs ($p_{\text{perm}} < 0.001$), confirming specificity of the rescue index. 
(B) Scatter plot illustrating the profound reciprocal correlation ($r = -0.873$) between meta-analytic depletion and in vivo SCFA repletion.}
\label{fig:figs6}
\end{figure}

\end{document}
'''
    with open("paper/supplementary.tex", "w") as f:
        f.write(tex)
    print("Wrote paper/supplementary.tex with resizebox tables and clean layout")

def build_makefile_and_build_script():
    makefile = r'''PDFLATEX = /Library/TeX/texbin/pdflatex
BIBTEX = /Library/TeX/texbin/bibtex

all: manuscript.pdf supplementary.pdf

manuscript.pdf: manuscript.tex references.bib
	$(PDFLATEX) -interaction=nonstopmode manuscript.tex
	$(BIBTEX) manuscript
	$(PDFLATEX) -interaction=nonstopmode manuscript.tex
	$(PDFLATEX) -interaction=nonstopmode manuscript.tex

supplementary.pdf: supplementary.tex references.bib
	$(PDFLATEX) -interaction=nonstopmode supplementary.tex
	$(BIBTEX) supplementary || true
	$(PDFLATEX) -interaction=nonstopmode supplementary.tex
	$(PDFLATEX) -interaction=nonstopmode supplementary.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lot *.lof *.synctex.gz *.pdf

.PHONY: all clean
'''
    with open("paper/Makefile", "w") as f:
        f.write(makefile)
    print("Wrote paper/Makefile")

    build_sh = r'''#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"

PDFLATEX="/Library/TeX/texbin/pdflatex"
BIBTEX="/Library/TeX/texbin/bibtex"

echo "=== Compiling paper/manuscript.pdf ==="
$PDFLATEX -interaction=nonstopmode manuscript.tex
$BIBTEX manuscript
$PDFLATEX -interaction=nonstopmode manuscript.tex
$PDFLATEX -interaction=nonstopmode manuscript.tex

echo "=== Compiling paper/supplementary.pdf ==="
$PDFLATEX -interaction=nonstopmode supplementary.tex
$BIBTEX supplementary || true
$PDFLATEX -interaction=nonstopmode supplementary.tex
$PDFLATEX -interaction=nonstopmode supplementary.tex

echo "=== LaTeX Compilation Completed Successfully! ==="
ls -lh manuscript.pdf supplementary.pdf
'''
    with open("paper/build.sh", "w") as f:
        f.write(build_sh)
    os.chmod("paper/build.sh", 0o755)
    print("Wrote paper/build.sh")

if __name__ == '__main__':
    build_manuscript_tex()
    build_supplementary_tex()
    build_makefile_and_build_script()
    print("All LaTeX sources and build scripts created successfully.")
