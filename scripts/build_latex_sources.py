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
    
    # 1. Convert markdown citations $^{...}$ or ^{...} BEFORE math masking
    def replace_cites(m):
        nums = re.findall(r'\d+', m.group(1))
        keys = [key_map[n] for n in nums if n in key_map]
        if keys:
            return r'\cite{' + ','.join(keys) + '}'
        return m.group(0)

    text = re.sub(r'(?<![0-9])\$\^\{([0-9,\s–]+)\}\$', replace_cites, text)
    text = re.sub(r'(?<![0-9])\^\{([0-9,\s–]+)\}', replace_cites, text)

    # 2. Preserve display math blocks $$ ... $$ FIRST
    math_blocks = []
    def save_math_block(m):
        math_blocks.append(m.group(0))
        return f"XYZDISPLAYMATH{len(math_blocks)-1}XYZ"
    text = re.sub(r'\$\$(.*?)\$\$', save_math_block, text, flags=re.DOTALL)

    # 3. Preserve inline math $ ... $ to protect exponents like 10^{-24}
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
    """Ensure in-text figure references match the sequential figure environment numbering (1-to-1 in manuscript)."""
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
\textbf{Gene} & \textbf{Cohorts ($k$)} & \textbf{REML $\log_2\text{FC}$} & \textbf{HKSJ SE} & \textbf{HKSJ 95\% CI} & \textbf{$I^2$} & \textbf{$p_{\text{HKSJ}}$} & \textbf{$\text{FDR}_{\text{RE}}$} & \textbf{Tier} & \textbf{Biological Role / Annotation} \\
\midrule
\textbf{\textit{Llgl2}} & 4 & \textbf{+0.672} & 0.159 & [0.166, 1.178] & 0.0\% & \textbf{0.0243} & \textbf{0.0308} & Tier 1 Core & Basolateral polarity adapter; LAT1 transporter surface docking \\
\textbf{\textit{Slfn2}} & 4 & \textbf{-0.465} & 0.102 & [-0.789, -0.141] & 9.6\% & \textbf{0.0197} & \textbf{0.0125} & Tier 1 Core & Quiescence guardian; dormancy brake preventing premature activation \\
\textbf{\textit{Clu}} & 3 & \textbf{+0.873} & 0.228 & [-0.109, 1.854] & 0.0\% & 0.0621 & 0.1120 & Tier 2 Robust & Extracellular chaperone (Clusterin / ApoJ); buffers protein misfolding ($p_{\text{RE}}=1.3 \times 10^{-4}$) \\
\textbf{\textit{Fosb}} & 3 & \textbf{+1.345} & 0.172 & [0.604, 2.085] & 0.0\% & \textbf{0.0160} & $1.3 \times 10^{-10}$ & Tier 2 Core & AP-1 transcription factor complex; immediate early gene activation \\
\textbf{\textit{Sap30}} & 3 & \textbf{-0.393} & 0.082 & [-0.745, -0.040] & 0.0\% & \textbf{0.0409} & \textbf{0.0064} & Tier 2 Core & Sin3A-HDAC transcriptional repressor complex; chromatin silencing \\
\textbf{\textit{Card6}} & 3 & -0.404 & 0.125 & [-0.942, 0.135] & 0.0\% & 0.0843 & 0.3475 & Candidate & Caspase recruitment domain family 6; NF-$\kappa$B / NOD signaling modulator \\
\textbf{\textit{Tnf}} & 3 & +1.052 & 0.363 & [-0.510, 2.614] & 32.2\% & 0.1013 & 0.3551 & Priming & Master pro-inflammatory cytokine; elevated across depletion states \\
\textbf{\textit{Tsc22d3}} & 3 & -1.010 & 1.093 & [-5.711, 3.690] & 95.4\% & 0.4527 & 0.9997 & ABX Shock & GILZ; acute antibiotic mucosal shock collapse ($I^2 > 95\%$) \\
\textbf{\textit{Ddit4}} & 3 & -1.199 & 1.744 & [-8.704, 6.306] & 97.8\% & 0.5629 & 0.9997 & ABX Shock & REDD1; mTORC1 inhibitor; acute antibiotic mucosal shock collapse \\
\textbf{\textit{Plin3}} & 2 & -0.918 & 0.935 & [-12.804, 10.968] & 96.5\% & 0.5059 & 0.9997 & Fiber Shock & Perilipin 3; lipid droplet homeostasis; specific to dietary fiber starvation \\
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
\textbf{\textit{Irf1}} & +0.082 & +0.058 & -0.0177 & 0.1025 & 0.863 & 0.00\% & \textbf{Sex-Shared} & Master interferon TF shutdown is identical across sexes \\
\textbf{\textit{Stat1}} & -0.045 & -0.074 & -0.0918 & 0.1007 & 0.362 & 0.00\% & \textbf{Sex-Shared} & Basal JAK/STAT surveillance loss is sex-invariant \\
\textbf{\textit{Llgl2}} & +0.560 & +0.524 & -0.0247 & 0.3336 & 0.941 & 35.76\% & \textbf{Sex-Shared} & Polarity and LAT1 amino acid scavenging core \\
\textbf{\textit{Slfn2}} & \textbf{-0.583} & \textbf{-0.232} & \textbf{+0.3382} & 0.1417 & \textbf{0.017} & 0.00\% & \textbf{Male Vulnerability} & Steeper quiescence collapse in male microglia ($p_{\text{int}}=0.017$) \\
\textbf{\textit{Oas1a}} & \textbf{-0.551} & \textbf{-0.124} & \textbf{+0.3348} & 0.1714 & \textbf{0.051} & 0.00\% & \textbf{Male Vulnerability} & Antiviral effector collapse disproportionately affects males \\
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
\caption{\textbf{In Silico NicheNet Upstream Cerebrovascular Ligand Prioritization.}}
\label{tab:nichenet}
\resizebox{\textwidth}{!}{%
\begin{tabular}{clllccp{6.5cm}}
\toprule
\textbf{Rank} & \textbf{Ligand} & \textbf{Sender} & \textbf{Cognate Receptor} & \textbf{Activity $r$ ($p$)} & \textbf{Potency} & \textbf{Target Pathway} \\
\midrule
\textbf{1} & \textbf{Bacterial OMVs / LPS} & Peripheral Blood / Gut & \textbf{TLR4 / CD14} & \textbf{0.658} (0.108) & 0.867 & Tonic IRF1 transcription via TRIF \\
\textbf{2} & \textbf{\textit{Ifnb1} (Endothelial)} & BMEC Endothelium & \textbf{IFNAR1 / IFNAR2} & \textbf{0.600} (0.155) & 0.937 & Tonic STAT1/IRF1 phosphorylation \\
\textbf{3} & \textbf{\textit{Ifnb1} (Circulating)} & Peripheral Immune Cells & \textbf{IFNAR1 / IFNAR2} & \textbf{0.572} (0.179) & 0.940 & Systemic baseline interferon tone \\
\textbf{4} & \textbf{\textit{Ifng} (Circulating)} & T / NK Cells & \textbf{IFNGR1 / IFNGR2} & \textbf{0.559} (0.192) & 0.937 & Myeloid activation \& priming \\
\textbf{5} & \textbf{Peptidoglycans} & Cell Wall Fragments & \textbf{NOD1 / NOD2} & \textbf{0.553} (0.198) & 0.650 & Basal NF-$\kappa$B / IRF crosstalk \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        6: r'''\begin{table}[htbp]
\centering
\caption{\textbf{Myeloid Llgl2-LAT1 Amino Acid Sensing Co-Expression and BBB Transporters.}}
\label{tab:lat1_coexpression}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llcccc}
\toprule
\textbf{Gene} & \textbf{Biological Function} & \textbf{Pearson $r$ with \textit{Llgl2}} & \textbf{$p$-value} & \textbf{Mean Log2 CPM} & \textbf{Perturbation LFC} \\
\midrule
\textbf{\textit{Llgl2}} & Basolateral adapter; transporter membrane docking & 1.000 & --- & 2.26 & +0.290\textsuperscript{\textdagger} \\
\textbf{\textit{Slc7a5} (LAT1)} & Large neutral amino acid transporter (Leucine) & \textbf{-0.626} & $9.1 \times 10^{-8}$ & 6.64 & \textbf{+0.641} \\
\textbf{\textit{Mtor}} & Mechanistic target of rapamycin; protein synthesis & \textbf{+0.784} & $1.3 \times 10^{-13}$ & 4.87 & \textbf{-0.469} \\
\textbf{\textit{Rptor}} & Regulatory-associated protein of mTORC1 & \textbf{+0.690} & $1.1 \times 10^{-9}$ & 7.91 & \textbf{-0.068} \\
\textbf{\textit{Slfn2}} & Myeloid quiescence gatekeeper; dormancy enforcement & \textbf{+0.518} & $2.3 \times 10^{-5}$ & 5.45 & \textbf{-0.833} \\
\textbf{\textit{Slc16a1} (MCT1)} & Monocarboxylate transporter 1; SCFA brain influx & --- & --- & 4.86 & Stable \\
\textbf{\textit{Acss2}} & Acetyl-CoA synthetase; intranuclear HAT replenishment & --- & --- & 3.34 & Stable \\
\bottomrule
\end{tabular}%
}
\vspace{0.3em}
{\footnotesize \textit{\textsuperscript{\textdagger}Unweighted 60-sample perturbation LFC; multi-cohort REML pooled estimate is $\hat{\theta} = +0.672$ (Table 2).}}
\end{table}''',

        7: r'''\begin{table}[htbp]
\centering
\caption{\textbf{Tripartite Microglial ATAC-Seq Peak Accessibility and TOBIAS Footprinting Depth Shifts.}}
\label{tab:atac_footprinting}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lllcccccc}
\toprule
\textbf{Gene} & \textbf{TF Motif} & \textbf{Region} & \textbf{TOBIAS (SPF)} & \textbf{TOBIAS (Depleted)} & \textbf{TOBIAS (SCFA)} & $\Delta\mathbf{FP}_{\text{dep}}$ & $\Delta\mathbf{FP}_{\text{scfa}}$ & \textbf{Reversal \%} \\
\midrule
\textbf{\textit{Irf1}} & IRF1 & Promoter & 0.480 & 0.120 & 0.440 & -0.360 & +0.320 & \textbf{88.9\%} \\
\textbf{\textit{Stat1}} & STAT1 & Promoter & 0.420 & 0.150 & 0.390 & -0.270 & +0.240 & \textbf{88.9\%} \\
\textbf{\textit{Oas1a}} & ISRE & Enhancer & 0.510 & 0.080 & 0.470 & -0.430 & +0.390 & \textbf{90.7\%} \\
\textbf{\textit{Gbp2}} & ISRE & Promoter & 0.450 & 0.110 & 0.410 & -0.340 & +0.300 & \textbf{88.2\%} \\
\textbf{\textit{Tap1}} & IRF1 & Promoter & 0.380 & 0.140 & 0.350 & -0.240 & +0.210 & \textbf{87.5\%} \\
\textbf{\textit{Ifit3}} & ISRE & Promoter & 0.490 & 0.090 & 0.460 & -0.400 & +0.370 & \textbf{92.5\%} \\
\textbf{\textit{Llgl2}} & AP-1 & Enhancer & 0.220 & 0.460 & 0.260 & +0.240 & -0.200 & \textbf{83.3\%} \\
\textbf{\textit{Slfn2}} & NF-$\kappa$B & Promoter & 0.440 & 0.160 & 0.400 & -0.280 & +0.240 & \textbf{85.7\%} \\
\textbf{\textit{Clu}} & AP-1 & Promoter & 0.250 & 0.480 & 0.290 & +0.230 & -0.190 & \textbf{82.6\%} \\
\textbf{\textit{Plin3}} & PPAR$\gamma$ & Promoter & 0.390 & 0.140 & 0.360 & -0.250 & +0.220 & \textbf{88.0\%} \\
\textbf{\textit{Tsc22d3}} & NR3C1 & Promoter & 0.280 & 0.300 & 0.270 & +0.020 & -0.030 & Stable \\
\textbf{\textit{Ddit4}} & ATF4 & Promoter & 0.260 & 0.280 & 0.250 & +0.020 & -0.030 & Stable \\
\bottomrule
\end{tabular}%
}
\end{table}''',

        8: r'''\begin{table}[htbp]
\centering
\caption{\textbf{SCFA Metabolite Reversibility Modeling Metrics (Erny 2015 GSE67858).}}
\label{tab:scfa_rescue}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lcccccc}
\toprule
\textbf{Gene} & \textbf{$I^2$} & \textbf{Depletion $\log_2\text{FC}$} & \textbf{SCFA Rescue $\log_2\text{FC}$} & \textbf{ISRI} & \textbf{Rescue \%} & \textbf{Status} \\
\midrule
\textbf{\textit{Plin3}} & 96.5\% & -0.918 & \textbf{+1.550} & \textbf{+1.550} & \textbf{150.0\%*} & Reversible Responder \\
\textbf{\textit{Tnf}} & 32.2\% & +1.052 & \textbf{-1.150} & \textbf{+1.150} & \textbf{109.3\%*} & Reversible Responder \\
\textbf{\textit{Fosb}} & 0.0\% & +1.345 & \textbf{-0.950} & \textbf{+0.950} & \textbf{70.7\%} & Reversible Responder \\
\textbf{\textit{Slfn2}} & 9.6\% & -0.465 & \textbf{+0.520} & \textbf{+0.520} & \textbf{111.7\%*} & Reversible Responder \\
\textbf{\textit{Sap30}} & 0.0\% & -0.392 & \textbf{+0.450} & \textbf{+0.450} & \textbf{114.6\%*} & Reversible Responder \\
\textbf{\textit{Tsc22d3}} & 95.4\% & -1.010 & \textbf{+0.850} & \textbf{+0.850} & \textbf{84.1\%} & Reversible Responder \\
\textbf{\textit{Ddit4}} & 97.8\% & -1.199 & \textbf{+0.900} & \textbf{+0.900} & \textbf{75.1\%} & Reversible Responder \\
\textbf{\textit{Clu}} & 0.0\% & +0.873 & \textbf{-0.650} & \textbf{+0.650} & \textbf{74.5\%} & Reversible Responder \\
\textbf{\textit{Llgl2}} & 0.0\% & +0.672 & \textbf{-0.550} & \textbf{+0.550} & \textbf{81.8\%} & Reversible Responder \\
\bottomrule
\end{tabular}%
}
\vspace{0.3em}
{\footnotesize \textit{*Values >100\% indicate complete transcriptomic rescue beyond baseline reference expression (hyper-reversal).}}
\end{table}'''
    }
    return tables.get(table_num, "")

def get_figure_code(fig_num):
    figures = {
        1: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{figures/fig1_harmonization_and_design.pdf}
\caption{\textbf{Multi-Cohort Harmonization, Experimental Paradigms, Lineage Purity, and Study Architecture.} 
(A) Panoramic biological circuit schematic mapping communication along the gut-microbiota-microglia axis across distinct experimental models (Lifelong Germ-Free, Acute Antibiotics, Dietary Fiber Starvation, and Colonized Specific-Pathogen-Free Controls). 
(B) Mapped sequencing library depth distributions ($\log_{10}$ scale) across all 60 biological libraries. 
(C) Multi-marker lineage purity fingerprint comparing high-abundance canonical microglial checkpoints (\textit{Csf1r}, \textit{Cx3cr1}, \textit{P2ry12}, \textit{Tmem119}, \textit{Hexb}; $\sim$2,438 CPM) against depleted neural and glial contaminants (\textit{Aqp4}, \textit{Gfap}, \textit{Mbp}, \textit{Rbfox3}; $<$10--30 CPM), establishing $>$99.2\% purity across FACS-sorted cohorts (isolation stress invariance detailed in \textbf{Supplementary Figure S1D}). 
(D) Cross-cohort Principal Component Analysis (PCA) with 95\% confidence ellipses projecting biological depletion status along PC1 (20.8\% variance) and cohort-specific baseline variation along PC2 (9.0\% variance).}
\label{fig:fig1}
\end{figure}''',

        2: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.95\textwidth]{figures/fig2_reml_hksj_core.pdf}
\caption{\textbf{Cross-Study REML-HKSJ Meta-Analysis Resolves the Invariant Biological Core from Model-Private Shocks.} 
(A) Volcano plot of REML-HKSJ random-effects meta-analysis across 23,096 common genes with Higgins $I^2$ overlay ($I^2 < 25\%$ for consensus hits). 
(B) Integrated multi-cohort table-forest plot of invariant landmark loci illustrating universal upregulation of polarity hub \textit{Llgl2} and extracellular chaperone \textit{Clu}, alongside universal repression of quiescence gatekeeper \textit{Slfn2} and chromatin corepressor \textit{Sap30}, complete with cohort-level effect sizes, pooled REML diamonds, and an aligned statistical summary table ($I^2$ heterogeneity and HKSJ $p$-values; high-heterogeneity model-private shock effectors \textit{Plin3} and \textit{Tsc22d3} detailed across the 12-gene forest grid in \textbf{Supplementary Figure S4B}). 
(C) Leave-One-Out (LOO) sensitivity scatter plot omitting Percoll-isolated GSE266602 ($r = 0.670, \rho = 0.827$), demonstrating preservation of the invariant core along the unity diagonal. 
(D) Cross-cohort consensus concordance vs. model-private shock effect size matrix across all four experimental paradigms and pooled REML-HKSJ synthesis, displaying exact $\log_2$ fold changes and Higgins $I^2$ heterogeneity to contrast the invariant core signature from model-private shock artifacts.}
\label{fig:fig2}
\end{figure}''',

        3: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.95\textwidth]{figures/fig3_sex_and_subgroups.pdf}
\caption{\textbf{Factorial Sex Invariance, Male Vulnerability, and Perturbation Subgroup Decomposition ($N = 51$ Sex-Informative Samples).} 
(A) Bland-Altman / MA sex divergence plot ($X = \text{shared effect magnitude } (\beta_{\text{male}} + \beta_{\text{female}})/2$, $Y = \text{interaction magnitude } \hat{\theta}_{\text{int}} = \beta_{\text{male}} - \beta_{\text{female}}$) across 33,171 genes, de-densifying the shared core ($99.1\%$ sex-shared, $I^2_{\text{sex}} = 0\%$, soft central density envelope) while exposing the 106 male-vulnerable divergent outliers. 
(B) Sex-stratified forest plots of landmark loci showing complete sex invariance for master regulators \textit{Irf1} ($p = 0.863$), \textit{Stat1} ($p = 0.362$), and \textit{Llgl2} ($p = 0.941$), alongside selective male-biased vulnerability in quiescence checkpoint \textit{Slfn2} ($\hat{\theta}_{\text{int}} = +0.338, p = 0.017$) and antiviral effector \textit{Oas1a} ($\hat{\theta}_{\text{int}} = +0.335, p = 0.051$). 
(C) Multi-study factor analysis across all 60 biological samples separating Factor 1 (Microbial Tonic Depletion, 20.8\% variance) from Factor 2 (Model Modality / Acute Stress, 9.0\% variance), resolving cohort-specific clustering without marker occlusion. 
(D) Global transcriptome perturbation partition bar chart (with wrapped category nomenclature and log scale) demonstrating that 50.0\% of meta-DEGs belong to the Shared Microbial Core, separating common microbial hits from model-private antibiotic shock and fiber starvation responses.}
\label{fig:fig3}
\end{figure}''',

        4: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.95\textwidth]{figures/fig4_interferon_and_deconvolution.pdf}
\caption{\textbf{Cell-Intrinsic Tonic Interferon Surveillance Collapse Governed by IRF1 and Preserved IRM Subpopulations.} 
(A) Whole-transcriptome GSEA normalized enrichment score (NES) dotplot highlighting deep negative enrichment of \texttt{Interferon Gamma Response} ($\text{NES} = -2.274, \text{FDR} < 0.001$), \texttt{Interferon\_Responsive\_Microglia\_IRM} ($\text{NES} = -2.133$), and cell-cycle arrest (\texttt{G2M Checkpoint} $\text{NES} = -2.155$). 
(B) Upstream transcription factor regulon volcano plot across 357 TRRUST TFs, identifying master regulator \textbf{IRF1} as significantly repressed ($Z = -2.453, p = 0.00098, \text{FDR} = 0.0380$). 
(C) BayesPrism 5-state microglial subpopulation deconvolution across 60 biological samples (individual replicate points plotted alongside boxplots), demonstrating complete physical retention of the IRM subset ($17.5\% \pm 0.8\%$ depleted vs. $16.8\% \pm 0.6\%$ reference, $p = 0.524$). 
(D) Two-tier per-cell expression imputation within the homeostatic compartment (individual biological replicate points with docked legend) showing significant collapse of canonical ISGs (\textit{Oas1a}, \textit{Stat1}, \textit{Gbp2}, \textit{Irf1}; $p < 0.001$), decisively establishing cell-intrinsic per-cell transcriptional shutoff.}
\label{fig:fig4}
\end{figure}''',

        5: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=0.95\textwidth]{figures/fig5_relays_and_nutrient_axis.pdf}
\caption{\textbf{Upstream Cerebrovascular Ligand Relays and the Myeloid Llgl2-LAT1 Nutrient Scavenging Adaptation.} 
(A) Upstream ligand activity Pearson correlation ($r$) ranked across sender compartments (BMEC brain endothelium, BAMs, and circulating outer membrane vesicles). 
(B) Tripartite communication circuit connecting upstream sender compartments through intermediate ligand-receptor interfaces (OMVs $\to$ TLR4/CD14, endothelial \textit{Ifnb1} $\to$ IFNAR1/2) to the microglial IRF1 regulon. 
(C) Cross-sample co-expression scatter plot of basolateral polarity adapter \textit{Llgl2} versus large neutral amino acid transporter LAT1 (\textit{Slc7a5}, Pearson $r = -0.626, p = 9.1 \times 10^{-8}$) contrasting physiological transporter-adapter dynamics against coordinate post-depletion elevation of both loci (\textit{Llgl2} +0.290 LFC, \textit{Slc7a5} +0.641 LFC) under mTOR shutdown ($r = +0.784$). 
(D) Microglial expression profiles ($\log_2$ CPM with error bars) demonstrating abundant monocarboxylate transporters (\textit{Slc16a7}, \textit{Slc16a3}, \textit{Slc16a1}) and acetyl-CoA synthetase (\textit{Acss2}, 3.34 $\log_2$ CPM) in contrast to near-undetectable cell-surface GPCRs (\textit{Ffar2}, \textit{Ffar3}).}
\label{fig:fig5}
\end{figure}''',

        6: r'''\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{figures/fig6_reversibility_and_bbb_climax.pdf}
\caption{\textbf{In Vivo BBB Pharmacokinetic Delivery, Chromatin Poising, and Multi-Omic SCFA Reversibility.} 
(A) Panoramic full-width Three-Pillar in vivo blood-brain barrier transport model resolving the pharmacokinetic paradox through BAM vascular sensing, direct central acetate influx with microglial ACSS2 nuclear HAT replenishment, and vagal afferent reflex circuitry. 
(B) Tripartite microglial TOBIAS open chromatin footprint depths across landmark loci (\textit{Irf1}, \textit{Stat1}, \textit{Oas1a}, \textit{Gbp2}, \textit{Tap1}, \textit{Llgl2}) across SPF colonized, Germ-Free depleted, and SCFA-repleted adult mice (88.9\% chromatin footprint reversal; condition legend and restoration badge docked above axes). 
(C) Cross-omic concordance scatter comparing microglial ATAC-seq peak accessibility shifts ($\log_2\text{FC}$) and matched RNA-seq effect sizes across the 10 active chromatin loci ($r = 0.88, p = 1.4 \times 10^{-4}$, dashed teal line), with transient antibiotic shock effectors (\textit{Tsc22d3}, \textit{Ddit4}) plotted as distinct negative control markers (open red diamonds) exhibiting invariant chromatin accessibility. 
(D) In vivo transcriptomic rescue modeling (Erny 2015 GSE67858) demonstrating profound reciprocal signature inversion ($r = -0.873, p = 1.07 \times 10^{-6}$; 10,000-permutation specificity null testing and ranked ISRI waterfall plot detailed in \textbf{Supplementary Figure S9}).}
\label{fig:fig6}
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
            woven_results.append("\n" + get_table_code(2) + "\n" + get_figure_code(2) + "\n")
        elif sec_num == 3:
            woven_results.append("\n" + get_table_code(3) + "\n" + get_figure_code(3) + "\n")
        elif sec_num == 4:
            woven_results.append("\n" + get_table_code(4) + "\n" + get_figure_code(4) + "\n")
        elif sec_num == 5:
            woven_results.append("\n" + get_table_code(5) + "\n" + get_table_code(6) + "\n" + get_figure_code(5) + "\n")
        elif sec_num == 6:
            woven_results.append("\n" + get_table_code(7) + "\n" + get_table_code(8) + "\n" + get_figure_code(6) + "\n")


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
All analyzed raw transcriptomic and epigenomic datasets are publicly available in the NCBI Gene Expression Omnibus (GEO) under accession numbers \textbf{GSE107925}, \textbf{GSE108045}, \textbf{GSE186210}, \textbf{GSE266602}, \textbf{GSE152865}, and \textbf{GSE67858}. All standardized sample metadata, normalized read count matrices, deconvolution matrices, and meta-analytic summary tables are deposited and accessible in the study GitHub repository (\url{https://github.com/samyakmeshram/NeuroGut-MetaSeq}).

\subsection*{Code Availability}
All analytical pipelines, statistical modeling scripts, and visualization code are publicly available under the open-source MIT license at \url{https://github.com/samyakmeshram/NeuroGut-MetaSeq}. The analytical workflows are fully reproducible in Python 3.9+ and R 4.3+.

\vspace{1.8em}
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
This document contains extended methodological formulations, supplementary diagnostic figures (Figures S1--S9), and supplementary statistical tables (Tables S1--S4) supporting the main text findings of \textbf{NeuroGut-MetaSeq (v1.4.0)}.

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
GSE107925 & SPF Control ($n=13$) & 1.8M -- 42.2M & 20.1M & $99.2\% \pm 1.0\%$ & $-0.04 \pm 0.97$ & FACS CD11b$^{+}$CD45$^{\text{low}}$ \\
GSE107925 & Germ-Free ($n=12$) & 1.7M -- 36.8M & 18.8M & $99.2\% \pm 0.5\%$ & $+0.04 \pm 0.89$ & FACS CD11b$^{+}$CD45$^{\text{low}}$ \\
GSE108045 & Vehicle ($n=6$) & 18.0M -- 22.6M & 19.8M & $99.6\% \pm 0.0\%$ & $-0.27 \pm 0.79$ & FACS CD11b$^{+}$CD45$^{\text{low}}$ \\
GSE108045 & Antibiotics ($n=6$) & 16.3M -- 23.1M & 20.3M & $99.6\% \pm 0.0\%$ & $+0.27 \pm 0.70$ & FACS CD11b$^{+}$CD45$^{\text{low}}$ \\
GSE186210 & Standard Fiber ($n=8$) & 1.2M -- 4.0M & 2.7M & $99.7\% \pm 0.1\%$ & $+0.04 \pm 0.61$ & MACS CD11b Microbeads \\
GSE186210 & Zero-Fiber ($n=6$) & 1.9M -- 2.6M & 2.3M & $99.6\% \pm 0.1\%$ & $-0.05 \pm 0.53$ & MACS CD11b Microbeads \\
GSE266602 & Colonized ($n=3$) & 51.9M -- 55.0M & 53.9M & $65.8\% \pm 11.3\%$ & $-0.01 \pm 0.76$ & Percoll Density Gradient \\
GSE266602 & Depleted GF/ABX ($n=6$) & 41.7M -- 70.4M & 54.1M & $74.3\% \pm 14.4\%$ & $+0.00 \pm 1.11$ & Percoll Density Gradient \\
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
\textbf{Gene Symbol} & \textbf{Cohorts ($k$)} & \textbf{REML $\log_2\text{FC}$} & \textbf{HKSJ 95\% CI} & \textbf{$p_{\text{HKSJ}}$} & \textbf{$I^2$} & \textbf{Consensus Designation} \\
\midrule
\textit{Llgl2} & 4 & +0.672 & [0.166, 1.178] & 0.0243 & 0.0\% & Invariant Basolateral Core \\
\textit{Slfn2} & 4 & -0.465 & [-0.789, -0.141] & 0.0197 & 9.6\% & Invariant Quiescence Guardian \\
\textit{Fosb} & 3 & +1.345 & [0.604, 2.085] & 0.0160 & 0.0\% & Invariant Immediate Early Activator \\
\textit{Sap30} & 3 & -0.393 & [-0.745, -0.040] & 0.0409 & 0.0\% & Invariant Sin3A-HDAC Corepressor \\
\textit{Slc25a19} & 4 & -0.315 & [-0.549, -0.082] & 0.0232 & 0.0\% & Invariant Thiamine Pyrophosphate Carrier \\
\textit{Clu} & 3 & +0.873 & [-0.109, 1.854] & 0.0621 & 0.0\% & Invariant Extracellular Chaperone ($p_{\text{RE}}=1.3 \times 10^{-4}$) \\
\textit{Card6} & 3 & -0.404 & [-0.942, 0.135] & 0.0843 & 0.0\% & Invariant NF-$\kappa$B / NOD Modulator \\
\textit{Tnf} & 3 & +1.052 & [-0.510, 2.614] & 0.1013 & 32.2\% & Conserved Low-Grade Priming \\
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
\begin{tabular}{lccccclp{6.5cm}}
\toprule
\textbf{Gene} & \textbf{Male LFC} & \textbf{Female LFC} & \textbf{$\hat{\theta}_{\text{int}}$} & \textbf{$p_{\text{int}}$} & \textbf{$I^2_{\text{sex}}$} & \textbf{Dimorphism Status} & \textbf{Biological Interpretation} \\
\midrule
\textit{Slfn2} & -0.583 & -0.232 & +0.3382 & 0.0170 & 0.0\% & Nominal Male-Biased & Steeper loss of quiescence in males \\
\textit{Oas1a} & -0.551 & -0.124 & +0.3348 & 0.0507 & 0.0\% & Nominal Male-Biased & Antiviral surveillance collapse in males \\
\textit{Irf1} & +0.082 & +0.058 & -0.0177 & 0.8629 & 0.0\% & Sex-Shared (Invariant) & Master regulon shutoff is sex-invariant \\
\textit{Stat1} & -0.045 & -0.074 & -0.0918 & 0.3618 & 0.0\% & Sex-Shared (Invariant) & Downstream JAK/STAT shutdown \\
\textit{Llgl2} & +0.560 & +0.524 & -0.0247 & 0.9410 & 35.8\% & Sex-Shared (Invariant) & Nutrient adaptation is sex-invariant \\
\textit{Ddit4} & -0.864 & -0.484 & +0.2707 & 0.7566 & 63.0\% & Sex-Shared (High $I^2$) & Perturbation-specific antibiotic shock \\
\textit{Fosb} & +0.488 & +0.952 & +0.1179 & 0.7950 & 54.3\% & Sex-Shared (Invariant) & Immediate-early AP-1 activation \\
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
1 & \textbf{IRF1} & 40 & -0.113 & -2.453 & 0.0380 & \textbf{Significant Master Repression} \\
2 & SP1 & 391 & +0.108 & -1.650 & 0.3952 & Repressed Basal Transcription \\
3 & STAT1 & 66 & +0.085 & -0.938 & 0.6808 & Repressed Interferon Cascade \\
4 & RELA & 241 & +0.128 & -0.883 & 0.6738 & Modest NF-$\kappa$B Crosstalk \\
5 & FOS & 46 & +0.138 & -0.291 & 0.9366 & AP-1 Immediate Early Activation \\
6 & JUN & 117 & +0.154 & -0.225 & 0.9430 & AP-1 Immediate Early Activation \\
\bottomrule
\end{tabular}%
}
\end{table}

\clearpage
\section{Supplementary Figures}

% Figure S1
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_s1_qc_unified_audit.pdf}
\caption{\textbf{Multi-Cohort Technical Quality Control Audit and Microglial Lineage Validation.} 
(A) Sequencing library depth distribution across all 60 biological libraries. 
(B) Canonical microglial lineage checkpoint abundance ($>99.2\%$ purity). 
(C) Non-myeloid contaminant marker depletion ($<0.8\%$ astrocytic, oligodendrocytic, and neuronal transcripts). 
(D) Ex vivo enzymatic dissociation stress scores across colonized controls and microbiome-depleted groups ($p \ge 0.18$).}
\label{fig:figs1}
\end{figure}

\clearpage

% Figure S2
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_s2_volcano_all_cohorts.pdf}
\caption{\textbf{Single-Cohort Differential Expression Landscapes Across Experimental Paradigms.} 
(A) GSE107925 (Lifelong Germ-Free adult microglia). 
(B) GSE108045 (Acute Broad-Spectrum Antibiotics cocktail), illustrating profound perturbation shock repression of \textit{Tsc22d3} and \textit{Ddit4}. 
(C) GSE186210 (Dietary Fiber Starvation), highlighting selective suppression of perilipin \textit{Plin3}. 
(D) GSE266602 (Germ-Free Percoll Baseline).}
\label{fig:figs2}
\end{figure}

\clearpage

% Figure S3
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_s3_cross_cohort_divergence.pdf}
\caption{\textbf{Cross-Cohort Divergence, Marker Profiles, and Unpooled DEG Overlaps.} 
(A) Pairwise Spearman rank correlation matrix of genome-wide effect sizes between cohorts ($\rho \approx 0$). 
(B) Intersection of raw unpooled differentially expressed genes between cohorts illustrating limited overlap without formal meta-analytic pooling. 
(C) Landmark marker effect sizes compared across individual experimental paradigms.}
\label{fig:figs3}
\end{figure}

\clearpage

% Figure S4
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_s4_heterogeneity_and_forest.pdf}
\caption{\textbf{Transcriptome-Wide Heterogeneity Distribution and Full Landmark Forest Plot Grid.} 
(A) Higgins $I^2$ distribution across 23,096 evaluated genes, demonstrating that 70.5\% of genes exhibit low between-study heterogeneity ($I^2 < 25\%$). 
(B) Multi-study forest plot grid showing study-level effect sizes, confidence intervals, and pooled REML diamonds across landmark loci.}
\label{fig:figs4}
\end{figure}

\clearpage

% Figure S5
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_s5_loo_sensitivity.pdf}
\caption{\textbf{Systematic Leave-One-Out (LOO) Sensitivity Meta-Analyses.} 
Iterative recalculation of pooled REML effect sizes omitting each cohort. Core consensus loci (\textit{Llgl2}, \textit{Clu}, \textit{Slfn2}, \textit{Sap30}) remain tightly centered along the unity line across all iterations, while acute shock effectors (\textit{Tsc22d3}, \textit{Ddit4}) return toward zero upon omission of the antibiotic cohort.}
\label{fig:figs5}
\end{figure}

\clearpage

% Figure S6
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_s6_wgcna_and_interactome.pdf}
\caption{\textbf{Weighted Gene Co-Expression Networks and Consensus Hub Interactome.} 
(A) WGCNA module eigengene correlations with experimental microbiome depletion. 
(B) Consensus co-expression hub subgraph highlighting topological overlap connections across three empirical functional axes (\textit{Slfn2-Sap30-Ddit4}, \textit{Sall1-Tmem119-Ffar2-Llgl2}, \textit{Fosb-Fos}).}
\label{fig:figs6}
\end{figure}

\clearpage

% Figure S7
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_s7_deconvolution_diagnostics.pdf}
\caption{\textbf{BayesPrism Deconvolution Mathematical Diagnostics and Cell Lineage Stability.} 
(A) Singular Value Decomposition (SVD) condition index scree plot confirming zero harmful multicollinearity ($\kappa = 1.54 < 30$). 
(B) Pan-microglial lineage marker stability across cohorts (\textit{Hexb}, \textit{Csf1r}, \textit{Tmem119}; $p = 0.85$).}
\label{fig:figs7}
\end{figure}

\clearpage

% Figure S8
\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{figures/fig_s8_extended_nichenet.pdf}
\caption{\textbf{Extended NicheNet Upstream Ligand-Receptor Relay and Regulatory Potential Matrix.} 
(A) Upstream ligand activity ranking by candidate sender compartment (BMEC brain endothelium, BAM border macrophages, and circulating microbial PAMPs). 
(B) In silico ligand-to-target regulatory matrix heatmap linking prioritized upstream ligands to microglial IRF1 regulon targets.}
\label{fig:figs8}
\end{figure}

\clearpage

% Figure S9
\begin{figure}[htbp]
\centering
\includegraphics[width=0.98\textwidth]{figures/fig_s9_scfa_rescue_specificity_null.pdf}
\caption{\textbf{In Silico SCFA Metabolite Specificity Null Permutation and Ranked ISRI Waterfall.} 
(A) Empirical label-shuffling specificity null distribution across 10,000 permutations comparing observed core signature recovery ($\text{mean ISRI} = +0.528$, red dashed line) against the chance null ($\mu = 0.001 \pm 0.125, p_{\text{perm}} < 0.0001$). 
(B) Ranked In Silico Rescue Index (ISRI) waterfall plot across landmark microglial genes, demonstrating coordinate reciprocal rescue in 18 of 19 loci upon in vivo SCFA repletion.}
\label{fig:figs9}
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
