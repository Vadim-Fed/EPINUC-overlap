Multiplexed Single-Molecule Epigenetic Analysis of Plasma-Isolated Nucleosomes for Cancer Diagnostics
The analysis of cell-free DNA (cfDNA) in plasma represents a rapidly advancing field in medicine, providing information on pathological processes in the body. 
Blood cfDNA is in the form of nucleosomes, which maintain their tissue- and cancer-specific epigenetic modifications. 
We developed EPINuc, a single-molecule multi-parametric-based assay to comprehensively profile the Epigenetics of Plasma Isolated Nucleosomes, 
DNA methylation and cancer-specific protein biomarkers. Our system allows high resolution detection of six active and repressive histone modifications, 
their ratios and combinatorial patterns, on millions of individual nucleosomes by single-molecule imaging. 
In addition, it provides sensitive and quantitative data on plasma proteins, including detection of non-secreted tumor-specific proteins such as mutant p53. 
Applying this analysis to a cohort of plasma samples detected colorectal cancer at high accuracy and sensitivity, even at early stages. 
Finally, combining EPINuc with single-molecule DNA sequencing revealed the tissue-of-origin of the tumor. 
EPINuc generates multi-layered information encompassing various cellular pathways from very limited liquid biopsy material, establishing a novel approach for cancer diagnostics.

# EPINuc-overlap

This is the script used to assess overlap significance between plasma sequenced reads positive for 5hmC/antibody of interest with a given tissue chip-seq Peaks. 
The script performs bootstrapping simulations to analyze the significance of the single-molecule reads overlap with various tissue peaks.

Analysis was ran on R vs 4.0

The files needed to run the analysis are:
1. Chromosome_length.csv - supplied here at the main branch.
2. BED file of plasma sequenced reads, positive for 5hmC/antibody of interest - available upon request.
3. BED files for the tested tissues - Available at https://www.encodeproject.org/ with accession numbers for relevant chip-seq datasets listed at Table S4.
