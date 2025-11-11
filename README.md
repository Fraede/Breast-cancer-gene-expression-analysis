Computational Analysis of Gene Expression Progression in Early-Stage Breast Cancer (GSE16873)

Author: Fareed Dubiure

Language: Python

Tools: pandas, NumPy, SciPy, matplotlib, statsmodels

Project Overview: This project investigates gene expression changes across progressive breast tissue stages—normal, simple ductal hyperplasia (SDH), atypical ductal hyperplasia (ADH), and ductal carcinoma in situ (DCIS)
using the publicly available microarray dataset GSE16873 from the NCBI Gene Expression Omnibus.

Objective: To identify genes that show consistent expression changes during early breast cancer development and highlight potential molecular markers of early transformation.

Methods:
Processed and normalized microarray data using Python.
Performed differential expression analysis with t-tests and log₂ fold changes.
Applied FDR correction for multiple testing.
Visualized results using volcano plots and expression trend analysis across four tissue stages.

Key Findings:
The majority of significantly altered genes were downregulated in DCIS compared to normal tissue.
ACTA2 showed a consistent negative progression, suggesting early loss of myoepithelial integrity and cytoskeletal organization.
Results indicate early tumorigenesis may involve transcriptional suppression rather than activation of oncogenic programs.

Figures:
Volcano plot of differentially expressed genes.
Trend visualization of the top 10 significant genes.

Dataset Source: GSE16873 – NCBI GEO

Platform: Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)

Citation
If you use this repository or dataset, please cite:

Allred, D.C. et al. (2008). Ductal carcinoma in situ and the emergence of diversity during breast cancer evolution. Clinical Cancer Research, 7(12), 3703–3712.

NCBI GEO Accession: GSE16873.
