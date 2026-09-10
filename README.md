# Single-cell transcriptomic atlas of human inferior rectus muscle in inactive thyroid-associated orbitopathy

Analysis code for the study: "Single-cell analysis reveals that endothelial–fibroblast crosstalk and midkine signaling drive tissue remodeling in inactive thyroid-associated orbitopathy".

## Contents

- `TAO_analysis_final_script.R` — self-contained analysis pipeline (R / Seurat)

The script reproduces all analyses reported in the manuscript:

1. 10X data loading, quality control, and cell filtering
2. Sample merging, normalization, and Harmony batch correction
3. Clustering and cell-type annotation (11 major cell types)
4. Cell-proportion statistics with Fisher's exact tests (FDR-adjusted)
5. Subclustering of 10 compartments (fibroblast, B cell, endothelial, myeloid, NK/T, pericyte, myocyte, muscle satellite, smooth muscle)
6. Differential expression (TED vs. control) and GO / KEGG / GSEA enrichment per subcluster
7. Fibroblast pseudotime trajectory (Monocle2, Fib1 root)
8. Endothelial EndoMT pseudotime trajectory (Monocle2, Endo1 root)
9. CellChat cell–cell communication analysis with MK pathway focus

## Requirements

R ≥ 4.1 with the following packages: Seurat (≥ 5.0, required for `JoinLayers()`
and the v5 assay/layer API), harmony, tidyverse, plyr, scales,
RColorBrewer, patchwork, CellChat, NMF, ggalluvial, pheatmap, ComplexHeatmap,
cowplot, gridExtra, clusterProfiler, enrichplot, org.Hs.eg.db, monocle (Monocle2),
Biobase, viridis, data.table.

## Usage

1. Download the raw 10X data from GEO (accession **GSE346666**).
2. Edit `base_dir` in Section 0 of the script to point to your local directory,
   which must contain:
   - `表达矩阵文件/总/` — subfolders of 10X output (one per sequencing library)
   - `表达矩阵文件/clinic.csv` — donor clinical metadata
3. Run the script from start to end.

All outputs (tables, figures, and `.RData` objects) are written to the
`results` directory created under `base_dir`.

## Data availability

Raw and processed scRNA-seq data have been deposited in the NCBI Gene
Expression Omnibus (GEO) under accession number **GSE346666** (currently under
controlled access; the records will be released for public access upon journal
publication, in accordance with GEO policy).

## Citation

If you use this code, please cite the associated article.
