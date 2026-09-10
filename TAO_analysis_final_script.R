################################################################################
# TAO Single-Cell RNA-seq Analysis - Final Script for GEO Submission
# Thyroid Associated Orbitopathy (TAO) scRNA-seq comprehensive analysis pipeline
# Author: Jinwei Cheng, MD (on behalf of all authors)
# Date: 2026-06-22
# Description: This script reproduces all analyses from the TAO scRNA-seq study,
#              including QC, Harmony integration, cell-type annotation, subclustering,
#              differential expression & enrichment (GO/KEGG/GSEA), fibroblast and
#              endothelial (EndoMT) pseudotime trajectories, and CellChat analysis.
#              It is self-contained and can be run from start to end,
#              assuming raw 10X data is available at the specified directory.
################################################################################

# ============================================================================
# SECTION 0: Global Configuration / 全局配置
# ============================================================================
# Modify base_dir to point to your data directory
# base_dir should contain a subfolder "总/" with 10X sample folders
# clinic.csv should be at base_dir/clinic.csv

base_dir <- "F:/OneDrive/桌面/单细胞"
data_root <- file.path(base_dir, "表达矩阵文件/总")
clinic_file <- file.path(base_dir, "表达矩阵文件/clinic.csv")
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# SECTION 1: Package Loading / 包加载
# ============================================================================
library(Seurat)
library(tidyverse)
library(dplyr)
library(harmony)
library(patchwork)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(scales)
library(stringr)

# ============================================================================
# SECTION 2: Data Loading & QC / 数据加载与质控
# ============================================================================

# 2.1 Read 10X data for all samples / 读取10X数据
dir_name <- list.dirs(data_root, full.names = FALSE, recursive = FALSE)
dir_name <- dir_name[dir_name != ""]
cat("Detected samples:", paste(dir_name, collapse = ", "), "\n")

scrnaseq_list <- lapply(dir_name, function(sample) {
  data_path <- file.path(data_root, sample)
  counts <- Read10X(data.dir = data_path)
  CreateSeuratObject(
    counts = counts,
    project = sample,
    min.cells = 3,
    min.features = 300
  )
})
names(scrnaseq_list) <- dir_name

# 2.2 Calculate QC metrics / 计算质控指标
calc_qc_metrics <- function(sc_obj) {
  sc_obj[["mt_percent"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
  hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2",
                 "HBM", "HBQ1", "HBZ")
  hb_genes <- intersect(hb_genes, rownames(sc_obj))
  if (length(hb_genes) > 0) {
    sc_obj[["hb_percent"]] <- PercentageFeatureSet(sc_obj, features = hb_genes)
  } else {
    sc_obj[["hb_percent"]] <- 0
  }
  return(sc_obj)
}
scrnaseq_list <- lapply(scrnaseq_list, calc_qc_metrics)

# 2.3 Filter cells / 过滤细胞
scrnaseq_list <- lapply(scrnaseq_list, function(sc) {
  subset(sc,
    subset = nFeature_RNA > 300 & nFeature_RNA < 5000 &
            mt_percent < 25 &
            hb_percent < 3 &
            nCount_RNA > 1000 &
            nCount_RNA < quantile(nCount_RNA, 0.97)
  )
})
cat("Cells after filtering:\n")
print(sapply(scrnaseq_list, ncol))

# ============================================================================
# SECTION 3: Merge, Normalize, & Harmony / 合并、归一化与Harmony批次校正
# ============================================================================

# 3.1 Merge samples / 合并样本
scRNA_merged <- merge(
  x = scrnaseq_list[[1]],
  y = scrnaseq_list[-1],
  add.cell.ids = names(scrnaseq_list)
)
cat("Total cells after merge:", ncol(scRNA_merged), "\n")

# 3.2 Normalize & feature selection / 归一化与特征选择
scRNA_merged <- scRNA_merged %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("mt_percent", "nCount_RNA")) %>%
  RunPCA(npcs = 50, verbose = FALSE)

# 3.3 Harmony batch correction / Harmony批次校正
scRNA_harmony <- RunHarmony(scRNA_merged, group.by.vars = "orig.ident")

# ============================================================================
# SECTION 4: Clustering & Annotation / 聚类与注释
# ============================================================================

# 4.1 Clustering (dims=1:40, resolution=0.2) / 聚类
scRNA_harmony <- scRNA_harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:40) %>%
  FindClusters(resolution = 0.2)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:40)

# JoinLayers for downstream analysis
scRNA_harmony <- JoinLayers(scRNA_harmony)

# 4.2 Marker gene identification / Marker基因鉴定
markers <- FindAllMarkers(
  scRNA_harmony,
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.5
)
sig.markers <- markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))
top10_markers <- sig.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
write.csv(top10_markers, file.path(results_dir, "top10_markers_per_cluster.csv"),
          row.names = FALSE)

# 4.3 Cell type annotation / 细胞类型注释
celltype_anno <- data.frame(
  ClusterID = 0:max(as.integer(scRNA_harmony$seurat_clusters)),
  celltype = "unknown",
  stringsAsFactors = FALSE
)
celltype_anno$celltype[celltype_anno$ClusterID %in% c(0, 8, 9)] <- "Fibroblast"
celltype_anno$celltype[celltype_anno$ClusterID == 1]  <- "Endothelial cell"
celltype_anno$celltype[celltype_anno$ClusterID == 3]  <- "Myeloid"
celltype_anno$celltype[celltype_anno$ClusterID == 2]  <- "NK/T cell"
celltype_anno$celltype[celltype_anno$ClusterID %in% c(10, 12)] <- "B cell"
celltype_anno$celltype[celltype_anno$ClusterID == 4]  <- "Muscle satellite cell"
celltype_anno$celltype[celltype_anno$ClusterID == 5]  <- "Smooth Muscle Cell"
celltype_anno$celltype[celltype_anno$ClusterID == 7]  <- "Pericyte"
celltype_anno$celltype[celltype_anno$ClusterID == 11] <- "neuro"
celltype_anno$celltype[celltype_anno$ClusterID == 6]  <- "Myocyte"
celltype_anno$celltype[celltype_anno$ClusterID == 13] <- "epi"

scRNA_harmony@meta.data$celltype <- factor(
  plyr::mapvalues(
    scRNA_harmony$seurat_clusters,
    from = celltype_anno$ClusterID,
    to = celltype_anno$celltype
  )
)

# ============================================================================
# SECTION 5: Clinical Data Integration / 临床数据整合
# ============================================================================
patients_metadata <- data.table::fread(clinic_file, header = TRUE)
metadata <- FetchData(scRNA_harmony, "orig.ident") %>%
  mutate(cell_id = rownames(.)) %>%
  left_join(patients_metadata, by = "orig.ident")
rownames(metadata) <- metadata$cell_id
scRNA_harmony <- AddMetaData(scRNA_harmony, metadata = metadata)

# ============================================================================
# SECTION 6: Cell Proportion Statistics & Fisher Test / 细胞比例统计与Fisher检验
# ============================================================================
Idents(scRNA_harmony) <- "celltype"

group_celltype_counts <- table(
  CellType = Idents(scRNA_harmony),
  Group = scRNA_harmony$group
)

celltype_proportions <- as.data.frame(prop.table(group_celltype_counts, margin = 2))
colnames(celltype_proportions) <- c("celltype", "group", "proportion")
celltype_proportions$proportion <- round(celltype_proportions$proportion * 100, 2)

# Fisher exact test for each cell type / Fisher精确检验
counts_df <- as.data.frame.matrix(group_celltype_counts)
if ("normal" %in% colnames(counts_df) && "TED" %in% colnames(counts_df)) {
  fisher_results <- data.frame(
    celltype = rownames(counts_df),
    normal_count = counts_df$normal,
    ted_count = counts_df$TED,
    p_value = NA,
    odds_ratio = NA,
    stringsAsFactors = FALSE
  )
  for (i in 1:nrow(counts_df)) {
    celltype <- rownames(counts_df)[i]
    other_count <- colSums(counts_df[-i, c("normal", "TED")])
    contingency_table <- matrix(
      c(counts_df[i, "normal"], counts_df[i, "TED"],
        other_count["normal"], other_count["TED"]),
      nrow = 2,
      dimnames = list(c(celltype, "Others"), c("normal", "TED"))
    )
    fisher_test <- fisher.test(contingency_table)
    fisher_results$p_value[i] <- fisher_test$p.value
    fisher_results$odds_ratio[i] <- fisher_test$estimate
  }
  fisher_results$adj_p_value <- p.adjust(fisher_results$p_value, method = "fdr")
  write.csv(fisher_results,
            file.path(results_dir, "celltype_fisher_test_results.csv"),
            row.names = FALSE)
}

# Save main object / 保存主对象
save(scRNA_harmony, celltype_proportions,
     file = file.path(results_dir, "scRNA_analysis_with_proportions_final.RData"))

# ============================================================================
# SECTION 7: Subcluster Analyses / 亚群再聚类分析
# ============================================================================
# Generic subcluster pipeline function / 通用亚群再聚类流程函数
run_subcluster <- function(parent_obj, celltype_name, obj_name,
                           nfeatures = 2000, dims_range, resolution,
                           annotations, scale_all_genes = TRUE) {
  # Subset target cell type / 提取目标细胞类型
  target_cells <- rownames(subset(parent_obj@meta.data, celltype == celltype_name))
  sub_obj <- subset(parent_obj, cells = target_cells)

  # Preprocessing / 数据预处理
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst",
                                   nfeatures = nfeatures)
  if (scale_all_genes) {
    sub_obj <- ScaleData(sub_obj, features = rownames(sub_obj))
  } else {
    sub_obj <- ScaleData(sub_obj)
  }
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(sub_obj), verbose = FALSE)

  # Harmony / 批次校正
  sub_obj <- RunHarmony(sub_obj, group.by.vars = "orig.ident")

  # Clustering / 聚类
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = dims_range)
  sub_obj <- FindClusters(sub_obj, resolution = resolution)
  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = dims_range)

  # FindAllMarkers / Marker基因
  sub_markers <- FindAllMarkers(
    sub_obj,
    test.use = "wilcox",
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.25
  )
  write.csv(sub_markers,
            file.path(results_dir, paste0(obj_name, "_markers.csv")),
            row.names = FALSE)

  # Annotation / 注释
  df_celltype_anno <- data.frame(
    ClusterID = 0:max(as.integer(sub_obj$seurat_clusters)),
    celltype = "unknown",
    stringsAsFactors = FALSE
  )
  for (i in seq_along(annotations)) {
    df_celltype_anno$celltype[df_celltype_anno$ClusterID == annotations[[i]]$cluster] <-
      annotations[[i]]$name
  }
  sub_obj@meta.data$celltype <- factor(
    plyr::mapvalues(
      sub_obj$seurat_clusters,
      from = df_celltype_anno$ClusterID,
      to = df_celltype_anno$celltype
    )
  )

  # Proportion stats & Fisher test / 比例统计与Fisher检验
  Idents(sub_obj) <- "celltype"
  sub_group_counts <- table(CellType = Idents(sub_obj), Group = sub_obj$group)
  sub_counts_df <- as.data.frame.matrix(sub_group_counts)

  sub_fisher_results <- data.frame(
    celltype = rownames(sub_counts_df),
    stringsAsFactors = FALSE
  )
  if ("normal" %in% colnames(sub_counts_df) && "TED" %in% colnames(sub_counts_df)) {
    sub_fisher_results$normal_count <- sub_counts_df$normal
    sub_fisher_results$ted_count <- sub_counts_df$TED
    sub_fisher_results$p_value <- NA
    sub_fisher_results$odds_ratio <- NA
    for (i in 1:nrow(sub_counts_df)) {
      ct <- rownames(sub_counts_df)[i]
      other <- colSums(sub_counts_df[-i, c("normal", "TED")])
      cont_table <- matrix(
        c(sub_counts_df[i, "normal"], sub_counts_df[i, "TED"],
          other["normal"], other["TED"]),
        nrow = 2,
        dimnames = list(c(ct, "Others"), c("normal", "TED"))
      )
      ft <- fisher.test(cont_table)
      sub_fisher_results$p_value[i] <- ft$p.value
      sub_fisher_results$odds_ratio[i] <- ft$estimate
    }
    sub_fisher_results$adj_p_value <- p.adjust(sub_fisher_results$p_value, method = "fdr")
    write.csv(sub_fisher_results,
              file.path(results_dir, paste0(obj_name, "_fisher_test.csv")),
              row.names = FALSE)
  }

  # Save / 保存
  save(list = c(obj_name, "df_celltype_anno", "sub_fisher_results"),
       file = file.path(results_dir, paste0(obj_name, "_analysis_final.RData")))
  assign(obj_name, sub_obj, envir = .GlobalEnv)
  return(sub_obj)
}

# ------------------------------------------------------------------
# 7.1 Fibroblast subcluster / 成纤维细胞再聚类
# ------------------------------------------------------------------
Fibsub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "Fibroblast",
  obj_name = "Fibsub",
  nfeatures = 2000,
  dims_range = 1:40,
  resolution = 0.2,
  annotations = list(
    list(cluster = 0, name = "Fib1"),
    list(cluster = 1, name = "Fib2"),
    list(cluster = 2, name = "Fib3"),
    list(cluster = 3, name = "Fib4"),
    list(cluster = 4, name = "Fib5")
  )
)

# ------------------------------------------------------------------
# 7.2 B cell subcluster / B细胞再聚类
# ------------------------------------------------------------------
B_sub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "B cell",
  obj_name = "B_sub",
  nfeatures = 2000,
  dims_range = 1:40,
  resolution = 0.05,
  annotations = list(
    list(cluster = 0, name = "B cell 1"),
    list(cluster = 1, name = "B cell 2")
  )
)

# ------------------------------------------------------------------
# 7.3 Endothelial cell subcluster / 内皮细胞再聚类
# ------------------------------------------------------------------
Endothelialsub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "Endothelial cell",
  obj_name = "Endothelialsub",
  nfeatures = 2000,
  dims_range = 1:20,
  resolution = 0.1,
  annotations = list(
    list(cluster = 0, name = "Endo1"),
    list(cluster = 1, name = "Endo2"),
    list(cluster = 2, name = "Endo3"),
    list(cluster = 3, name = "Endo4")
  )
)

# ------------------------------------------------------------------
# 7.4 Myeloid subcluster / 髓系细胞再聚类
# ------------------------------------------------------------------
Myeloidsub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "Myeloid",
  obj_name = "Myeloidsub",
  nfeatures = 2000,
  dims_range = 1:40,
  resolution = 0.05,
  annotations = list(
    list(cluster = 0, name = "Mac"),
    list(cluster = 1, name = "Mono"),
    list(cluster = 2, name = "DC")
  )
)

# ------------------------------------------------------------------
# 7.5 NK/T cell subcluster / NK/T细胞再聚类
# ------------------------------------------------------------------
NK_Tsub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "NK/T cell",
  obj_name = "NK_Tsub",
  nfeatures = 2000,
  dims_range = 1:35,
  resolution = 0.05,
  annotations = list(
    list(cluster = 0, name = "T cell"),
    list(cluster = 1, name = "NK cell"),
    list(cluster = 2, name = "T cell")
  )
)

# ------------------------------------------------------------------
# 7.6 Pericyte subcluster / 周细胞再聚类
# Bug fix: Pericyte_subb → Pericyte_sub (original had typo with extra 'b')
# ------------------------------------------------------------------
Pericyte_sub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "Pericyte",
  obj_name = "Pericyte_sub",
  nfeatures = 2000,
  dims_range = 1:30,
  resolution = 0.05,
  annotations = list(
    list(cluster = 0, name = "pericyte1"),
    list(cluster = 1, name = "pericyte2"),
    list(cluster = 2, name = "pericyte3")
  )
)

# ------------------------------------------------------------------
# 7.7 Myocyte subcluster / 肌细胞再聚类
# ------------------------------------------------------------------
Myo_sub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "Myocyte",
  obj_name = "Myo_sub",
  nfeatures = 2000,
  dims_range = 1:25,
  resolution = 0.05,
  annotations = list(
    list(cluster = 0, name = "myo1"),
    list(cluster = 1, name = "myo2"),
    list(cluster = 2, name = "myo3"),
    list(cluster = 3, name = "myo4")
  )
)

# ------------------------------------------------------------------
# 7.8 Muscle satellite cell subcluster / 肌卫星细胞再聚类
# ------------------------------------------------------------------
Satellite_sub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "Muscle satellite cell",
  obj_name = "Satellite_sub",
  nfeatures = 2000,
  dims_range = 1:20,
  resolution = 0.1,
  annotations = list(
    list(cluster = 0, name = "satellite1"),
    list(cluster = 1, name = "satellite2"),
    list(cluster = 2, name = "satellite3")
  )
)
# Note: Satellite subcluster uses logfc.threshold=0.2, min.pct=0.2 for FindAllMarkers
# The run_subcluster function uses 0.25/0.25 by default; adjust if needed.

# ------------------------------------------------------------------
# 7.9 Smooth Muscle Cell subcluster / 平滑肌细胞再聚类
# ------------------------------------------------------------------
SMC_sub <- run_subcluster(
  parent_obj = scRNA_harmony,
  celltype_name = "Smooth Muscle Cell",
  obj_name = "SMC_sub",
  nfeatures = 2000,
  dims_range = 1:40,
  resolution = 0.05,
  annotations = list(
    list(cluster = 0, name = "smc1"),
    list(cluster = 1, name = "smc2"),
    list(cluster = 2, name = "smc3")
  )
)

# ------------------------------------------------------------------
# 7.10 Subpopulation Integration Back to Main Object / 亚群整合回大群
# ------------------------------------------------------------------
# Write subcluster labels (Fib/Endo/Myeloid) back to scRNA_harmony
# This enables downstream CellChat analysis with subcluster-level resolution

# Verify matching between subcluster objects and main object
fib_cell_ids <- WhichCells(scRNA_harmony, idents = "Fibroblast")
endo_cell_ids <- WhichCells(scRNA_harmony, idents = "Endothelial cell")
myeloid_cell_ids <- WhichCells(scRNA_harmony, idents = "Myeloid")

cat("Fibroblast subcluster matching rate:",
    mean(colnames(Fibsub) %in% fib_cell_ids) * 100, "%\n")
cat("Endothelial subcluster matching rate:",
    mean(colnames(Endothelialsub) %in% endo_cell_ids) * 100, "%\n")
cat("Myeloid subcluster matching rate:",
    mean(colnames(Myeloidsub) %in% myeloid_cell_ids) * 100, "%\n")

# Update main object identities with subcluster labels
Idents(scRNA_harmony, cells = colnames(Fibsub)) <- Idents(Fibsub)
Idents(scRNA_harmony, cells = colnames(Endothelialsub)) <- Idents(Endothelialsub)
Idents(scRNA_harmony, cells = colnames(Myeloidsub)) <- Idents(Myeloidsub)

# Write integrated labels to metadata for persistent tracking
scRNA_harmony@meta.data$integrated_celltype <- as.character(Idents(scRNA_harmony))

# Visualize integrated result
pdf(file.path(results_dir, "scRNA_harmony_subclusters_integrated.pdf"),
    width = 10, height = 8)
DimPlot(scRNA_harmony, label = TRUE, repel = TRUE) +
  ggtitle("Subcluster Labels Integrated into Main Object")
dev.off()

cat("Subpopulation integration complete.\n")
cat("integrated_celltype column added to scRNA_harmony metadata.\n")

# ============================================================================
# SECTION 8: GSEA & Enrichment Analysis for Each Subcluster / GSEA与富集分析
# ============================================================================
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(Cairo)

# Generic enrichment pipeline function / 通用富集分析流程函数
run_enrichment <- function(sub_obj, obj_name, ident_col = "celltype") {
  # Differential analysis: TED vs normal / 差异分析
  Idents(sub_obj) <- ident_col

  group_degs <- FindMarkers(
    sub_obj,
    group.by = "group",
    ident.1 = "TED",
    ident.2 = "normal",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = FALSE,
    test.use = "wilcox",
    verbose = FALSE
  ) %>% mutate(gene = rownames(.))

  # Significant DEGs / 显著差异基因
  group_degs_sig <- group_degs %>%
    filter(
      !is.na(p_val_adj),
      p_val_adj < 0.05,
      abs(avg_log2FC) > 0.5,
      pct.1 > 0.1,
      pct.2 > 0.1
    ) %>%
    arrange(desc(abs(avg_log2FC)))

  write.csv(group_degs_sig,
            file.path(results_dir, paste0(obj_name, "_group_DEGs.csv")),
            row.names = FALSE)

  # Celltype markers & top100 / 细胞类型Marker及top100
  celltype_markers <- FindAllMarkers(
    sub_obj,
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.25,
    verbose = FALSE
  )

  top100_markers <- celltype_markers %>%
    group_by(cluster) %>%
    top_n(n = 100, wt = avg_log2FC) %>%
    ungroup()
  write.csv(top100_markers,
            file.path(results_dir, paste0(obj_name, "_top100_markers.csv")),
            row.names = FALSE)

  # Gene ID conversion / 基因ID转换
  unique_genes <- unique(group_degs_sig$gene)
  valid_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  valid_genes <- unique_genes[unique_genes %in% valid_keys]
  ids <- bitr(valid_genes, fromType = 'SYMBOL', toType = 'ENTREZID',
              OrgDb = 'org.Hs.eg.db')

  group_degs_sig <- merge(group_degs_sig, ids, by.x = 'gene', by.y = 'SYMBOL')
  group_degs_sig <- group_degs_sig[order(group_degs_sig$avg_log2FC, decreasing = TRUE), ]

  # GO enrichment / GO富集分析
  go_genes <- group_degs_sig %>%
    filter(abs(avg_log2FC) > 1) %>%
    pull(ENTREZID)

  if (length(go_genes) > 0) {
    go_result <- enrichGO(
      go_genes,
      OrgDb = "org.Hs.eg.db",
      ont = "All",
      readable = TRUE
    )
    pdf(file.path(results_dir, paste0(obj_name, "_GO_dotplot.pdf")),
        width = 12, height = 14)
    print(dotplot(go_result, showCategory = 10, split = "ONTOLOGY",
                  title = paste0("TED vs normal ", obj_name, " GO")) +
          facet_grid(ONTOLOGY ~ ., scale = 'free'))
    dev.off()
    write.csv(as.data.frame(go_result),
              file.path(results_dir, paste0(obj_name, "_GO_results.csv")),
              row.names = FALSE)
  }

  # KEGG enrichment / KEGG富集分析
  if (length(go_genes) > 0) {
    kegg_result <- enrichKEGG(
      gene = go_genes,
      organism = "hsa",
      pvalueCutoff = 0.05
    )
    if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
      pdf(file.path(results_dir, paste0(obj_name, "_KEGG_dotplot.pdf")),
          width = 8, height = 5)
      print(dotplot(kegg_result, showCategory = 10,
                    title = paste0("TED vs normal ", obj_name, " KEGG")))
      dev.off()
      write.csv(as.data.frame(kegg_result),
                file.path(results_dir, paste0(obj_name, "_KEGG_results.csv")),
                row.names = FALSE)
    }
  }

  # Celltype-specific GO / 细胞类型特异性GO
  if ("cluster" %in% colnames(top100_markers)) {
    bp_result <- compareCluster(
      gene ~ cluster,
      data = top100_markers,
      fun = 'enrichGO',
      OrgDb = 'org.Hs.eg.db',
      keyType = 'SYMBOL',
      ont = "All"
    )
    pdf(file.path(results_dir, paste0(obj_name, "_celltype_GO.pdf")),
        width = 7, height = 12)
    print(dotplot(bp_result, showCategory = 10, font.size = 7) +
          facet_grid(ONTOLOGY ~ ., scale = 'free') +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    dev.off()
  }

  # Celltype-specific KEGG / 细胞类型特异性KEGG
  if ("cluster" %in% colnames(top100_markers)) {
    marker_ids <- bitr(top100_markers$gene, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    top100_markers_km <- left_join(top100_markers, marker_ids,
                                   by = c('gene' = "SYMBOL"))
    kegg_cluster <- compareCluster(
      ENTREZID ~ cluster,
      data = top100_markers_km,
      fun = "enrichKEGG",
      organism = 'hsa',
      pvalueCutoff = 0.05
    )
    if (!is.null(kegg_cluster)) {
      pdf(file.path(results_dir, paste0(obj_name, "_celltype_KEGG.pdf")),
          width = 12, height = 10)
      print(dotplot(kegg_cluster, showCategory = 10, font.size = 10) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)))
      dev.off()
    }
  }

  # GSEA (KEGG) / GSEA分析
  # Prepare ranked gene list for GSEA / 准备GSEA排序基因列表
  group_degs_gsea <- group_degs %>%
    mutate(gene = rownames(.))

  ids_gsea <- bitr(group_degs_gsea$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db')
  group_degs_gsea <- merge(group_degs_gsea, ids_gsea,
                            by.x = 'gene', by.y = 'SYMBOL')
  group_degs_gsea <- group_degs_gsea[order(group_degs_gsea$avg_log2FC,
                                             decreasing = TRUE), ]

  gsea_gene_list <- group_degs_gsea$avg_log2FC
  names(gsea_gene_list) <- group_degs_gsea$ENTREZID

  gsea_kegg <- gseKEGG(
    gsea_gene_list,
    organism = "hsa",
    pvalueCutoff = 0.05
  )

  if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
    pdf(file.path(results_dir, paste0(obj_name, "_GSEA_KEGG.pdf")),
        width = 15, height = 8)
    print(gseaplot2(gsea_kegg, 1:min(10, nrow(gsea_kegg)),
                    pvalue_table = FALSE, base_size = 14))
    dev.off()
    write.csv(as.data.frame(gsea_kegg),
              file.path(results_dir, paste0(obj_name, "_GSEA_KEGG_results.csv")),
              row.names = FALSE)
  }
}

# Run enrichment for each subcluster / 对各亚群执行富集分析
run_enrichment(Fibsub, "Fibsub")
run_enrichment(B_sub, "B_sub")
run_enrichment(Endothelialsub, "Endothelialsub")
run_enrichment(Myeloidsub, "Myeloidsub")
run_enrichment(NK_Tsub, "NK_Tsub")
run_enrichment(Pericyte_sub, "Pericyte_sub")
run_enrichment(Myo_sub, "Myo_sub")
run_enrichment(Satellite_sub, "Satellite_sub")
run_enrichment(SMC_sub, "SMC_sub")

# ============================================================================
# SECTION 9: Fibroblast Pseudotime Analysis (Monocle2) / 成纤维细胞拟时序分析
# ============================================================================
library(monocle)
library(Biobase)
library(viridis)

# 9.1 Sample cells for pseudotime analysis / 采样用于拟时序分析
set.seed(42)
sample_size <- 8000
cells_to_sample <- Cells(Fibsub)
if (length(cells_to_sample) > sample_size) {
  cells_to_sample <- sample(cells_to_sample, size = sample_size)
}
random_subset <- subset(Fibsub, cells = cells_to_sample)

# Set identities to celltype labels (Fib1-Fib5)
Idents(random_subset) <- random_subset$celltype

# Subset to include only Fib1-Fib5 subtypes
seurat_pseudo <- subset(random_subset, idents = c("Fib1", "Fib2", "Fib3", "Fib4", "Fib5"))

# 9.2 Create CellDataSet object / 创建CellDataSet对象
expr_matrix <- as.matrix(GetAssayData(seurat_pseudo, assay = "RNA", slot = "counts"))
sample_sheet <- seurat_pseudo@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(seurat_pseudo))
rownames(gene_annotation) <- rownames(seurat_pseudo)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

cds <- newCellDataSet(
  expr_matrix,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

# 9.3 Differential gene test by celltype / 按celltype差异基因检测
# Use celltype (Fib1-Fib5) as the ordering variable
diff_celltype <- differentialGeneTest(
  cds[expressed_genes, ],
  fullModelFormulaStr = "~celltype",
  cores = 4
)
write.csv(diff_celltype, file.path(results_dir, "pseudotime_degForCellOrdering.csv"))

# Select top 1000 ordering genes
diff_celltype <- diff_celltype[order(diff_celltype$qval), ]
ordering_genes <- row.names(diff_celltype)[1:min(1000, nrow(diff_celltype))]
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)

pseudotime_dir <- file.path(results_dir, "pseudotime")
dir.create(pseudotime_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(pseudotime_dir, "monocle2_ordering_genes.pdf"), width = 8, height = 6)
plot_ordering_genes(cds)
dev.off()

# 9.4 Dimensionality reduction & ordering / DDRTree降维与排序
cds <- reduceDimension(cds, method = "DDRTree")
cds <- orderCells(cds)

# 9.5 Set root state (Fib1 as root) / 设置起点（Fib1为根节点）
# Define helper function to find the state containing the most Fib1 cells
GM_state <- function(cds) {
  if (length(unique(pData(cds)$State)) > 1) {
    T0_counts <- table(pData(cds)$State, pData(cds)$celltype)[, "Fib1"]
    return(as.numeric(names(T0_counts)[which.max(T0_counts)]))
  } else {
    return(1)
  }
}

cds <- orderCells(cds, root_state = GM_state(cds), reverse = TRUE)

# 9.6 Trajectory visualization / 拟时序轨迹可视化
# State trajectory
p1 <- plot_cell_trajectory(cds, color_by = "State") +
  theme(text = element_text(size = 18))
ggsave(file.path(pseudotime_dir, "monocle2_state_trajectory.pdf"),
       plot = p1, width = 12, height = 9)

# Celltype trajectory (Fib1-Fib5)
p2 <- plot_cell_trajectory(cds, color_by = "celltype") +
  theme(text = element_text(size = 18))
ggsave(file.path(pseudotime_dir, "monocle2_celltype_trajectory.pdf"),
       plot = p2, width = 12, height = 9)

# Pseudotime trajectory
p3 <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  theme(text = element_text(size = 18))
ggsave(file.path(pseudotime_dir, "monocle2_Pseudotime.pdf"),
       plot = p3, width = 12, height = 9)

# Combined plot (celltype + pseudotime)
ggsave(file.path(pseudotime_dir, "monocle2_combined.pdf"),
       plot = p2 + p3, width = 20, height = 9)

# Facet by celltype
p_facet <- plot_cell_trajectory(cds, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 2) +
  theme(text = element_text(size = 18))
ggsave(file.path(pseudotime_dir, "monocle2_celltype_facet.pdf"),
       plot = p_facet, width = 16, height = 10)

# 9.7 Pseudotime-dependent genes / 拟时序依赖基因
diff_test_res <- differentialGeneTest(
  cds[expressed_genes[1:min(1000, length(expressed_genes))], ],
  fullModelFormulaStr = "~sm.ns(Pseudotime)",
  cores = 4
)
write.csv(diff_test_res, file.path(results_dir, "pseudotime_diff_test_res.csv"))

# Top 5 pseudotime-dependent genes visualization
top5g <- rownames(diff_test_res[order(diff_test_res$qval), ])[1:5]
p_top5 <- plot_genes_in_pseudotime(cds[top5g, ], color_by = "celltype", ncol = 1) +
  theme(text = element_text(size = 14))
ggsave(file.path(pseudotime_dir, "monocle2_top5_pseudotime.pdf"),
       plot = p_top5, width = 10, height = 12)

# 9.8 Pseudotime heatmap (top 20 genes) / 拟时序热图
diff_test_res <- diff_test_res[order(diff_test_res$qval), ]
top20g <- rownames(diff_test_res[1:20, ])

pdf(file.path(pseudotime_dir, "monocle2_pseudotime_heatmap.pdf"),
    width = 14, height = 10)
plot_pseudotime_heatmap(
  cds[top20g, ],
  num_clusters = 6,
  show_rownames = TRUE,
  hmcols = colorRampPalette(viridis(10))(100)
)
dev.off()

# Save pseudotime analysis objects
save(cds, ordering_genes, diff_test_res,
     file = file.path(results_dir, "pseudotime_analysis_final.RData"))

cat("Fibroblast pseudotime analysis (Monocle2) complete.\n")

# ============================================================================
# SECTION 9B: Endothelial Pseudotime Analysis (EndoMT trajectory) / 内皮细胞拟时序分析
# ============================================================================
# Reconstructs the Endo1 -> Endo4 EndoMT trajectory used in Figure 3.
# Root state is defined as the state containing the most Endo1 cells
# (Endo1 = quiescent capillary endothelium, start of trajectory).

# 9B.1 Sample cells for pseudotime analysis / 采样用于拟时序分析
set.seed(42)
sample_size_endo <- 8000
if (ncol(Endothelialsub) > sample_size_endo) {
  cell_indices_endo <- sample(Cells(Endothelialsub), size = sample_size_endo)
} else {
  cell_indices_endo <- Cells(Endothelialsub)
}
random_subset_endo <- subset(Endothelialsub, cells = cell_indices_endo)

# Set identities to celltype labels (Endo1-Endo4)
Idents(random_subset_endo) <- random_subset_endo$celltype

# Subset to include only Endo1-Endo4 subtypes
seurat_pseudo_endo <- subset(random_subset_endo,
                             idents = c("Endo1", "Endo2", "Endo3", "Endo4"))

# 9B.2 Create CellDataSet object / 创建CellDataSet对象
expr_matrix_endo <- as.matrix(GetAssayData(seurat_pseudo_endo,
                                          assay = "RNA", slot = "counts"))
sample_sheet_endo <- seurat_pseudo_endo@meta.data
gene_annotation_endo <- data.frame(gene_short_name = rownames(seurat_pseudo_endo))
rownames(gene_annotation_endo) <- rownames(seurat_pseudo_endo)

pd_endo <- new("AnnotatedDataFrame", data = sample_sheet_endo)
fd_endo <- new("AnnotatedDataFrame", data = gene_annotation_endo)

cds_endo <- newCellDataSet(
  expr_matrix_endo,
  phenoData = pd_endo,
  featureData = fd_endo,
  expressionFamily = negbinomial.size()
)

cds_endo <- estimateSizeFactors(cds_endo)
cds_endo <- estimateDispersions(cds_endo)
cds_endo <- detectGenes(cds_endo, min_expr = 0.1)
expressed_genes_endo <- row.names(subset(fData(cds_endo),
                                         num_cells_expressed >= 10))

# 9B.3 Differential gene test by celltype / 按celltype差异基因检测
diff_celltype_endo <- differentialGeneTest(
  cds_endo[expressed_genes_endo, ],
  fullModelFormulaStr = "~celltype",
  cores = 4
)
write.csv(diff_celltype_endo,
          file.path(results_dir, "pseudotime_endo_degForCellOrdering.csv"))

# Select top 1000 ordering genes
diff_celltype_endo <- diff_celltype_endo[order(diff_celltype_endo$qval), ]
ordering_genes_endo <- row.names(diff_celltype_endo)[1:min(1000, nrow(diff_celltype_endo))]
cds_endo <- setOrderingFilter(cds_endo, ordering_genes = ordering_genes_endo)

pseudotime_endo_dir <- file.path(results_dir, "pseudotime_endo")
dir.create(pseudotime_endo_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(pseudotime_endo_dir, "monocle2_endo_ordering_genes.pdf"),
    width = 8, height = 6)
plot_ordering_genes(cds_endo)
dev.off()

# 9B.4 Dimensionality reduction & ordering / DDRTree降维与排序
cds_endo <- reduceDimension(cds_endo, method = "DDRTree")
cds_endo <- orderCells(cds_endo)

# 9B.5 Set root state (Endo1 as root) / 设置起点（Endo1为根节点）
# Define helper function to find the state containing the most Endo1 cells
# (quiescent capillary endothelium = trajectory origin)
endo_root_state <- function(cds_obj) {
  if (length(unique(pData(cds_obj)$State)) > 1) {
    T0_counts_endo <- table(pData(cds_obj)$State, pData(cds_obj)$celltype)[, "Endo1"]
    return(as.numeric(names(T0_counts_endo)[which.max(T0_counts_endo)]))
  } else {
    return(1)
  }
}

cds_endo <- orderCells(cds_endo, root_state = endo_root_state(cds_endo))

# 9B.6 Trajectory visualization / 拟时序轨迹可视化
# State trajectory
p1_endo <- plot_cell_trajectory(cds_endo, color_by = "State") +
  theme(text = element_text(size = 18))
ggsave(file.path(pseudotime_endo_dir, "monocle2_endo_state_trajectory.pdf"),
       plot = p1_endo, width = 12, height = 9)

# Celltype trajectory (Endo1-Endo4)
p2_endo <- plot_cell_trajectory(cds_endo, color_by = "celltype") +
  theme(text = element_text(size = 18))
ggsave(file.path(pseudotime_endo_dir, "monocle2_endo_celltype_trajectory.pdf"),
       plot = p2_endo, width = 12, height = 9)

# Pseudotime trajectory
p3_endo <- plot_cell_trajectory(cds_endo, color_by = "Pseudotime") +
  theme(text = element_text(size = 18))
ggsave(file.path(pseudotime_endo_dir, "monocle2_endo_Pseudotime.pdf"),
       plot = p3_endo, width = 12, height = 9)

# Combined plot (celltype + pseudotime)
ggsave(file.path(pseudotime_endo_dir, "monocle2_endo_combined.pdf"),
       plot = p2_endo + p3_endo, width = 20, height = 9)

# Facet by celltype
p_facet_endo <- plot_cell_trajectory(cds_endo, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 2) +
  theme(text = element_text(size = 18))
ggsave(file.path(pseudotime_endo_dir, "monocle2_endo_celltype_facet.pdf"),
       plot = p_facet_endo, width = 16, height = 10)

# 9B.7 Pseudotime-dependent genes / 拟时序依赖基因
diff_test_res_endo <- differentialGeneTest(
  cds_endo[expressed_genes_endo[1:min(1000, length(expressed_genes_endo))], ],
  fullModelFormulaStr = "~sm.ns(Pseudotime)",
  cores = 4
)
write.csv(diff_test_res_endo,
          file.path(results_dir, "pseudotime_endo_diff_test_res.csv"))

# Top 5 pseudotime-dependent genes visualization
top5g_endo <- rownames(diff_test_res_endo[order(diff_test_res_endo$qval), ])[1:5]
p_top5_endo <- plot_genes_in_pseudotime(cds_endo[top5g_endo, ],
                                        color_by = "celltype", ncol = 1) +
  theme(text = element_text(size = 14))
ggsave(file.path(pseudotime_endo_dir, "monocle2_endo_top5_pseudotime.pdf"),
       plot = p_top5_endo, width = 10, height = 12)

# 9B.8 Pseudotime heatmap (top 20 genes) / 拟时序热图
diff_test_res_endo <- diff_test_res_endo[order(diff_test_res_endo$qval), ]
top20g_endo <- rownames(diff_test_res_endo[1:20, ])

pdf(file.path(pseudotime_endo_dir, "monocle2_endo_pseudotime_heatmap.pdf"),
    width = 14, height = 10)
plot_pseudotime_heatmap(
  cds_endo[top20g_endo, ],
  num_clusters = 6,
  show_rownames = TRUE,
  hmcols = colorRampPalette(viridis(10))(100)
)
dev.off()

# Save endothelial pseudotime analysis objects
save(cds_endo, ordering_genes_endo, diff_test_res_endo,
     file = file.path(results_dir, "pseudotime_endo_analysis_final.RData"))

cat("Endothelial pseudotime analysis (EndoMT trajectory) complete.\n")

# ============================================================================
# SECTION 10: CellChat Cell Communication Analysis / CellChat细胞通讯分析
# ============================================================================
library(CellChat)
library(NMF)
library(ggalluvial)
library(pheatmap)
library(cowplot)
library(gridExtra)

# 10.1 Data preparation / 数据准备
# Remove epi and neuro cell types for CellChat (as in original analysis)
ifnb <- subset(scRNA_harmony, celltype != "epi")
ifnb <- subset(ifnb, celltype != "neuro")

# Use integrated subcluster labels (Fib1-5, Endo1-4, Mac/Mono/DC) from Section 7
# so that CellChat can reference subcluster-level sender/receiver groups below
ifnb$celltype <- factor(ifnb$integrated_celltype)

# Split by group / 按分组拆分
ifnb.list <- SplitObject(ifnb, split.by = "group")

# Normalize and find variable features for integration
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Integration / 数据整合
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(
  object.list = ifnb.list,
  anchor.features = features,
  reduction = "rpca",
  k.anchor = 20
)
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
immune.combined <- RunHarmony(immune.combined, group.by.vars = "orig.ident")
immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:50)
immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:50)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
DefaultAssay(immune.combined) <- "RNA"
Idents(immune.combined) <- "celltype"

saveRDS(immune.combined, file.path(results_dir, "immune.combined.rds"))

# 10.2 TED group CellChat / TED组CellChat分析
stim.object <- subset(immune.combined, group == "TED")
stim.data.input <- GetAssayData(stim.object, assay = "RNA", slot = "data")
stim.meta <- stim.object@meta.data[, c("celltype", "group")]

stim.cellchat <- createCellChat(object = stim.data.input)
stim.cellchat <- addMeta(stim.cellchat, meta = stim.meta)
stim.cellchat <- setIdent(stim.cellchat, ident.use = "celltype")
stim.cellchat@DB <- CellChatDB.human

stim.cellchat <- subsetData(stim.cellchat, features = NULL)
future::plan("multisession", workers = 4)
stim.cellchat <- identifyOverExpressedGenes(stim.cellchat)
stim.cellchat <- identifyOverExpressedInteractions(stim.cellchat)
stim.cellchat <- projectData(stim.cellchat, PPI.human)

stim.cellchat <- computeCommunProb(stim.cellchat, raw.use = TRUE)
stim.cellchat <- filterCommunication(stim.cellchat, min.cells = 10)
stim.cellchat <- computeCommunProbPathway(stim.cellchat)
stim.cellchat <- aggregateNet(stim.cellchat)
stim.cellchat <- netAnalysis_computeCentrality(stim.cellchat, slot.name = "netP")

group1.net <- subsetCommunication(stim.cellchat)
write.csv(group1.net,
          file.path(results_dir, "TED_cellchat_interactions.csv"),
          row.names = FALSE)
saveRDS(stim.cellchat, file.path(results_dir, "stim.cellchat.rds"))

# 10.3 Normal group CellChat / Normal组CellChat分析
ctrl.object <- subset(immune.combined, group == "normal")
ctrl.data.input <- GetAssayData(ctrl.object, assay = "RNA", slot = "data")
ctrl.meta <- ctrl.object@meta.data[, c("celltype", "group")]

ctrl.cellchat <- createCellChat(object = ctrl.data.input)
ctrl.cellchat <- addMeta(ctrl.cellchat, meta = ctrl.meta)
ctrl.cellchat <- setIdent(ctrl.cellchat, ident.use = "celltype")
ctrl.cellchat@DB <- CellChatDB.human

ctrl.cellchat <- subsetData(ctrl.cellchat)
future::plan("multisession", workers = 4)
ctrl.cellchat <- identifyOverExpressedGenes(ctrl.cellchat)
ctrl.cellchat <- identifyOverExpressedInteractions(ctrl.cellchat)
ctrl.cellchat <- projectData(ctrl.cellchat, PPI.human)

ctrl.cellchat <- computeCommunProb(ctrl.cellchat)
ctrl.cellchat <- filterCommunication(ctrl.cellchat, min.cells = 10)
ctrl.cellchat <- computeCommunProbPathway(ctrl.cellchat)
ctrl.cellchat <- aggregateNet(ctrl.cellchat)
ctrl.cellchat <- netAnalysis_computeCentrality(ctrl.cellchat, slot.name = "netP")

saveRDS(ctrl.cellchat, file.path(results_dir, "ctrl.cellchat.rds"))

# 10.4 Comparison & Visualization / 比较与可视化
object.list <- list(CTRL = ctrl.cellchat, STIM = stim.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat_dir <- file.path(results_dir, "cellchat_comparison")
dir.create(cellchat_dir, showWarnings = FALSE, recursive = TRUE)

# Interaction count & strength comparison / 互作次数与强度比较
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2),
                           measure = "weight")
ggsave(file.path(cellchat_dir, "interaction_comparison_bar.pdf"),
       plot = gg1 + gg2, width = 12, height = 5)

# Difference network / 差异网络图
pdf(file.path(cellchat_dir, "interaction_diff_network.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2), xpd = TRUE, mar = c(1, 1, 3, 1))
netVisual_diffInteraction(cellchat, weight.scale = TRUE,
                          title = "Number of interactions")
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight",
                          title = "Interaction weights/strength")
dev.off()

# Pathway comparison / 信号通路比较
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
ggsave(file.path(cellchat_dir, "pathway_comparison.pdf"),
       plot = gg1 + gg2, width = 14, height = 6)

# 10.5 Expanded CellChat Comparison Visualizations / 扩展CellChat对比可视化
library(ComplexHeatmap)

# 10.5.1 Cell type interaction count circles / 细胞类型互作次数环形图
pdf(file.path(cellchat_dir, "celltype_interaction_counts_circle.pdf"),
    width = 12, height = 6)
par(mfrow = c(1, 2))
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE,
                   label.edge = FALSE,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# 10.5.2 Selected cell types interaction (Fib1-5 + Pericyte) / 特定细胞类型互作
s.cell <- c("Fib1", "Fib2", "Fib3", "Fib4", "Fib5", "Pericyte")
count1 <- object.list[[1]]@net$count[s.cell, s.cell]
count2 <- object.list[[2]]@net$count[s.cell, s.cell]
weight.max_select <- max(max(count1), max(count2))

pdf(file.path(cellchat_dir, "selected_celltype_Fib_Pericyte_interaction.pdf"),
    width = 10, height = 6.5)
par(mfrow = c(1, 2))
netVisual_circle(count1, weight.scale = TRUE, label.edge = TRUE,
                 edge.weight.max = weight.max_select, edge.width.max = 12,
                 title.name = paste0("Interactions - ", names(object.list)[1]))
netVisual_circle(count2, weight.scale = TRUE, label.edge = TRUE,
                 edge.weight.max = weight.max_select, edge.width.max = 12,
                 title.name = paste0("Interactions - ", names(object.list)[2]))
dev.off()

# 10.5.3 Interaction difference heatmap / 互作差异热图
if (all(c("STIM", "CTRL") %in% names(cellchat@net))) {
  diff.count <- cellchat@net$STIM$count - cellchat@net$CTRL$count
  write.csv(cellchat@net$STIM$count,
            file.path(cellchat_dir, "STIM_interaction_count.csv"), quote = FALSE)
  write.csv(cellchat@net$CTRL$count,
            file.path(cellchat_dir, "CTRL_interaction_count.csv"), quote = FALSE)

  pdf(file.path(cellchat_dir, "interaction_diff_heatmap.pdf"), width = 10, height = 8)
  pheatmap(diff.count, treeheight_row = 0, treeheight_col = 0)
  dev.off()
}

# 10.5.4 Signaling role heatmap (ComplexHeatmap) / 信号角色热图
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

# Outgoing signaling role
ht1_out <- netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "outgoing", signaling = pathway.union,
  title = names(object.list)[1], width = 5, height = 6)
ht2_out <- netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "outgoing", signaling = pathway.union,
  title = names(object.list)[2], width = 5, height = 6)

pdf(file.path(cellchat_dir, "signaling_role_outgoing.pdf"), width = 10, height = 6)
ComplexHeatmap::draw(ht1_out + ht2_out, ht_gap = unit(0.5, "cm"))
dev.off()

# All signaling role
ht1_all <- netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "all", signaling = pathway.union,
  title = names(object.list)[1], width = 8, height = 10)
ht2_all <- netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "all", signaling = pathway.union,
  title = names(object.list)[2], width = 8, height = 10)

pdf(file.path(cellchat_dir, "signaling_role_all.pdf"), width = 10, height = 6)
ComplexHeatmap::draw(ht1_all + ht2_all, ht_gap = unit(0.5, "cm"))
dev.off()

# Incoming signaling role
ht1_in <- netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "incoming", signaling = pathway.union,
  title = names(object.list)[1], width = 8, height = 10)
ht2_in <- netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "incoming", signaling = pathway.union,
  title = names(object.list)[2], width = 8, height = 10)

pdf(file.path(cellchat_dir, "signaling_role_incoming.pdf"), width = 10, height = 6)
ComplexHeatmap::draw(ht1_in + ht2_in, ht_gap = unit(0.5, "cm"))
dev.off()

# 10.5.5 Pathway similarity analysis / 通路相似性分析
# Note: functional similarity manifold may not render; structural works reliably
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

p_sim <- rankSimilarity(cellchat, type = "structural") +
  ggtitle("Structural similarity of pathway")
ggsave(file.path(cellchat_dir, "pathway_similarity.pdf"),
       plot = p_sim, width = 8, height = 5)

# 10.5.6 MK signaling pathway specific comparison / MK信号通路特异性对比
pathways.show <- "MK"

# MK network (circle layout)
weight.max_mk <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
pdf(file.path(cellchat_dir, "MK_pathway_network.pdf"), width = 10, height = 6.5)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max_mk[1], edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# MK heatmap
pdf(file.path(cellchat_dir, "MK_pathway_heatmap.pdf"), width = 12, height = 6.5)
par(mfrow = c(1, 2), xpd = TRUE)
ht_mk <- list()
for (i in 1:length(object.list)) {
  ht_mk[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show,
                                  color.heatmap = "Reds",
                                  title.name = paste(pathways.show, "signaling", names(object.list)[i]))
}
ComplexHeatmap::draw(ht_mk[[1]] + ht_mk[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

# MK chord diagram
pdf(file.path(cellchat_dir, "MK_pathway_chord.pdf"), width = 10, height = 6.5)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord",
                      pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# 10.5.7 MK gene expression visualization / MK基因表达可视化
cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("STIM", "CTRL"))

pdf(file.path(cellchat_dir, "MK_geneExpression_STIM_CTRL.pdf"), width = 10, height = 8)
plotGeneExpression(cellchat, signaling = "MK", split.by = "datasets", colors.ggplot = TRUE)
dev.off()

# MK pathway: MDK gene only
pdf(file.path(cellchat_dir, "MK_MDK_geneExpression.pdf"), width = 10, height = 6)
plotGeneExpression(cellchat, signaling = "MK", split.by = "datasets",
                   colors.ggplot = TRUE, features = "MDK")
dev.off()

# 10.5.8 Ligand-receptor bubble comparison (MK) / 配体-受体气泡对比图(MK)
# Endo1-4 + Fib1-5 as both sources and targets
p_mk_bubble <- netVisual_bubble(
  cellchat,
  sources.use = c("Endo1", "Endo2", "Endo3", "Endo4",
                   "Fib1", "Fib2", "Fib3", "Fib4", "Fib5"),
  targets.use = c("Endo1", "Endo2", "Endo3", "Endo4",
                   "Fib1", "Fib2", "Fib3", "Fib4", "Fib5"),
  comparison = c(1, 2),
  angle.x = 45,
  signaling = c("MK"),
  thresh = 0.0001
)

# Filter columns with no signal
keep_cols <- sapply(unique(p_mk_bubble$data$x_label), function(col) {
  sub <- p_mk_bubble$data[p_mk_bubble$data$x_label == col, ]
  any(sub$prob > 0.0001)
})
keep_names <- unique(p_mk_bubble$data$x_label)[keep_cols]
p_mk_bubble$data <- p_mk_bubble$data[p_mk_bubble$data$x_label %in% keep_names, ]
p_mk_bubble <- p_mk_bubble + scale_radius(limits = c(0, NA))

ggsave(file.path(cellchat_dir, "MK_LR_bubble_comparison.pdf"),
       plot = p_mk_bubble, width = 14, height = 9)

# Save CellChat merged object
saveRDS(cellchat, file.path(cellchat_dir, "cellchat_merged.rds"))

# ============================================================================
# SECTION 11: Final Save / 最终保存
# ============================================================================
# Save all major objects / 保存所有主要对象
save(
  scRNA_harmony,
  Fibsub, B_sub, Endothelialsub, Myeloidsub, NK_Tsub,
  Pericyte_sub, Myo_sub, Satellite_sub, SMC_sub,
  celltype_proportions, fisher_results,
  celltype_anno,
  cds, ordering_genes, diff_test_res,
  cds_endo, ordering_genes_endo, diff_test_res_endo,
  immune.combined,
  stim.cellchat, ctrl.cellchat, cellchat,
  file = file.path(results_dir, "TAO_analysis_final.RData")
)

cat("\n=== Analysis Complete ===\n")
cat("All results saved to:", results_dir, "\n")
cat("Key output files:\n")
cat("  - Fibroblast pseudotime: results/pseudotime/\n")
cat("  - Endothelial pseudotime (EndoMT): results/pseudotime_endo/\n")
cat("  - CellChat comparison: results/cellchat_comparison/\n")
cat("  - Subcluster integration: scRNA_harmony_subclusters_integrated.pdf\n")

# Record software environment for reproducibility / 记录软件环境
sessionInfo()

# ============================================================================
# END OF SCRIPT
# ============================================================================
