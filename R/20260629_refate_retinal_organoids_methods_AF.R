# scRNA-seq analysis for pilocarpine-treated and control retinal organoids

# load libraries

library(Seurat)
library(SeuratObject)
library(scClassify)
library(ggplot2)
library(dplyr)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(BiocParallel)


# load organoid data

## control sample

matrix_data <- readMM("~/refate_retinal_organoids/Anh_GEMX_FLEX_output/outs/per_sample_outs/sample1/sample_filtered_feature_bc_matrix/matrix.mtx")
feature_names <- read.delim("~/refate_retinal_organoids/Anh_GEMX_FLEX_output/outs/per_sample_outs/sample1/sample_filtered_feature_bc_matrix/features.tsv", header = FALSE, stringsAsFactors = FALSE)
barcode_names <- read.delim("~/refate_retinal_organoids/Anh_GEMX_FLEX_output/outs/per_sample_outs/sample1/sample_filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)

colnames(matrix_data) <- barcode_names$V1
rownames(matrix_data) <- feature_names$V2

sample1 <- CreateSeuratObject(counts = matrix_data, assay = "RNA")
sample1$condition <- "Control_10wk"

## pilocarpine-treated sample

matrix_data <- readMM("~/refate_retinal_organoids/Anh_GEMX_FLEX_output/outs/per_sample_outs/sample2/sample_filtered_feature_bc_matrix/matrix.mtx")
feature_names <- read.delim("~/refate_retinal_organoids/Anh_GEMX_FLEX_output/outs/per_sample_outs/sample2/sample_filtered_feature_bc_matrix/features.tsv", header = FALSE, stringsAsFactors = FALSE)
barcode_names <- read.delim("~/refate_retinal_organoids/Anh_GEMX_FLEX_output/outs/per_sample_outs/sample2/sample_filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)

colnames(matrix_data) <- barcode_names$V1
rownames(matrix_data) <- feature_names$V2

sample2 <- CreateSeuratObject(counts = matrix_data, assay = "RNA") 
sample2$condition <- "Pilocarpine_10wk"


# combine samples

combined_seu <- merge(
  x = sample1,
  y = sample2,
  add.cell.ids = c("S1", "S2")
)
combined_seu <- JoinLayers(combined_seu)


# filter out low quality cells

mito.genes <- rownames(sample1)[grep("MT-", rownames(sample1))]
total_counts <- colSums(GetAssayData(combined_seu, layer = "counts"))
mito_counts <-colSums(GetAssayData(combined_seu, layer = "counts")[mito.genes, ]) 
combined_seu$mito_prop <- mito_counts / total_counts

combined_seu_filtered <- subset(
  combined_seu,
  subset =
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &   
    nFeature_RNA > 200 &
    nFeature_RNA < 10000 &
    mito_prop < 0.1
)


# RNA data processing

combined_seu_filtered <- NormalizeData(combined_seu_filtered)
combined_seu_filtered <- FindVariableFeatures(combined_seu_filtered)
combined_seu_filtered <- ScaleData(combined_seu_filtered)
combined_seu_filtered <- RunPCA(combined_seu_filtered) 
combined_seu_filtered <- RunUMAP(combined_seu_filtered, reduction = "pca", dims = 1:30, reduction.name = "umap") 


# plot UMAP

plot1 <- DimPlot(combined_seu_filtered, reduction = "umap", group.by = "condition")


# load and process retina reference datasets

set.seed(123)

## HANI RETINA ATLAS

retina_ref_hani <- readRDS("~/refate_retinal_organoids/human_retinal_final.rds")

# cap cell type counts
cells_to_keep <- c()
atlas.celltypes <- unique(retina_ref_hani$atlas_celltype_scRcls)

Idents(retina_ref_hani) <- retina_ref_hani$atlas_celltype_scRcls

for (celltype in atlas.celltypes) {
  
  cells_of_type <- WhichCells(retina_ref_hani, idents = celltype)
  
  if (length(cells_of_type) > 2000) {
    cells_of_type <- sample(cells_of_type, 2000)
  }
  cells_to_keep <- c(cells_to_keep, cells_of_type)
  
}
if (length(cells_to_keep) > 0) {
  retina_ref_hani_capped <- subset(retina_ref_hani, cells = cells_to_keep)
}

## HRCA RETINA ATLAS

retina_ref_hca <- readRDS("~/refate_retinal_organoids/seu_HCA_allcells.RDS")

# convert ensembl IDs to gene symbols
ensembl_ids <- rownames(retina_ref_hca)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

valid_idx <- !is.na(gene_symbols)
ensembl_ids <- ensembl_ids[valid_idx]
gene_symbols <- gene_symbols[valid_idx]

keep_idx <- !duplicated(gene_symbols)
ensembl_ids <- ensembl_ids[keep_idx]
gene_symbols <- gene_symbols[keep_idx]

ref_counts <- GetAssayData(retina_ref_hca, layer = "counts")
ref_counts <- ref_counts[ensembl_ids, ]

rownames(ref_counts) <- gene_symbols

metadata <- retina_ref_hca@meta.data
retina_ref_hca <- CreateSeuratObject(counts = ref_counts, meta.data = metadata)

# processing 
retina_ref_hca <- NormalizeData(retina_ref_hca)
retina_ref_hca <- FindVariableFeatures(retina_ref_hca)
retina_ref_hca <- ScaleData(retina_ref_hca)
retina_ref_hca <- RunPCA(retina_ref_hca)

# cap cell type counts
cells_to_keep <- c()
atlas.celltypes <- unique(retina_ref_hca$majorclass)

Idents(retina_ref_hca) <- retina_ref_hca$majorclass

for (celltype in atlas.celltypes) {
  
  cells_of_type <- WhichCells(retina_ref_hca, idents = celltype)
  
  if (length(cells_of_type) > 2000) {
    cells_of_type <- sample(cells_of_type, 2000)
  }
  cells_to_keep <- c(cells_to_keep, cells_of_type)
  
}
if (length(cells_to_keep) > 0) {
  retina_ref_hca_capped <- subset(retina_ref_hca, cells = cells_to_keep)
}

# DORGAU FOETAL RETINA ATLAS

retina_ref_dorgau <- readRDS("~/refate_retinal_organoids/seurat_list_dorgau_allsample_reference.RDS")

retina_ref_dorgau <- retina_ref_dorgau[grepl("Retina", names(retina_ref_dorgau))]

retina_ref_dorgau <- retina_ref_dorgau[[9]] # 10 weeks only
remove_these <- c(
  "Corneal stroma", 
  "Lens cells"
)
cells_to_keep <- rownames(retina_ref_dorgau@meta.data)[!(retina_ref_dorgau$CellType %in% remove_these)]
retina_ref_dorgau <- subset(retina_ref_dorgau, cells = cells_to_keep)

retina_ref_dorgau <- NormalizeData(retina_ref_dorgau)
retina_ref_dorgau <- FindVariableFeatures(retina_ref_dorgau)
retina_ref_dorgau <- ScaleData(retina_ref_dorgau)
retina_ref_dorgau <- RunPCA(retina_ref_dorgau)

# cap cell type counts
cells_to_keep <- c()
atlas.celltypes <- unique(retina_ref_dorgau$CellType)

Idents(retina_ref_dorgau) <- retina_ref_dorgau$CellType

for (celltype in atlas.celltypes) {
  
  cells_of_type <- WhichCells(retina_ref_dorgau, idents = celltype)
  
  if (length(cells_of_type) > 500) {
    cells_of_type <- sample(cells_of_type, 500)
  }
  cells_to_keep <- c(cells_to_keep, cells_of_type)
  
}
if (length(cells_to_keep) > 0) {
  retina_ref_dorgau_capped <- subset(retina_ref_dorgau, cells = cells_to_keep)
}


# train scClassify model on each capped reference 

hani_model_capped <- train_scClassify(
  exprsMat_train = retina_ref_hani_capped@assays$integrated$data,
  cellTypes_train = retina_ref_hani_capped$atlas_celltype_scRcls,
  tree = "HOPACH",
  selectFeatures = "limma",
  topN = 50,
  hopach_kmax = 5,
  pSig = 0.05,
  cellType_tree = NULL,
  weightsCal = FALSE,
  parallel = TRUE,
  BPPARAM = MulticoreParam(workers=15),
  verbose = TRUE,
  returnList = FALSE
)

hca_model_capped <- train_scClassify(
  exprsMat_train = retina_ref_hca_capped@assays$RNA$data,
  cellTypes_train = retina_ref_hca_capped$majorclass,
  tree = "HOPACH",
  selectFeatures = "limma",
  topN = 50,
  hopach_kmax = 5,
  pSig = 0.05,
  cellType_tree = NULL,
  weightsCal = FALSE,
  parallel = TRUE,
  BPPARAM = MulticoreParam(workers=15),
  verbose = TRUE,
  returnList = FALSE
)

dorgau_model_capped <- train_scClassify(
  exprsMat_train = retina_ref_dorgau_capped@assays$RNA$data,
  cellTypes_train = retina_ref_dorgau_capped$CellType,
  tree = "HOPACH",
  selectFeatures = "limma",
  topN = 50,
  hopach_kmax = 5,
  pSig = 0.05,
  cellType_tree = NULL,
  weightsCal = FALSE,
  parallel = TRUE,
  BPPARAM = MulticoreParam(workers=15),
  verbose = TRUE,
  returnList = FALSE
)


# run scClassify on the combined samples

hani_predictions <- predict_scClassify(
  combined_seu_filtered@assays$RNA$data,
  hani_model_capped,
  k = 10,
  prob_threshold = 0.7,
  cor_threshold_static = 0.5,
  cor_threshold_high = 0.7,
  features = "limma",
  algorithm = "WKNN",
  similarity = "pearson",
  cutoff_method = c("dynamic", "static"),
  parallel = TRUE,
  BPPARAM = MulticoreParam(workers=15),
  verbose = TRUE
)

combined_seu_filtered <- AddMetaData(combined_seu_filtered, metadata = hani_predictions$pearson_WKNN_limma$predLabelMat[,4], col.name = "scClassify.hani.capped.predicted.id")

hca_predictions <- predict_scClassify(
  combined_seu_filtered@assays$RNA$data,
  hca_model_capped,
  k = 10,
  prob_threshold = 0.7,
  cor_threshold_static = 0.5,
  cor_threshold_high = 0.7,
  features = "limma",
  algorithm = "WKNN",
  similarity = "pearson",
  cutoff_method = c("dynamic", "static"),
  parallel = TRUE,
  BPPARAM = MulticoreParam(workers=15),
  verbose = TRUE
)

combined_seu_filtered <- AddMetaData(combined_seu_filtered, metadata = hca_predictions$pearson_WKNN_limma$predLabelMat[,4], col.name = "scClassify.hca.capped.predicted.id")

dorgau_predictions <- predict_scClassify(
  combined_seu_filtered@assays$RNA$data,
  dorgau_model_capped,
  k = 10,
  prob_threshold = 0.7,
  cor_threshold_static = 0.5,
  cor_threshold_high = 0.7,
  features = "limma",
  algorithm = "WKNN",
  similarity = "pearson",
  cutoff_method = c("dynamic", "static"),
  parallel = TRUE,
  BPPARAM = MulticoreParam(workers=15),
  verbose = TRUE
)

combined_seu_filtered <- AddMetaData(combined_seu_filtered, metadata = dorgau_predictions$pearson_WKNN_limma$predLabelMat[,3], col.name = "scClassify.dorgau.pcw10.capped.predicted.id")


# harmonise cell type labels

combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("AC", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "Amacrine Cells"
combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("BC", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "Bipolar Cells"
combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("Cone", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "Cones"
combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("HC", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "Horizontal Cells"
combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("MG", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "Muller Glia"
combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("RGC", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "Retinal Ganglion Cells"
combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("Rod", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "Rods"

combined_seu_filtered$scClassify.hani.capped.predicted.id[grep("Microglia", combined_seu_filtered$scClassify.hani.capped.predicted.id)] <- "unassigned"
combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("Microglia", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "unassigned"
combined_seu_filtered$scClassify.hca.capped.predicted.id[grep("Astrocyte", combined_seu_filtered$scClassify.hca.capped.predicted.id)] <- "unassigned"


# re-annotate progenitor cells according to Dorgau reference

rpc_cells <- combined_seu_filtered[["scClassify.dorgau.pcw10.capped.predicted.id"]][, 1] == "RPCs"
combined_seu_filtered[["scClassify.hca.capped.predicted.id"]][rpc_cells, 1] <- "RPCs"
combined_seu_filtered[["scClassify.hani.capped.predicted.id"]][rpc_cells, 1] <- "RPCs"


# set colour palette

celltype_cols <- c(
  "Amacrine Cells" = "#BFE28A",
  "unassigned" = "#BDBDBD",
  "Muller Glia" = "#C9D8D0",
  "Cones" = "#F03B20",
  "Retinal Ganglion Cells" = "#B7E3EF",
  "Horizontal Cells" = "#F7B6C8",
  "Bipolar Cells" = "#006400",
  "RPE" = "#5B4BA3",
  "Rods" = "#E65AD9",
  "RPCs" = "#FFF59D"
)


# plot cell type proportions

df1 <- data.frame(
  Condition = combined_seu_filtered$condition,
  CellType = combined_seu_filtered$scClassify.hani.capped.predicted.id,
  reference = "Hani"
)
df2 <- data.frame(
  Condition = combined_seu_filtered$condition,
  CellType = combined_seu_filtered$scClassify.hca.capped.predicted.id,
  reference = "HCA"
)
df <- rbind(df1, df2)

df_plot <- df %>%
  mutate(
    combined = paste0(Condition, "_", reference)
  )

# set order
df_plot$combined <- factor(
  df_plot$combined,
  levels = c(
    "Control_10wk_Hani",
    "Control_10wk_HCA",
    "Pilocarpine_10wk_Hani",
    "Pilocarpine_10wk_HCA"
  )
)

df_plot$CellType <- factor(
  df_plot$CellType,
  levels = c("Cones","Rods","Horizontal Cells", "Amacrine Cells", "Muller Glia", "Bipolar Cells", "Retinal Ganglion Cells", "RPE", "RPCs", "unassigned"
  )
)

df_prop <- df_plot %>%
  dplyr::group_by(combined, CellType) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(combined) %>%
  dplyr::mutate(prop = n / sum(n))

plot2 <- ggplot(df_prop, aes(x = combined, y = prop, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.8) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    y = "Cell proportion",
    fill = "Cell type"
  ) +
  scale_fill_manual(values = celltype_cols, drop = FALSE)


# plot cell type counts

df_counts <- df_plot %>%
  dplyr::group_by(CellType, combined) %>%
  dplyr::summarise(n = n(), .groups = "drop")

plot3 <- ggplot(df_counts, aes(x = combined, y = n, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = celltype_cols, drop = FALSE)


# plot UMAPs coloured by cell type annotations

plot4 <- DimPlot(
  combined_seu_filtered,
  reduction = "umap",
  group.by = "scClassify.hani.capped.predicted.id",
  cols = celltype_cols
)

plot5 <- DimPlot(
  combined_seu_filtered,
  reduction = "umap",
  group.by = "scClassify.hca.capped.predicted.id",
  cols = celltype_cols
)


# plot marker expression dot plot

# main figure 5E markers:
genes <- c("CRX", "RXRG", "OTX2", "GNB3", "THRB", "GNGT2", "GNAT2")
# alternatively, supplementary figure 6G markers:
genes <- c("NR2E3", "NRL", "PDE6B")

expr <- FetchData(combined_seu_filtered, vars = genes)

df <- data.frame(
  condition = combined_seu_filtered$condition,
  expr
)
df_long <- df %>%
  pivot_longer(cols = all_of(genes),
               names_to = "gene",
               values_to = "expression")

df_long$expressing <- df_long$expression > 0

bubble_df <- df_long %>%
  dplyr::group_by(condition, gene) %>%
  dplyr::summarise(
    total_cells = n(),
    expressing_cells = sum(expressing),
    proportion = expressing_cells / total_cells,
    mean_expression = mean(expression),
    .groups = "drop"
  )

bubble_df$condition <- factor(bubble_df$condition, levels = c("Pilocarpine_10wk","Control_10wk"))

plot6 <- ggplot(bubble_df,
             aes(x = gene,
                 y = condition,
                 size = proportion,
                 color = mean_expression)) +
  geom_point() +
  scale_color_distiller(palette = "PuBu", direction = 1) +
  scale_size(range = c(2, 40)) +
  theme_classic()

