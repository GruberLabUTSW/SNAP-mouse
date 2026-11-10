#!/bin/bash
library(scDblFinder, lib.loc = conda_lib) 
library(dplyr)
library(hdf5r, lib.loc = conda_lib)
set.seed(42)
sample_name <- "Control"
sample_name1 <- "Cdki"
sample_name2 <- "CCL4"
sample_name3 <- "CCl4_Cdki"

cellranger_output <- ".../CONTROL/outs"
cellranger_output1 <- "...CDKi/outs"
cellranger_output2 <- ".../CCL4/outs"
cellranger_output3 <- ".../CCL4-CDK/outs"

counts <- Read10X(data.dir = file.path(cellranger_output, "raw_feature_bc_matrix"))
counts_cdki <- Read10X(data.dir = file.path(cellranger_output1, "raw_feature_bc_matrix"))
counts_ccl4 <- Read10X(data.dir = file.path(cellranger_output2, "raw_feature_bc_matrix"))
counts_ccl4_cdki <- Read10X(data.dir = file.path(cellranger_output3, "raw_feature_bc_matrix"))
cellbender_output <- ".../Cellbender_Output" 

counts_cellbender <- scCustomize::Read_CellBender_h5_Mat("Cellbender_Output/Control/control_cellbender_output_filtered.h5")
counts_cellbender_cdki <- scCustomize::Read_CellBender_h5_Mat("Cellbender_Output/CDKi/CDKi_cellbender_output_filtered.h5")
counts_cellbender_CCL4 <- scCustomize::Read_CellBender_h5_Mat("Cellbender_Output/CCL4/ccl4_cellbender_output_filtered.h5")
counts_cellbender_CCL4_CDK <- scCustomize::Read_CellBender_h5_Mat("Cellbender_Output/CCL4-CDK/ccl4-cdk_cellbender_output_filtered.h5")

seurat_obj <- CreateSeuratObject(counts = counts, project = "Control")
seurat_obj_cdki <- CreateSeuratObject(counts = counts_cdki, project = "CDKi")
seurat_obj_ccl4 <- CreateSeuratObject(counts = counts_ccl4, project = "CCL4")
seurat_obj_ccl4_cdki <- CreateSeuratObject(counts = counts_ccl4_cdki, project = "CCL4-CDK")

seurat_obj[["FILTERED"]] <- CreateAssayObject(counts = counts_cellbender)
seurat_obj_cdki[["FILTERED"]] <- CreateAssayObject(counts = counts_cellbender_cdki)
seurat_obj_ccl4[["FILTERED"]] <- CreateAssayObject(counts = counts_cellbender_CCL4)
seurat_obj_ccl4_cdki[["FILTERED"]] <- CreateAssayObject(counts = counts_cellbender_CCL4_CDK)
Assays(seurat_obj)

seurat_clean <- subset(seurat_obj, cells = colnames(counts_cellbender))
seurat_obj_cdki_clean      <- subset(seurat_obj_cdki, cells = colnames(counts_cellbender_cdki))
seurat_obj_ccl4_clean      <- subset(seurat_obj_ccl4, cells = colnames(counts_cellbender_CCL4))
seurat_obj_ccl4_cdki_clean <- subset(seurat_obj_ccl4_cdki, cells = colnames(counts_cellbender_CCL4_CDK)) 

DefaultAssay(seurat_clean) <- "FILTERED" 
DefaultAssay(seurat_obj_cdki_clean) <- "FILTERED"  
DefaultAssay(seurat_obj_ccl4_clean) <- "FILTERED"
DefaultAssay(seurat_obj_ccl4_cdki_clean) <- "FILTERED"

# At this point, data has only been cleaned with cellbender. Next, I am cleaning with doublet finder 
seurat_clean <- subset(seurat_obj, cells = colnames(counts_cellbender))
seurat_obj_cdki_clean      <- subset(seurat_obj_cdki, cells = colnames(counts_cellbender_cdki))
seurat_obj_ccl4_clean      <- subset(seurat_obj_ccl4, cells = colnames(counts_cellbender_CCL4))
seurat_obj_ccl4_cdki_clean <- subset(seurat_obj_ccl4_cdki, cells = colnames(counts_cellbender_CCL4_CDK)) 

DefaultAssay(seurat_clean) <- "FILTERED"
DefaultAssay(seurat_obj_cdki_clean) <- "FILTERED"
DefaultAssay(seurat_obj_ccl4_clean) <- "FILTERED"
DefaultAssay(seurat_obj_ccl4_cdki_clean) <- "FILTERED"

sce_control <- SingleCellExperiment(assays = list(counts = GetAssayData(seurat_clean, assay = "FILTERED", layer = "counts")))
sce_cdki <- SingleCellExperiment(assays = list(counts = GetAssayData(seurat_obj_cdki_clean, assay = "FILTERED", layer = "counts")))
sce_ccl4 <- SingleCellExperiment(assays = list(counts = GetAssayData(seurat_obj_ccl4_clean, assay = "FILTERED", layer = "counts")))
sce_ccl4_cdki <- SingleCellExperiment(assays = list(counts = GetAssayData(seurat_obj_ccl4_cdki_clean, assay = "FILTERED", layer = "counts")))

sce_control <- scDblFinder(sce_control)
sce_cdki <- scDblFinder(sce_cdki)
sce_ccl4 <- scDblFinder(sce_ccl4)
sce_ccl4_cdki <- scDblFinder(sce_ccl4_cdki)

seurat_clean$scDblFinder.class <- sce_control$scDblFinder.class[colnames(seurat_clean)]
seurat_obj_cdki_clean$scDblFinder.class <- sce_cdki$scDblFinder.class[colnames(seurat_obj_cdki_clean)]
seurat_obj_ccl4_clean$scDblFinder.class <- sce_ccl4$scDblFinder.class[colnames(seurat_obj_ccl4_clean)]
seurat_obj_ccl4_cdki_clean$scDblFinder.class <- sce_ccl4_cdki$scDblFinder.class[colnames(seurat_obj_ccl4_cdki_clean)]

seurat_clean <- subset(seurat_clean, subset = scDblFinder.class == "singlet")
seurat_obj_cdki_clean <- subset(seurat_obj_cdki_clean, subset = scDblFinder.class == "singlet")
seurat_obj_ccl4_clean <- subset(seurat_obj_ccl4_clean, subset = scDblFinder.class == "singlet")
seurat_obj_ccl4_cdki_clean <- subset(seurat_obj_ccl4_cdki_clean, subset = scDblFinder.class == "singlet")

dbl_class_control <- setNames(sce_control$scDblFinder.class, colnames(sce_control))
dbl_class_cdki <- setNames(sce_cdki$scDblFinder.class, colnames(sce_cdki))
dbl_class_ccl4 <- setNames(sce_ccl4$scDblFinder.class, colnames(sce_ccl4))
dbl_class_ccl4_cdki <- setNames(sce_ccl4_cdki$scDblFinder.class, colnames(sce_ccl4_cdki))

seurat_clean$scDblFinder.class <- dbl_class_control[colnames(seurat_clean)]
seurat_obj_cdki_clean$scDblFinder.class <- dbl_class_cdki[colnames(seurat_obj_cdki_clean)]
seurat_obj_ccl4_clean$scDblFinder.class <- dbl_class_ccl4[colnames(seurat_obj_ccl4_clean)]
seurat_obj_ccl4_cdki_clean$scDblFinder.class <- dbl_class_ccl4_cdki[colnames(seurat_obj_ccl4_cdki_clean)]  

table(seurat_clean$scDblFinder.class, useNA = "always")
table(seurat_obj_cdki_clean$scDblFinder.class, useNA = "always")
table(seurat_obj_ccl4_clean$scDblFinder.class, useNA = "always")
table(seurat_obj_ccl4_cdki_clean$scDblFinder.class, useNA = "always")

seurat_clean_control <- subset(seurat_clean, subset = scDblFinder.class == "singlet")
seurat_clean_cdki <- subset(seurat_obj_cdki_clean, subset = scDblFinder.class == "singlet")
seurat_clean_ccl4 <- subset(seurat_obj_ccl4_clean, subset = scDblFinder.class == "singlet")
seurat_clean_ccl4_cdki <- subset(seurat_obj_ccl4_cdki_clean, subset = scDblFinder.class == "singlet")

all_markers <- readRDS("2026_5_5_scRNA_seq/stst_all_markers.rds")

library(dplyr)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 25) %>%
  ungroup()
marker_list <- split(top_markers$gene, top_markers$cluster)

seurat_clean <- AddModuleScore(
  seurat_clean,
  features = marker_list,
  name = "CellTypeScore",
  nbin = 10,    # reduce from default 24
  ctrl = 25     # number of control genes per bin; lower if you still get errors
)

seurat_clean_cdki <- AddModuleScore(
  seurat_clean_cdki,
  features = marker_list,
  name = "CellTypeScore",
  nbin = 10,    # reduce from default 24
  ctrl = 25     # number of control genes per bin; lower if you still get errors
)

seurat_clean_ccl4 <- AddModuleScore(
  seurat_clean_ccl4,
  features = marker_list,
  name = "CellTypeScore",
  nbin = 10,    # reduce from default 24
  ctrl = 25     # number of control genes per bin; lower if you still get errors
)

seurat_clean_ccl4_cdki <- AddModuleScore(
  seurat_clean_ccl4_cdki,
  features = marker_list,
  name = "CellTypeScore",
  nbin = 10,    # reduce from default 24
  ctrl = 25     # number of control genes per bin; lower if you still get errors
)

score_cols <- paste0("CellTypeScore", seq_along(marker_list))
marker_list_control <- lapply(marker_list, function(genes) intersect(genes, rownames(seurat_clean)))
marker_list_cdki <- lapply(marker_list, function(genes) intersect(genes, rownames(seurat_clean_cdki)))
marker_list_ccl4 <- lapply(marker_list, function(genes) intersect(genes, rownames(seurat_clean_ccl4)))
marker_list_ccl4_cdki <- lapply(marker_list, function(genes) intersect(genes, rownames(seurat_clean_ccl4_cdki)))

score_cols <- paste0("CellTypeScore", seq_along(marker_list))
names(score_cols) <- names(marker_list)
                                                                                               
print(score_cols)

new_names <- setNames(names(score_cols), score_cols)  # reverse lookup: column -> celltype
colnames(seurat_clean@meta.data)[colnames(seurat_clean@meta.data) %in% score_cols] <- 
  new_names[colnames(seurat_clean@meta.data)[colnames(seurat_clean@meta.data) %in% score_cols]]
colnames(seurat_clean_cdki@meta.data)[colnames(seurat_clean_cdki@meta.data) %in% score_cols] <- 
  new_names[colnames(seurat_clean_cdki@meta.data)[colnames(seurat_clean_cdki@meta.data) %in% score_cols]]
colnames(seurat_clean_ccl4@meta.data)[colnames(seurat_clean_ccl4@meta.data) %in% score_cols] <- 
  new_names[colnames(seurat_clean_ccl4@meta.data)[colnames(seurat_clean_ccl4@meta.data) %in% score_cols]]
colnames(seurat_clean_ccl4_cdki@meta.data)[colnames(seurat_clean_ccl4_cdki@meta.data) %in% score_cols] <- 
  new_names[colnames(seurat_clean_ccl4_cdki@meta.data)[colnames(seurat_clean_ccl4_cdki@meta.data) %in% score_cols]]

DefaultAssay(seurat_clean)                                                                
DefaultAssay(seurat_clean_cdki)
DefaultAssay(seurat_clean_ccl4)                                                                
DefaultAssay(seurat_clean_ccl4_cdki)

seurat_clean <- NormalizeData(seurat_clean)                                                                
seurat_clean_cdki <- NormalizeData(seurat_clean_cdki)
seurat_clean_ccl4 <- NormalizeData(seurat_clean_ccl4)                                                                
seurat_clean_ccl4_cdki <- NormalizeData(seurat_clean_ccl4_cdki)

seurat_clean <- FindVariableFeatures(seurat_clean)
seurat_clean_cdki <- FindVariableFeatures(seurat_clean_cdki)
seurat_clean_ccl4 <- FindVariableFeatures(seurat_clean_ccl4)
seurat_clean_ccl4_cdki <- FindVariableFeatures(seurat_clean_ccl4_cdki)

length(VariableFeatures(seurat_clean))                                                                
length(VariableFeatures(seurat_clean_cdki))
length(VariableFeatures(seurat_clean_ccl4))                                                                
length(VariableFeatures(seurat_clean_ccl4_cdki))

seurat_clean <- ScaleData(seurat_clean)                                                                
seurat_clean_cdki <- ScaleData(seurat_clean_cdki)
seurat_clean_ccl4 <- ScaleData(seurat_clean_ccl4)                                                                
seurat_clean_ccl4_cdki <- ScaleData(seurat_clean_ccl4_cdki)

seurat_clean <- RunPCA(seurat_clean)                                                                
seurat_clean_cdki <- RunPCA(seurat_clean_cdki)
seurat_clean_ccl4 <- RunPCA(seurat_clean_ccl4)                                                                
seurat_clean_ccl4_cdki <- RunPCA(seurat_clean_ccl4_cdki)

Reductions(seurat_clean)                                                                 
Reductions(seurat_clean_cdki)
Reductions(seurat_clean_ccl4)                                                                
Reductions(seurat_clean_ccl4_cdki)

seurat_clean <- RunUMAP(seurat_clean, dims = 1:30)                                       
seurat_clean_cdki <- RunUMAP(seurat_clean_cdki, dims = 1:30)
seurat_clean_ccl4 <- RunUMAP(seurat_clean_ccl4, dims = 1:30)                                                                
seurat_clean_ccl4_cdki <- RunUMAP(seurat_clean_ccl4_cdki, dims = 1:30)

seurat_clean <- FindNeighbors(seurat_clean, dims = 1:30)
seurat_clean_cdki <- FindNeighbors(seurat_clean_cdki, dims = 1:30)
seurat_clean_ccl4 <- FindNeighbors(seurat_clean_ccl4, dims = 1:30)
seurat_clean_ccl4_cdki <- FindNeighbors(seurat_clean_ccl4_cdki, dims = 1:30)

seurat_clean <- FindClusters(seurat_clean, resolution = 0.5)
seurat_clean_cdki <- FindClusters(seurat_clean_cdki, resolution = 0.5)
seurat_clean_ccl4 <- FindClusters(seurat_clean_ccl4, resolution = 0.5)
seurat_clean_ccl4_cdki <- FindClusters(seurat_clean_ccl4_cdki, resolution = 0.5)

colnames(seurat_clean@meta.data)
colnames(seurat_clean_cdki@meta.data)
colnames(seurat_clean_ccl4@meta.data)
colnames(seurat_clean_ccl4_cdki@meta.data)

celltype_cols <- c("pDCs", "Kupffer cells", "T cells", "Monocytes & Monocyte-derived cells",
                    "ILC1s", "B cells", "cDC2s", "NK cells", "Neutrophils", "Mig. cDCs",
                    "cDC1s", "Endothelial cells", "Hepatocytes", "Cholangiocytes",
                    "Fibroblasts", "Basophils", "HsPCs")

cluster_scores_control <- seurat_clean@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(celltype_cols), mean)) %>%
  tibble::column_to_rownames("seurat_clusters")
cluster_assignment <- apply(cluster_scores_control, 1, function(x) names(x)[which.max(x)])
print(cluster_assignment)

cluster_scores_cdki <- seurat_clean_cdki@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(celltype_cols), mean)) %>%
  tibble::column_to_rownames("seurat_clusters")
cluster_assignment_cdki <- apply(cluster_scores_cdki, 1, function(x) names(x)[which.max(x)])

cluster_scores_ccl4 <- seurat_clean_ccl4@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(celltype_cols), mean)) %>%
  tibble::column_to_rownames("seurat_clusters")
cluster_assignment_ccl4 <- apply(cluster_scores_ccl4, 1, function(x) names(x)[which.max(x)])

cluster_scores_ccl4_cdki <- seurat_clean_ccl4_cdki@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(celltype_cols), mean)) %>%
  tibble::column_to_rownames("seurat_clusters")
cluster_assignment_ccl4_cdki <- apply(cluster_scores_ccl4_cdki, 1, function(x) names(x)[which.max(x)])

seurat_clean@meta.data$cell_type <- cluster_assignment[as.character(seurat_clean$seurat_clusters)]
seurat_clean_cdki@meta.data$cell_type <- cluster_assignment_cdki[as.character(seurat_clean_cdki$seurat_clusters)]
seurat_clean_ccl4@meta.data$cell_type <- cluster_assignment_ccl4[as.character(seurat_clean_ccl4$seurat_clusters)]
seurat_clean_ccl4_cdki@meta.data$cell_type <- cluster_assignment_ccl4_cdki[as.character(seurat_clean_ccl4_cdki$seurat_clusters)]

> table(seurat_clean_ccl4_cdki@meta.data$cell_type)

                       Hepatocytes Monocytes & Monocyte-derived cells 
                              2244                                 21 

> table(seurat_clean_cdki@meta.data$cell_type)

                       Hepatocytes Neutrophils 
       2867         193 
       
> table(seurat_clean_ccl4@meta.data$cell_type)

                    Cholangiocytes                  Endothelial cells 
                                34                                319 
                       Fibroblasts                        Hepatocytes 
                               587                              16886 
Monocytes & Monocyte-derived cells                        Neutrophils 
                               540                                423 
                          NK cells 
                               220 

> table(seurat_clean@meta.data$cell_type)

                             cDC1s                  Endothelial cells 
                               360                                187 
                       Fibroblasts                        Hepatocytes 
                               209                               3207 
Monocytes & Monocyte-derived cells 
                                31 


seurat_clean@meta.data$cell_type <- unname(cluster_assignment[as.character(seurat_clean$seurat_clusters)])
seurat_clean_cdki@meta.data$cell_type <- unname(cluster_assignment_cdki[as.character(seurat_clean_cdki$seurat_clusters)])
seurat_clean_ccl4@meta.data$cell_type <- unname(cluster_assignment_ccl4[as.character(seurat_clean_ccl4$seurat_clusters)])
seurat_clean_ccl4_cdki@meta.data$cell_type <- unname(cluster_assignment_ccl4_cdki[as.character(seurat_clean_ccl4_cdki$seurat_clusters)])

table(seurat_clean$cell_type, useNA = "always")

saveRDS(seurat_clean, file = "2026_7_17_seurat_obj/2026_7_20_Control_filtered.rds")
saveRDS(seurat_clean_cdki, file = "2026_7_17_seurat_obj/2026_7_20_CDKI_filtered.rds")
saveRDS(seurat_clean_ccl4, file = "2026_7_17_seurat_obj/2026_7_20_CCl4_filtered.rds")
saveRDS(seurat_clean_ccl4_cdki, file = "2026_7_17_seurat_obj/2026_7_20_CCl4_CDKI_filtered.rds")


## Preparing for integration 
library(Seurat)
library(patchwork)
dir <- getwd()
dir(file.path("2026_7_17_seurat_obj"))
         
sample_metadata <- data.frame(
  sample_id = c(
    "CCl4_CDKI", "CCl4", "CDKI", "Control"),
  stringsAsFactors = FALSE
)

sample_colors <- c(
  "#984EA3", "#E41A1C", "#377EB8", "#4DAF4A")  # Healthy 1-4: purple - ccl4_cdki, red - ccl4, blue - cdki, green - control

 
# Load QC-filtered Seurat objects 
qc_data_path <- "2026_7_17_seurat_obj"
 
seurat_list <- lapply(sample_metadata$sample_id, function(sample_id) {
  obj <- readRDS(file.path(qc_data_path, paste0(sample_id, "_filtered.rds")))
 
  # Ensure metadata is consistent
  obj$sample_id <- sample_id

  return(obj)
})
 
names(seurat_list) <- sample_metadata$sample_id
 
> all_genes <- lapply(seurat_list, rownames)
common_genes <- Reduce(intersect, all_genes)
 
# Subset all samples to common genes
seurat_list <- lapply(seurat_list, function(obj) {
  obj[common_genes, ]
})
 
# Report dimensions
cat("\nLoaded", length(seurat_list), "filtered samples:\n")
for (sample_id in names(seurat_list)) {
  cat(sprintf("  %s: %d cells × %d genes\n", 
              sample_id, 
              ncol(seurat_list[[sample_id]]), 
              nrow(seurat_list[[sample_id]])))
}
 
cat("\nLoaded", length(seurat_list), "filtered samples:\n")
for (sample_id in names(seurat_list)) {
  cat(sprintf("  %s: %d cells × %d genes\n", 
              sample_id, 
              ncol(seurat_list[[sample_id]]), 
              nrow(seurat_list[[sample_id]])))
}

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "Liver_CCL4"
)
cell_names <- colnames(merged_seurat)
merged_seurat$sample_id <- gsub("_[ACGT]{16}-1$", "", cell_names)

# 2. Verify the names are now completely distinct
unique(merged_seurat$sample_id)
# Expected: "CCl4_CDKI" "CCl4" "CDKI" "Control"

# Convert FILTERED to Assay5
merged_seurat[["FILTERED"]] <- as(merged_seurat[["FILTERED"]], "Assay5")

merged_seurat[["FILTERED"]] <- split(merged_seurat[["FILTERED"]], f = merged_seurat$sample_id)
Layers(merged_seurat[["FILTERED"]])   # should now show counts.CCl4_CDKI, counts.CCl4, counts.CDKI, counts.Control 

DefaultAssay(merged_seurat) <- "FILTERED"
Layers(merged_seurat, search = "counts")   # confirm 4 layers

merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)
merged_seurat <- FindVariableFeatures(merged_seurat, verbose = FALSE)
merged_seurat <- ScaleData(merged_seurat, verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, npcs = 50, verbose = FALSE)

integrated_cca <- IntegrateLayers(
  object = merged_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:30,
  verbose = FALSE
)
 
# Standard downstream workflow
integrated_cca <- FindNeighbors(integrated_cca, reduction = "integrated.cca", dims = 1:30, verbose = FALSE)
integrated_cca <- FindClusters(integrated_cca, resolution = 0.4, verbose = FALSE)
integrated_cca <- RunUMAP(integrated_cca, reduction = "integrated.cca", dims = 1:30, 
                          reduction.name = "umap.cca", verbose = FALSE)



p_celltype_cca <- DimPlot(integrated_cca,
                      reduction = "umap.cca",
                      group.by = "cell_type",
                      label = TRUE,
                      label.size = 3,
                      pt.size = 0.05) +
  theme(legend.position = "right",
        aspect.ratio = 1) +  # forces square plot area
  ggtitle("Cell Type")

p_condition_cca <- DimPlot(integrated_cca,
                       reduction = "umap.cca",
                       group.by = "orig.ident",
                       pt.size = 0.05) +
  theme(legend.position = "right",
        aspect.ratio = 1) +  # forces square plot area
  ggtitle("Condition")

# Combine side by side
combined <- p_celltype_cca | p_condition_cca

ggsave("2026_7_15_UMAP_Plots_cellbender_scdblfinder/2026_8_31_FINAL_integration_by_cca_umap_cell_clusters_condition.png", combined,
       width = 20, height = 8, dpi = 300)

colnames(integrated_cca@meta.data) <- gsub(".*_snn_res\\.", "clusters_res_", 
                                            colnames(integrated_cca@meta.data))

integrated_cca$seurat_clusters <- integrated_cca@meta.data[[paste0("clusters_res_", 0.4)]]
Idents(integrated_cca) <- "seurat_clusters"

# Calculate proportions for each sample
library(dplyr)

# Calculate proportions for each sample
proportion_data <- integrated_cca@meta.data %>%
  group_by(sample_id, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(
    total_cells = sum(count),
    proportion = count / total_cells,
    percentage = proportion * 100
  ) %>%
  ungroup()

 
# Save proportion data
write.csv(proportion_data, "results/cell_type_proportions_per_sample.csv", row.names = FALSE)
 
# Display summary
proportion_summary <- proportion_data %>%
  group_by(condition, cell_type) %>%
  summarise(
    mean_percentage = mean(percentage),
    sd_percentage = sd(percentage),
    .groups = "drop"
  )

  # Stacked barplot showing composition of each sample (percentages) -- note, this was the previous way of annotating (from 7_16), not the one done more recently (7_29) only because i have RDS from 7_16 already
p_stacked <- ggplot(proportion_data, 
                    aes(x = sample_id, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Cell Type Composition by Sample",
       x = "Sample", y = "Percentage (%)",
       fill = "Cell Type") +
  scale_fill_brewer(palette = "Set3")
 
ggsave("2026_7_15_UMAP_Plots_cellbender_scdblfinder/2026_7_31_01_stacked_composition.png", p_stacked,
       width = 12, height = 6, dpi = 300)

#2026_9_1_Downsampled_stacked_barplot
p_stacked <- ggplot(proportion_data, aes(x = sample_id, y = percentage, fill = cell_type)) + 
  geom_bar(stat = "identity", width = 0.6) + # Adjust this value (e.g., 0.4 to 0.7) to change thickness
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") + 
  labs(title = "Cell Type Composition by Sample", x = "Sample", y = "Percentage (%)", fill = "Cell Type") + 
  scale_fill_brewer(palette = "Set3")

ggsave(
  filename = "2026_7_15_UMAP_Plots_cellbender_scdblfinder/2026_9_1_Downsampled_stacked_barplot.pdf", 
  plot = p_stacked,
  device = cairo_pdf,
  width = 12, height = 6, dpi = 300)

# stacked barplot showing composition of each sample (counts)
p_stacked_counts <- ggplot(proportion_data, aes(x = sample_id, y = count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Cell Type Composition by Sample",
    x = "Sample",
    y = "Number of Cells",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

p_stacked_counts

ggsave("2026_7_15_UMAP_Plots_cellbender_scdblfinder/2026_7_31_02_stacked_composition_counts.png",
       p_stacked_counts, width = 8, height = 6, dpi = 300)
 
# Grouped barplot with error bars
p_grouped <- ggplot(proportion_summary, 
                    aes(x = cell_type, y = mean, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.9), width = 0.2) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Cell Type Proportions by Condition",
       x = "Cell Type", y = "Mean Percentage (%)") +
  scale_fill_manual(values = c("Healthy" = "#2E86AB", "Post_Treatment" = "#F18F01"))
 
ggsave("plots/proportions/02_grouped_proportions.png", p_grouped,
       width = 12, height = 6, dpi = 300)





#-----------------------------------------------
# STEP 6: Cell-level DE analysis with FindMarkers
#-----------------------------------------------
 
# Get unique cell types
cell_types <- unique(integrated_cca$cell_type)
integrated_cca <- JoinLayers(integrated_cca)
options(future.globals.maxSize = 4 * 1024^3)
 
 
# Function to run FindMarkers for one cell type
run_findmarkers <- function(celltype, integrated_cca, ident1, ident2) {
  cat("  Processing:", celltype, "-", ident1, "vs", ident2, "\n")
  
  seurat_obj_sub <- subset(integrated_cca, subset = cell_type == celltype)
  
  markers <- FindMarkers(
    seurat_obj_sub,
    ident.1 = ident1,
    ident.2 = ident2,
    group.by = "sample_id",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1,
    verbose = FALSE
  )
  
  markers$gene <- rownames(markers)
  markers$cell_type <- celltype
  markers$comparison <- paste0(ident1, "_vs_", ident2)
  
  return(markers)
}

comparisons <- list(
  c("CCl4", "Control"),
  c("CDKI", "Control"),
  c("CCl4_CDKI", "Control")
)

all_results <- lapply(comparisons, function(comp) {
  res_by_celltype <- lapply(
    cell_types, run_findmarkers,
    integrated_cca = integrated_cca,
    ident1 = comp[1],
    ident2 = comp[2]
  )
  names(res_by_celltype) <- cell_types
  res_by_celltype
})
names(all_results) <- sapply(comparisons, function(c) paste0(c[1], "_vs_", c[2]))
 


##### JUST HEPATOCYTES ########
hepatocytes <- subset(integrated_cca, subset = cell_type == "Hepatocytes")

# Sanity check cell counts per condition before testing
table(hepatocytes$sample_id)

comparisons <- list(
  c("CCl4", "Control"),
  c("CDKI", "Control"),
  c("CCl4_CDKI", "Control")
)

hepatocyte_results <- lapply(comparisons, function(comp) {
  cat("Comparing:", comp[1], "vs", comp[2], "\n")
  
  markers <- FindMarkers(
    hepatocytes,
    ident.1 = comp[1],
    ident.2 = comp[2],
    group.by = "sample_id",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1,
    verbose = FALSE
  )
  
  markers$gene <- rownames(markers)
  markers$comparison <- paste0(comp[1], "_vs_", comp[2])
  
  return(markers)
})

names(hepatocyte_results) <- sapply(comparisons, function(c) paste0(c[1], "_vs_", c[2]))

# Add significance flag
findmarkers_all$significant <- findmarkers_all$p_val_adj < 0.05 & 
                               abs(findmarkers_all$avg_log2FC) > 0.5
 
# Save results
write.csv(findmarkers_all, "results/findmarkers_celllevel_results.csv", row.names = FALSE)

# ORIGINAL CODE:
hepatocyte_results_df <- bind_rows(hepatocyte_results, .id = "comparison")
hepatocyte_results_df <- hepatocyte_results_df %>%
  group_by(comparison) %>%
  mutate(q_value = p.adjust(p_val, method = "BH")) %>%
  ungroup()
> hepatocyte_results_df$significant <- hepatocyte_results_df$q_value < 0.05 & 
                                      abs(hepatocyte_results_df$avg_log2FC) > 0.25

> sig_summary <- hepatocyte_results_df %>%
  group_by(comparison) %>%
  summarise(
    total_genes_tested = n(),
    sig_genes = sum(significant),
    sig_up = sum(avg_log2FC > 0.25 & q_value < 0.05 & significant),
    sig_down = sum(avg_log2FC < -0.25 & q_value < 0.05 & significant),
    .groups = "drop"
  )

dir.create("2026_8_4_diff_expression_results", showWarnings = FALSE)  # in case it doesn't exist yet
write.csv(hepatocyte_results_df, "2026_8_4_diff_expression_results/2026_8_4_hepatocyte_findmarkers_results.csv", row.names = FALSE)
write.csv(sig_summary, "2026_8_4_diff_expression_results/2026_8_4_hepatocyte_sig_summary.csv", row.names = FALSE)
