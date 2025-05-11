# sketch.R
# --------------------------------------------
# Title: Sketching Analysis for Seurat Object
# Author: Ibrahim
# Description: This script provides a function `run_sketch_analysis()` 
#              that takes in a Seurat object and runs SketchData-based 
#              downsampling followed by PCA, UMAP embedding, and generates 
#              diagnostic plots, as part of the main tutorials of IEEE sessions.
# 
# Input:
#   - seurat_obj: A Seurat object.
#   - output_dir: Directory path where output figures will be saved.
#
# Output:
#   - sk_seurat: The sketched Seurat object.
#   - Plots:
#       - UMAP (before and after)
#       - PCA (before and after)
#       - Elbow plots (before and after)
#       - Variance explained comparison plot
#       - Patient-wise cell proportion comparison plot
# --------------------------------------------

run_sketch_analysis <- function(seurat_obj, output_dir = "./figures", sketch_cells = 50000) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  library(Seurat)
  library(ggplot2)
  library(dplyr)

  # Before sketching: PCA and UMAP plots
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  ggsave(file.path(output_dir, "before_UMAP_plot_celltype_major.png"),
         plot = DimPlot(seurat_obj, group.by ="celltype_major"), width = 6, height = 4, dpi = 300)
  
  ggsave(file.path(output_dir, "before_PCA_plot_celltype_major.png"), 
         plot = DimPlot(seurat_obj, reduction = "pca", group.by = "celltype_major"), width = 6, height = 4, dpi = 300)

  ggsave(file.path(output_dir, "before_PCA_plot_Patients.png"), 
         plot = DimPlot(seurat_obj, reduction = "pca"), width = 6, height = 4, dpi = 300)

  ggsave(file.path(output_dir, "before_elbow_plot.png"), 
         plot = ElbowPlot(seurat_obj), width = 6, height = 4, dpi = 300)

  # Perform sketching
  sketched <- SketchData(seurat_obj, ncells = sketch_cells, method = "LeverageScore",
                         sketched.assay = "sketch", over.write = TRUE)

  sk_seurat <- seurat_obj[, colnames(sketched@assays$sketch)]
  sk_seurat <- RunPCA(sk_seurat)
  ggsave(file.path(output_dir,"after_UMAP_celltype_major.png"),
         plot = DimPlot(sk_seurat, group.by ="celltype_major"), width = 6, height = 4, dpi = 300)
         
  ggsave(file.path(output_dir, "after_PCA_plot_celltype_major.png"), 
         plot = DimPlot(sk_seurat, reduction = "pca", group.by = "celltype_major"), width = 6, height = 4, dpi = 300)

  ggsave(file.path(output_dir, "after_PCA_plot_Patients.png"), 
         plot = DimPlot(sk_seurat, reduction = "pca"), width = 6, height = 4, dpi = 300)

  ggsave(file.path(output_dir, "after_elbow_plot.png"), 
         plot = ElbowPlot(sk_seurat), width = 6, height = 4, dpi = 300)

  # Variance explained comparison
  pca_before <- seurat_obj[["pca"]]
  var_exp_before <- pca_before@stdev^2
  var_exp_before_ratio <- var_exp_before / sum(var_exp_before)

  pca_after <- sk_seurat[["pca"]]
  var_exp_after <- pca_after@stdev^2
  var_exp_after_ratio <- var_exp_after / sum(var_exp_after)

  png(file.path(output_dir, "variance_explained_comparison.png"), width = 800, height = 600)
  plot(cumsum(var_exp_before_ratio), type = "l", col = "blue", ylim = c(0,1),
       ylab = "Cumulative Variance Explained", xlab = "PCs",
       main = "Variance Explained Before vs After Sketching")
  lines(cumsum(var_exp_after_ratio), col = "red")
  legend("bottomright", legend = c("Before", "After"), col = c("blue", "red"), lty = 1)
  dev.off()

  # Print cumulative variance for first 20 PCs
  message("Cumulative variance (first 20 PCs): ", 
          round(sum(var_exp_after_ratio[1:20]) / sum(var_exp_before_ratio[1:20]), 3))

  # Cell proportion comparison
  before_cell_props <- seurat_obj@meta.data %>%
    group_by(Patient) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(prop = n / sum(n) * 100)

  after_cell_props <- sk_seurat@meta.data %>%
    group_by(Patient) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(prop = n / sum(n) * 100)

  comparison_df <- merge(before_cell_props, after_cell_props, by = "Patient", suffixes = c("_before", "_after"))

  p <- ggplot(comparison_df, aes(x = Patient)) +
    geom_bar(aes(y = prop_before, fill = "Before"), stat = "identity", position = "dodge") +
    geom_bar(aes(y = prop_after, fill = "After"), stat = "identity", position = "dodge") +
    labs(y = "Percentage of Cells", x = "Patient", fill = "Condition") +
    theme_minimal() +
    scale_fill_manual(values = c("Before" = "blue", "After" = "red")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  ggsave(file.path(output_dir, "cell_proportion_comparison.png"), plot = p, width = 10, height = 6, dpi = 300)

  return(sk_seurat)
}
