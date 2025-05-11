# 📊 Sketching a Seurat Object Using Leverage Score Sampling

This script performs sketching (downsampling) of a Seurat object using **Leverage Score Sampling**, followed by PCA analysis and plotting. It helps reduce the number of cells in your dataset while retaining biological diversity and data structure.

---

## 📁 File Structure

- `sketch.R` – The main R script to perform sketching and generate visualizations.
- `figures/` – Output directory that will store the generated plots.
- `README.md` – This documentation file.

---

## 🚀 How to Use

### 1. Prerequisites

Install the required R packages:

```r
install.packages(c("ggplot2", "dplyr"))
# For Seurat:
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
