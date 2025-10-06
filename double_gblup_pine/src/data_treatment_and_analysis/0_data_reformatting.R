# script meant to reformat data for genomic prediction
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(devtools)
install_other_requirements <- F
if (install_other_requirements) {
  install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install("snpStats")
  BiocManager::install("mixOmicsTeam/mixOmics")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
}
library(mixOmics)
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(Matrix)
library(graphics)
library(htmlwidgets)
library(rstudioapi)
library(stringr)
library(tidyr)
library(dplyr)
library(lsmeans)

# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../functions.R")

# set options to increase memory
options(expressions = 5e5)
options(warn = -1)
emm_options(rg.limit = 10e6)

# set path for genomic data and phenotype data
genom_dir_path <- "../../data/genomic_data/"
pheno_dir_path <- "../../data/phenotype_data/"

# set output result path for genomic graphics
output_genom_graphics_path <- "../../results/genomic_prediction_graphics/"

# set maximum number of principal components to be tested using akaike
# information criterion
max_n_comp_ <- 10

# umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 15
min_dist_ <- 0.1

# define traits_
traits_ <- c("H", "I", "D", "T4", "T5", "T6")

# get phenotypic data
pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "PhenotypeField.txt"
)))
# drop last empty column
pheno_df <- pheno_df[, -ncol(pheno_df)]
colnames(pheno_df)[match("ID", colnames(pheno_df))] <- "Genotype"

# transform data with pivot_longer to
# group columns H, I, D, T4, T5, T6 based on years
pheno_df <- as.data.frame(pheno_df %>%
  pivot_longer(
    cols = starts_with(c("H", "I", "D", "T4", "T5", "T6")),
    names_to = c(".value", "Year"),
    names_pattern = "([A-Za-z0-9]+)[-]?(\\d{2})"
  ) %>%
  mutate(Year = as.integer(paste0("20", Year))))

# create environment variable
pheno_df$Envir <- paste0(
  pheno_df$Site, "_",
  pheno_df$Year, "_",
  pheno_df$Block
)

# generate latitude and longitude variables per environment
pheno_df <- generate_latitude_longitude_variables_by_environment(
  pheno_df
)

# get genomic data
geno_df <- t(as.data.frame(fread(paste0(
  genom_dir_path,
  "genotypefield.txt"
))))
colnames_ <- geno_df[1, ]
colnames(geno_df) <- colnames_
geno_df <- geno_df[-1, ]
geno_df <- replace_nn_with_most_frequent(geno_df)
rownames_ <- rownames(geno_df)
geno_df <- recode_genotypes(geno_df) # recode as 0, 1 and 2
rownames(geno_df) <- rownames_

# write reformatted datasets
fwrite(pheno_df,
  file = paste0(pheno_dir_path, "phenotype_data.csv")
)
fwrite(as.data.frame(geno_df),
  file = paste0(genom_dir_path, "genomic_data.csv"),
  row.names = T
)

# convert geno_df to numeric matrix
geno_mat <- scale(apply(geno_df, 2, as.numeric), center = T, scale = F)
k_mat <- crossprod(t(geno_mat))

# create a heatmap for genomic covariance matrix
png(
  paste0(
    output_genom_graphics_path,
    "pine_genomic_covariance_heatmap.png"
  ),
  width = 1200, height = 1200, res = 150
)
pheatmap(k_mat,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Pine dataset genomic covariance matrix heatmap",
  color = colorRampPalette(c("blue", "red"))(100),
  show_rownames = T,
  show_colnames = T
)
dev.off()

# compute pca for genomic data
geno_pca <- mixOmics::pca(apply(geno_df, 2, as.numeric), 2,
  center = T,
  scale = F,
  ncomp = max_n_comp_
)
pc_coord_df_ <- as.data.frame(geno_pca$variates)[, 1:max_n_comp_]
colnames(pc_coord_df_) <- str_replace_all(
  colnames(pc_coord_df_),
  pattern = "X.", replacement = ""
)
pc_var_names_ <- colnames(pc_coord_df_)
pc_coord_df_$Genotype <- rownames(geno_df)
geno_pca_exp_var_ <- geno_pca$explained_variance

# create pca plot
pca_fig_x_y <- plot_ly(
  data = pc_coord_df_,
  x = ~PC1, y = ~PC2,
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for pine dataset genomic data (20795 SNP)",
    xaxis = list(title = paste0(
      names(geno_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(geno_pca_exp_var_)[1], 2), "%"
    )),
    yaxis = list(title = paste0(
      names(geno_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(geno_pca_exp_var_)[2], 2), "%"
    ))
  )

# save pca graphic
saveWidget(pca_fig_x_y, file = paste0(
  output_genom_graphics_path, "pine_pca_genomic_data.html"
))

# compute umap in 2d
geno_umap_2d_model <- umap(
  apply(geno_df, 2, as.numeric),
  n_components = 2,
  random_state = random_state_umap_,
  n_neighbors = n_neighbors_umap_,
  min_dist = min_dist_
)
geno_umap_2d <- data.frame(geno_umap_2d_model[["layout"]])

# create umap plot
umap_fig_x_y <- plot_ly(
  data = geno_umap_2d,
  x = ~X1, y = ~X2,
  type = "scatter", mode = "markers", color = "orange"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "UMAP 2D plot for pine dataset genomic data (20795 SNP)",
    xaxis = list(title = "First component"),
    yaxis = list(title = "Second component")
  )

# save umap graphic
saveWidget(umap_fig_x_y, file = paste0(
  output_genom_graphics_path, "pine_umap_genomic_data.html"
))
