# script meant to perform genomic prediction and analyses for pine
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
if ("pine_env" %in% conda_list()$name) {
  print("using pine_env")
  use_condaenv("pine_env")
}
# install other requirements from github if necessary
install_other_requirements <- F
if (install_other_requirements) {
  # reticulate::install_miniconda()
  conda_create("pine_env")
  use_condaenv("pine_env")
  library(devtools)
  devtools::install_github("ljacquin/KRMM")
  devtools::install_github("rstudio/tensorflow")
  library(tensorflow)
  install_tensorflow(envname = "pine_env")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
  install.packages("umap")
}
use_tensorflow_or_umap <- F
if (use_tensorflow_or_umap) {
  # leave tensorflow and keras for later use
  library(tensorflow)
  library(keras3)
  library(umap)
  tensorflow::tf$random$set_seed(0)
  py_module_available("keras") # must return TRUE
  py_module_available("tensorflow") # must return TRUE
  py_discover_config("keras") # more info on the python env, tf and keras
}
library(MASS)
library(data.table)
library(stringr)
library(lme4)
library(tidyr)
library(doParallel)
library(doRNG)
library(robustbase)
library(foreach)
library(parallel)
library(Matrix)
library(matrixcalc)
library(rgl)
library(Rfast)
library(cvTools)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(KRMM)
library(kernlab)
library(future)
library(future.apply)
library(grDevices)
library(matrixcalc)
library(Matrix)

# define computation mode, i.e. local or cluster
computation_mode <- "local"

# if comutations are local in rstudio, detect and set script path
# automatically using rstudioapi
if (identical(computation_mode, "local")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}

# source functions
source("../functions.R")

# set options
options(future.globals.maxSize = 60 * 1024^3)
options(expressions = 5e5)
options(warn = -1)

# set color gradients and color vector for predictive abilities (pa)
pa_colors_ <- "red"

# set color vector for computed genomic heritabilities (h2)
h2_colors_ <- "brown"

# define number of cores
nb_cores_ <- 12

# define function(s) and package(s) to export for parallelization
pkgs_to_export_ <- c(
  "kernlab",
  "KRMM",
  "foreach",
  "cvTools",
  "Matrix",
  "MASS"
)
# set input paths
geno_dir_path <- "../../data/genomic_data/"
pheno_dir_path <- "../../data/phenotype_data/"

# output result path for genotype graphics
output_pred_results_path <- "../../results/genomic_prediction/"
output_pred_graphics_path <- "../../results/genomic_prediction_graphics/"

# define kernel
kernel_ <- "linear"

# define traits for genomic prediction and analyses
traits_ <- c("H", "I", "D", "T4", "T5", "T6")

# define shift seed value by
mult_seed_by_ <- 100

# set k for K-folds cv
k_folds_ <- 5

# define number of shuffles
n_shuff_ <- 20

for (trait_ in traits_)
{
  print(trait_)

  # get ebvs
  ebv_df <- as.data.frame(fread(paste0(
    pheno_dir_path,
    "gblup_ebv.csv"
  )))[, c("Genotype", trait_)]

  # get genotype data
  omic_df <- as.data.frame(fread(paste0(
    geno_dir_path,
    "genomic_data.csv"
  )))
  colnames(omic_df)[1] <- "Genotype"

  # remove rows with na, inf or char associated to trait for ebv_df
  idx_na_inf_char <- which(is.na(ebv_df[, trait_]) |
    is.infinite(ebv_df[, trait_]) |
    is.na(suppressWarnings(as.numeric(ebv_df[, trait_]))))
  if (length(idx_na_inf_char) > 0) {
    ebv_df <- ebv_df[-idx_na_inf_char, ]
  }

  # merge ebv_df and omic_df for integrity of analyses
  merged_df <- merge(ebv_df, omic_df, by = "Genotype")
  ebv_df <- merged_df[, c("Genotype", trait_)]
  omic_df <- merged_df[, -match(c("Genotype", trait_), colnames(merged_df))]
  rownames(omic_df) <- merged_df$Genotype
  rm(merged_df)

  # get trait ebvs
  trait_ebv <- ebv_df[, trait_]

  # register parallel backend
  cl <- makeCluster(nb_cores_)
  registerDoParallel(cl)

  # get number of genotypes
  n <- nrow(omic_df)

  # create folds for k-fold cross-validation
  Folds <- cvFolds(n, K = k_folds_, type = "consecutive")

  df_result_ <- foreach(
    shuff_ = 1:n_shuff_,
    .packages = pkgs_to_export_,
    .combine = rbind
  ) %dopar% {
    # set seed, define a new set of indices,
    set.seed(shuff_ * mult_seed_by_)
    idx_ <- sample(1:n, size = n, replace = F)
    fold_results <- foreach(
      fold_ = 1:k_folds_,
      .packages = pkgs_to_export_,
      .combine = rbind
    ) %dopar% {
      # define indices for validation and training based on k-folds cv
      idx_val_fold <- which(Folds$which == fold_)
      idx_val <- idx_[idx_val_fold]
      idx_train <- idx_[-idx_val_fold]

      # initialize vector of results for the current fold
      fold_result <- c(
        "double_GBLUP_ebv_pa" = NA,
        "double_GBLUP_ebv_h2" = NA
      )

      # training and prediction based on ebvs

      # train and predict with GBLUP (linear kernel krmm)
      linear_krmm_model <- krmm(
        Y = trait_ebv[idx_train],
        X = rep(1, length(idx_train)),
        Matrix_covariates = omic_df[idx_train, ],
        method = "GBLUP",
      )
      f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
        Matrix_covariates = omic_df[idx_val, ],
        add_fixed_effects = F
      )
      fold_result["double_GBLUP_ebv_pa"] <- cor(
        f_hat_val_linear_krmm,
        trait_ebv[idx_val]
      )
      fold_result["double_GBLUP_ebv_h2"] <- compute_genomic_h2(
        linear_krmm_model$sigma2K_hat,
        linear_krmm_model$sigma2E_hat
      )
      fold_result
    }
    fold_results
  }
  # stop the parallel backend
  stopCluster(cl)
  registerDoSEQ()

  # create directory for trait_ graphics if it does not exist
  if (!dir.exists(paste0(output_pred_graphics_path, trait_, "/"))) {
    dir.create(paste0(output_pred_graphics_path, trait_, "/"))
  }

  # get methods pa
  df_pa_ <- as.data.frame(df_result_[, 1])
  colnames(df_pa_) <- "double_GBLUP_ebv_pa"
  method_names <- colnames(df_pa_)
  df_pa_ <- as.data.frame(apply(df_pa_, 2, as.numeric))

  # initialize plot_ly for violin + boxplots
  violin_box_plots_pa_ <- plot_ly()

  # add violin plots
  for (i in seq_along(method_names)) {
    violin_box_plots_pa_ <- add_trace(
      violin_box_plots_pa_,
      type = "violin", # specify violin plot
      y = df_pa_[[method_names[i]]],
      name = method_names[i],
      points = F, # show all points
      jitter = 0.3, # add jitter for spread
      pointpos = -1.8, # adjust point position relative to the violin plot
      marker = list(color = pa_colors_[i]),
      fillcolor = pa_colors_[i],
      line = list(color = pa_colors_[i]),
      meanline = list(visible = F), # do not show mean line
      scalemode = "width", # keep width constant for comparison
      opacity = 0.6 # slight transparency to see the boxplot behind it
    )
  }

  # add boxplots on top of the violin plots but hide from legend
  for (i in seq_along(method_names)) {
    violin_box_plots_pa_ <- add_boxplot(
      violin_box_plots_pa_,
      y = df_pa_[[method_names[i]]],
      name = method_names[i],
      marker = list(color = pa_colors_[i]),
      line = list(color = "black", width = 2), # black line for the boxplot
      fillcolor = "rgba(255,255,255,0)", # transparent fill to see the violin plot underneath
      width = 0.2, # narrower boxplot to fit inside the violin
      notchwidth = 0.4, # add notch for median
      showlegend = F # hide boxplot from legend
    )
  }

  # add layout
  violin_box_plots_pa_ <- violin_box_plots_pa_ %>%
    layout(
      title = paste0(
        "Genomic PA distributions of models for ",
        trait_, ", based on ", ncol(omic_df), " SNP across ",
        n_shuff_, " shuffling scenarios for ", k_folds_, "-folds cv"
      ),
      yaxis = list(
        title = "Predictive ability (PA)",
        range = c(-0.8, 1)
      ),
      legend = list(title = list(text = "Prediction model"))
    )

  # save violin_box_plots_pa_ graphics
  saveWidget(violin_box_plots_pa_, file = paste0(
    output_pred_graphics_path, trait_, "/genomic_pred_results_",
    trait_, "_", kernel_, "_kernel_", ncol(omic_df), "_SNP_",
    k_folds_, "_folds_CV.html"
  ))

  # add stats and save predictive ability results
  df_pa_ <- signif(apply(df_pa_, 2, as.numeric), 2)
  rownames(df_pa_) <- paste0("gp_scenario_", 1:nrow(df_pa_))

  df_stat <- as.data.frame(rbind(
    apply(df_pa_, 2, function(x) median(x, na.rm = TRUE)),
    apply(df_pa_, 2, function(x) IQR(x, na.rm = TRUE))
  ))
  df_stat <- signif(apply(df_stat, 2, as.numeric), 2)
  rownames(df_stat) <- c("pa_median", "pa_iqr")
  df_stat <- as.data.frame(df_stat)

  df_pa_ <- rbind(df_pa_, df_stat)

  fwrite(df_pa_,
    file = paste0(
      output_pred_results_path,
      "genomic_pred_results_", ncol(omic_df), "_SNP_",
      trait_, "_", kernel_, "_kernel_", k_folds_, "_folds_CV.csv"
    ), row.names = T
  )
  print(df_pa_)

  # get computed h2 for gblup
  df_h2_ <- as.data.frame(df_result_[, 2])
  colnames(df_h2_) <- "double_GBLUP_ebv_h2"
  method_names <- colnames(df_h2_)
  df_h2_ <- as.data.frame(apply(df_h2_, 2, as.numeric))

  # initialize plot_ly for violin + boxplots
  violin_box_plots_h2_ <- plot_ly()

  # add violin plots without points
  for (i in seq_along(method_names)) {
    violin_box_plots_h2_ <- add_trace(
      violin_box_plots_h2_,
      type = "violin", # specify violin plot
      y = df_h2_[[method_names[i]]],
      name = method_names[i],
      points = F, # remove points
      jitter = 0.3, # add jitter for spread (not used since points are removed)
      pointpos = -1.8, # adjust point position relative to the violin plot (not used here)
      marker = list(color = h2_colors_[i]),
      fillcolor = h2_colors_[i],
      line = list(color = h2_colors_[i]),
      meanline = list(visible = F), # do not show mean line
      scalemode = "width", # keep width constant for comparison
      opacity = 0.6 # slight transparency to see the boxplot behind it
    )
  }

  # add boxplots on top of the violin plots but hide from legend
  for (i in seq_along(method_names)) {
    violin_box_plots_h2_ <- add_boxplot(
      violin_box_plots_h2_,
      y = df_h2_[[method_names[i]]],
      name = method_names[i],
      marker = list(color = h2_colors_[i]),
      line = list(color = "black", width = 2), # black line for the boxplot
      fillcolor = "rgba(255,255,255,0)", # transparent fill to see the violin plot underneath
      width = 0.2, # narrower boxplot to fit inside the violin
      notchwidth = 0.4, # add notch for median
      showlegend = F # hide boxplot from legend
    )
  }

  # add layout
  violin_box_plots_h2_ <- violin_box_plots_h2_ %>%
    layout(
      title = paste0(
        "Genomic heritability (h2) distributions estimated from GBLUP prediction models for ",
        trait_, ", based on ", ncol(omic_df), " SNP across ",
        n_shuff_, " shuffling scenarios for ", k_folds_, "-folds cv"
      ),
      yaxis = list(
        title = "Genomic heritability (h2)",
        range = c(0, 1)
      ),
      legend = list(title = list(text = "Prediction model"))
    )

  # save violin_box_plots_h2_ graphics
  saveWidget(violin_box_plots_h2_, file = paste0(
    output_pred_graphics_path, trait_, "/genomic_h2_results_",
    trait_, "_", kernel_, "_kernel_", ncol(omic_df), "_SNP_",
    k_folds_, "_folds_CV.html"
  ))

  # add stats and save h2 results
  df_h2_ <- signif(apply(df_h2_, 2, as.numeric), 2)
  rownames(df_h2_) <- paste0("gp_scenario_", 1:nrow(df_h2_))

  df_stat <- as.data.frame(rbind(
    apply(df_h2_, 2, function(x) median(x, na.rm = TRUE)),
    apply(df_h2_, 2, function(x) IQR(x, na.rm = TRUE))
  ))
  df_stat <- signif(apply(df_stat, 2, as.numeric), 2)
  rownames(df_stat) <- c("h2_median", "h2_iqr")
  df_stat <- as.data.frame(df_stat)

  df_h2_ <- rbind(df_h2_, df_stat)

  fwrite(df_h2_,
    file = paste0(
      output_pred_results_path,
      "genomic_h2_results_", ncol(omic_df), "_SNP_",
      trait_, "_", kernel_, "_kernel_",
      k_folds_, "_folds_CV.csv"
    ), row.names = T
  )
  print(df_h2_)
}
