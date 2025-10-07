# script meant to compute gblup estimated breeding values (ebvs) for all traits
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
library(devtools)
if ("pine_env" %in% conda_list()$name) {
  use_condaenv("pine_env")
}
library(tidyverse)
library(tidyr)
library(data.table)
library(lubridate)
library(plotly)
library(htmlwidgets)
library(emmeans)
library(SpATS)
library(stringr)
library(lme4)
library(anytime)
library(foreach)
library(parallel)
library(doParallel)
library(MASS)
library(lme4)
library(vcd)
library(car)
# install_github('ljacquin/KRMM')
library(KRMM)

# define computation mode, i.e. "local" or "cluster"
computation_mode <- "local"

# if comutations are local in rstudio, detect and set script path
# automatically using rstudioapi
if (identical(computation_mode, "local")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}

# source functions
source("../functions.R")

# set options to increase memory and suppress warnings
options(expressions = 1e5)
emm_options(rg.limit = 5e5)
options(warn = -1)

# set paths
pheno_dir_path_ <- "../../data/phenotype_data/"
pheno_file_path_ <- paste0(
  pheno_dir_path_,
  "phenotype_data.csv"
)
spats_adj_pheno_path <- paste0(
  pheno_dir_path_,
  "spats_per_env_adjusted_phenotypes/"
)
genom_dir_path <- "../../data/genomic_data/"

# set vars to keep for gblup computation
vars_to_keep_ <- c("Envir", "Genotype")

# get file names for spats adjusted phenotypes and replace pattern
# "_spats_adjusted_.*" with "" to get trait names
files_names_spats_adj_pheno <- list.files(spats_adj_pheno_path)
trait_names_ <- str_replace_all(files_names_spats_adj_pheno,
  "_spats_adjusted_.*",
  replacement = ""
)
# initialize lists for gblup associated to all traits
list_gblup_ebv_per_geno <- vector("list", length(trait_names_))
names(list_gblup_ebv_per_geno) <- trait_names_

# get genomic data
geno_df <- as.data.frame(fread(paste0(genom_dir_path, "genomic_data.csv")))

# for each trait, compute gblups across all environments
for (file_ in files_names_spats_adj_pheno) {
  print(paste0("computation for file : ", file_))

  df_ <- as.data.frame(fread(paste0(spats_adj_pheno_path, file_)))
  df_ <- df_[, c(vars_to_keep_, colnames(df_)[str_detect(
    colnames(df_),
    "spats_adj_pheno"
  )])]
  trait_ <- colnames(df_)[str_detect(colnames(df_), "spats_adj_pheno")]
  # for isolated genotypes with single phenotypic values (i.e. non repeated),
  # repeat their unique phenotypic values one time to make blup computation
  # possible
  idx_non_repeat_geno <- which(table(df_$Genotype) < 2)
  if (length(idx_non_repeat_geno) > 0) {
    df_ <- rbind(df_, df_[idx_non_repeat_geno, ])
  }

  # remove any NA for trait adjusted phenotype before gblup
  df_ <- df_[!is.na(df_[, trait_]), ]

  if (length(unique(df_$Envir)) > 1) {
    # compute gblup estimated breeding values (ebv) to be used as phenotypes

    # make incidence for fixed effects
    x_mat_fixed_eff <- model.matrix(~ Envir - 1, data = df_)
    x_mat_fixed_eff <- cbind(
      rep(1, nrow(x_mat_fixed_eff)), x_mat_fixed_eff
    )
    colnames(x_mat_fixed_eff)[1] <- "Intercept"
    colnames(x_mat_fixed_eff) <- str_remove_all(
      colnames(x_mat_fixed_eff),
      pattern = "Envir"
    )
    x_mat_fixed_eff <- x_mat_fixed_eff[, c("Intercept", unique(df_$Envir))]

    # make incidence for random effects
    z_mat_rand_eff <- model.matrix(~ Genotype - 1, data = df_)
    colnames(z_mat_rand_eff) <- str_remove_all(
      colnames(z_mat_rand_eff),
      pattern = "Genotype"
    )
    z_mat_rand_eff <- z_mat_rand_eff[, unique(df_$Genotype)]

    # order marker data with respect to random effects
    geno_mat <- geno_df[match(unique(df_$Genotype), geno_df$V1), ]

    # compute gblup ebv
    linear_krmm_model <- krmm(
      Y = df_[, trait_], X = x_mat_fixed_eff, Z = z_mat_rand_eff,
      Matrix_covariates = geno_mat[, -1], method = "GBLUP"
    )

    gblup_ebv <- predict_krmm(linear_krmm_model,
      Matrix_covariates = geno_mat[, -1],
      X = x_mat_fixed_eff,
      Z = z_mat_rand_eff,
      add_fixed_effects = F
    )
    gblup_ebv_df <- data.frame("Genotype" = geno_mat$V1, "Trait" = gblup_ebv)
    colnames(gblup_ebv_df)[2] <- str_remove(trait_, "_spats_adj_pheno")

    # add data frame to list
    list_gblup_ebv_per_geno[[
      str_remove(trait_, "_spats_adj_pheno")
    ]] <- gblup_ebv_df
  }
}

# reduce gblup ebv list
gblup_ebv_df <- Reduce(
  function(x, y) {
    merge(x, y, by = "Genotype", all = T)
  },
  list_gblup_ebv_per_geno
)

# write gblup ebvs
colnames(gblup_ebv_df) <- c("Genotype", trait_names_)

fwrite(gblup_ebv_df,
  file = paste0(
    pheno_dir_path_,
    "gblup_ebv.csv"
  )
)
