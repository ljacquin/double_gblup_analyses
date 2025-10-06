# script meant to compute gblup estimated breeding values (ebvs)
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
library(devtools)
if ("rice_env" %in% conda_list()$name) {
  use_condaenv("rice_env")
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
genom_dir_path <- "../../data/genomic_data/"

# define traits
traits_ <- c("FL", "PH", "YLD", "ZN")

# set vars to keep for gblup computation
vars_to_keep_ <- c("Envir", "Genotype", traits_)

# read genomic data
geno_df <- as.data.frame(fread(paste0(genom_dir_path, "genomic_data.csv")))
colnames(geno_df)[1] <- "Genotype"
geno_df$Genotype <- as.factor(geno_df$Genotype)

# for each trait, compute gblups across all environments
# for each genotype
pheno_df <- as.data.frame(fread(paste0(pheno_dir_path_, "phenotype_data.csv")))
pheno_df <- pheno_df[, vars_to_keep_]

# initialize lists for gblup associated to all traits
list_gblup_ebv <- vector("list", length(traits_))
names(list_gblup_ebv) <- traits_

for (trait_ in traits_) {
  
  print(trait_)
  
  # remove any NA for trait
  df_ <- pheno_df[!is.na(as.numeric(pheno_df[, trait_])), ]
  df_$Genotype <- as.factor(df_$Genotype)

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
  geno_mat <- geno_df[match(unique(df_$Genotype), geno_df$Genotype), ]

  # compute gblup ebv
  linear_krmm_model <- krmm(
    Y = as.numeric(df_[, trait_]),
    X = x_mat_fixed_eff, Z = z_mat_rand_eff,
    Matrix_covariates = geno_mat[, -1],
    method = "GBLUP"
  )
  
  # get gblup ebvs
  gblup_ebv <- predict_krmm(
    linear_krmm_model,
    Matrix_covariates = geno_mat[, -1],
    X = x_mat_fixed_eff,
    Z = z_mat_rand_eff,
    add_fixed_effects = F
  )

  gblup_ebv_df <- data.frame(
    "Genotype" = geno_mat$Genotype,
    "Trait" = gblup_ebv
  )
  colnames(gblup_ebv_df)[2] <- trait_

  list_gblup_ebv[[trait_]] <- gblup_ebv_df
}

# reduce gblup ebv list
gblup_ebv_df <- Reduce(
  function(x, y) {
    merge(x, y, by = "Genotype", all = T)
  },
  list_gblup_ebv
)

# save computed gblup ebv
colnames(gblup_ebv_df) <- c("Genotype", traits_)

fwrite(
  gblup_ebv_df,
  file = paste0(pheno_dir_path_, "gblup_ebv.csv")
)
