#bnmf aortic traits example run 



#----
start=Sys.time()


# load requires packages
install.packages("pacman")
pacman::p_load(tidyverse, data.table, readxl, magrittr, dplyr, strex,
               rstudioapi, DT, kableExtra, GenomicRanges)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges", force=TRUE)
BiocManager::install("Homo.sapiens", force=TRUE)

install.packages("future.apply")
library(future.apply)

library(dplyr)


# load project scripts containing bNMF functions
source("/Users/adama/Downloads/Tchandjieu_Lab/NMF_tutorial/scripts/choose_variants.R")  # ld_pruning, count_traits_per_variant, fina_variants_needing_proxies, & choose_potential_proxies
source("/Users/adama/Downloads/Tchandjieu_Lab/NMF_tutorial/scripts/prep_bNMF_2025.R")  # fetch_summary_stats & prep_z_matrix
source("/Users/adama/Downloads/Tchandjieu_Lab/NMF_tutorial/scripts/run_bNMF_2025.R")  # run_bNMF & summarize_bNMF

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library(readr)

df <- read_csv("/Users/adama/Downloads/Tchandjieu_Lab/NMF_python/post_prune/adama_bnmf_file.csv")

#load zscore file

data_dir_new = "/Users/adama/Downloads/Tchandjieu_Lab/real_bnmf"
bmnf_file_new <- file.path(data_dir_new, "real_bnmf_file.csv")

bnmf_file<- read.csv('/Users/adama/Downloads/Tchandjieu_Lab/NMF_python/all_2606__bnmf_matrix.csv', row.names = 1)



#zscores 
bnmf_reps_new <- run_bNMF(bnmf_file,
                      n_reps=25,
                      tolerance = 1e-6)

class(bnmf_reps_new)

#bnmf_subset <- bnmf_reps_new[1:2]  # use only first 10 reps
summarize_bNMF(bnmf_reps_new, dir_save = data_dir_new)

#summarize_bNMF(bnmf_reps_new, dir_save=data_dir_new)

save.image(file = file.path(project_dir, "pipeline_data.RData"))

end=Sys.time()
print("Total pipeline runtime:")
print(end-start)