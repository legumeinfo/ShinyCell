# --------------------------------------------------------------
# First go to your local ShinyCell directory
setwd("~/Documents/NCGR/LIS/ShinyCell-glyma/")
# --------------------------------------------------------------
# If you modified the ShinyCell package itself, remove and reinstall it
# remove.packages("ShinyCell")
# local build:
# devtools::install() # pkg = "."
# or from repository:
# devtools::install_github("legumeinfo/ShinyCell")
# --------------------------------------------------------------
# To build your ShinyCell application (example)
library(Seurat)
library(ShinyCell)

dir_glyma <- "glyma.expr.Cervantes-Perez_Zogli_2024/"
title_glyma <- "Glycine max Nodules vs. Root Seedlings"

#seurat_nodules <- readRDS("notes/Nodules_new_clusters-renamed.rds") # after renaming genes
seurat_nodules <- readRDS("notes/Nodules_new_clusters-scrubbed.rds") # after renaming genes and clusters
scConf_nodules <- createConfig(seurat_nodules)
makeShinyFiles(seurat_nodules, scConf_nodules, shiny.prefix = "sc1", shiny.dir = dir_glyma)

#seurat_root <- readRDS("notes/root_annotatin_new_clusters-renamed.rds") # after renaming genes
seurat_root <- readRDS("notes/root_annotatin_new_clusters-scrubbed.rds") # after renaming genes and clusters
scConf_root <- createConfig(seurat_root)
makeShinyFiles(seurat_root, scConf_root, shiny.prefix = "sc2", shiny.dir = dir_glyma)

# Properties file, run only once unless changed
# props_glyma <- list(
#   gene_prefix = "glyma.Wm82.gnm2.ann1.",
#   gene_lookup_filename = "soybean_ensembl_gene_id_lookup.txt",
#   dataset_index = list(
#     nodules = 1,
#     root = 2
#   )
# )
# saveRDS(props_glyma, paste0(dir_glyma, "properties.rds"))

citation = list(
  author  = "Cervantes-Perez SA, Zogli P, Amini S, Thibivilliers S, et al.",
  title   = "Single-cell transcriptome atlases of soybean root and mature nodule reveal new regulatory programs that control the nodulation process",
  journal = "Plant Communication",
  volume  = "5(8)",
  page    = "100984",
  year    = "2024", 
  doi     = "10.1016/j.xplc.2024.100984",
  link    = "https://pubmed.ncbi.nlm.nih.gov/38845198/"
)

makeShinyCodesMulti(
  shiny.title = title_glyma,
  shiny.footnotes = citation,
  shiny.prefix = c("sc1", "sc2"),
  shiny.headers = c("Nodules", "Root Seedlings"),
  shiny.dir = dir_glyma
)

# --------------------------------------------------------------
# To run your ShinyCell application locally
library(shiny)

shinyAppDir(dir_glyma, options = list(launch.browser = TRUE))
# --------------------------------------------------------------
