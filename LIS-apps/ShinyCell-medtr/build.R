# --------------------------------------------------------------
# First go to your local ShinyCell directory
setwd("~/Documents/NCGR/LIS/ShinyCell-medtr/")
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

dir_medtr <- "medtr.A17.gnm5.ann1_6.expr.Cervantes-Perez_Thibivilliers_2022/"
title_medtr <- "Medicago truncatula Meliloti vs. Mock Inoculated Root"

seurat_medtr <- readRDS("notes/shinycell.rds")
scConf_medtr <- createConfig(seurat_medtr)

# Properties file, run only once unless changed
# props_medtr <- list(
#   gene_prefix = "medtr.A17.gnm5.ann1_6.",
#   gene_lookup_filename = "medicago_gene_id_lookup.txt"
# )
# saveRDS(props_medtr, paste0(dir_medtr, "properties.rds"))

makeShinyApp(seurat_medtr, scConf_medtr, shiny.title = title_medtr, shiny.dir = dir_medtr)
# --------------------------------------------------------------
# To run your ShinyCell application locally
library(shiny)

shinyAppDir(dir_medtr, options = list(launch.browser = TRUE))
# --------------------------------------------------------------
