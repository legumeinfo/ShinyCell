# --------------------------------------------------------------
# Extension of RenameGenesSeurat() from
# https://vertesy.github.io/Seurat.utils/
# https://github.com/satijalab/seurat/issues/1049
# --------------------------------------------------------------
# First go to your local ShinyCell directory
setwd("~/Documents/NCGR/LIS/ShinyCell-glyma/")

library(Seurat)
#source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/R/Seurat.Utils.R")
library(stringi)

RenameGenesSeurat2 <- function(obj, newnames, var_features, assay_name) {
  assay <- obj@assays[[assay_name]]
  if (nrow(assay) == length(newnames)) {
    if (length(assay@counts)) rownames(assay@counts) <- newnames
    if (length(assay@data)) rownames(assay@data) <- newnames
    if (length(assay@scale.data)) rownames(assay@scale.data) <- newnames
    if (length(assay@var.features)) assay@var.features <- var_features
    if (nrow(assay@meta.features)) rownames(assay@meta.features) <- newnames
  } else {"Unequal gene sets: nrow(assay) != nrow(newnames)"}
  obj@assays[[assay_name]] <- assay
  return(obj)
}

## Rename genes like "GLYMA-18G043200" to "glyma.Wm82.gnm2.ann1.Glyma.18G043200"
# Rename genes like "GLYMA-18G043200" to "Glyma.18G043200"
# TODO: rename shorter gene names from table
renameGenes <- function(seurat_filename) {
  seurat <- readRDS(seurat_filename)

  for (assay_name in Assays(seurat)) {
    gene_names <- rownames(seurat@assays[[assay_name]])
    gene_ids <- stri_match_first(gene_names, regex = "GLYMA[-_](.+)")[, 2]
    matches <- !is.na(gene_ids)
    #gene_names[matches] <- paste0("glyma.Wm82.gnm2.ann1.Glyma.", gene_ids[matches])
    gene_names[matches] <- paste0("Glyma.", gene_ids[matches])
    gene_names[!matches] <- gene_names[!matches]

    var_features <- seurat@assays[[assay_name]]@var.features
    gene_ids <- stri_match_first(var_features, regex = "GLYMA[-_](.+)")[, 2]
    matches <- !is.na(gene_ids)
    #var_features[matches] <- paste0("glyma.Wm82.gnm2.ann1.Glyma.", gene_ids[matches])
    var_features[matches] <- paste0("Glyma.", gene_ids[matches])
    var_features[!matches] <- var_features[!matches]

    seurat <- RenameGenesSeurat2(seurat, gene_names, var_features, assay_name)
  }
  new_filename <- paste0(stri_match_first(seurat_filename, regex = "(.+)\\.rds")[, 2], "-renamed.rds")
  saveRDS(seurat, new_filename)
}

renameGenes("notes/Nodules_new_clusters.rds")
renameGenes("notes/root_annotatin_new_clusters.rds")
# --------------------------------------------------------------
# Renaming clusters and deleting idents
library(dplyr)
seurat_nodules <- readRDS("notes/Nodules_new_clusters-renamed.rds")

# Drop redundant idents
seurat_nodules$ppod1ocol <- seurat_nodules$integrated_snn_res.0.3 <- seurat_nodules$seurat_clusters <- NULL
# Change 'newclusters' to 'clusters' and replace with biological cluster names
seurat_nodules$clusters <- seurat_nodules$newclusters
seurat_nodules$newclusters <- NULL
seurat_nodules$clusters <- recode(seurat_nodules$clusters,
  A = "A: Inner/Outer Cortex",
  B = "B: Inter/Outer Cortex",
  C = "C: Vascular Endodermis",
  D = "D: Sclereid Layer",
  E = "E: Vascular Bundle",
  F = "F: Infected",
  G = "G: Infected",
  H = "H: Infected, Nitrogen-Fixing",
  I = "I: Senescing Nodule",
  J = "J: Uninfected",
  K = "K: Uninfected"
)
saveRDS(seurat_nodules, "notes/Nodules_new_clusters-scrubbed.rds")

seurat_root <- readRDS("notes/root_annotatin_new_clusters-renamed.rds")
seurat_root$ppod1ocol <- seurat_root$integrated_snn_res.0.3 <- seurat_root$seurat_clusters <- NULL
seurat_root$clusters <- seurat_root$newclusters
seurat_root$newclusters <- NULL
seurat_root$clusters <- recode(seurat_root$clusters,
  `1` = "1: Epidermis",
  `2` = "2: Epidermis",
  `3` = "3: Epidermis",
  `4` = "4: Epidermis",
  `5` = "5: Dividing cells",
  `6` = "6: Dividing cells",
  `7` = "7: Cortex",
  `8` = "8: Cortex",
  `9` = "9: Cortex",
  `10` = "10: Endodermis",
  `11` = "11: Endodermis",
  `12` = "12: Cambium",
  `13` = "13: Pericycle",
  `14` = "14: Xylem",
  `15` = "15: Phloem fiber",
  `16` = "16: Phloem"
)
saveRDS(seurat_root, "notes/root_annotatin_new_clusters-scrubbed.rds")
# --------------------------------------------------------------
