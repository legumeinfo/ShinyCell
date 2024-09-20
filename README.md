# &#127793; ShinyCell (LIS version)

Welcome to the [Legume Information System](https://www.legumeinfo.org) fork of `ShinyCell`. The legume seedling icon &#127793; indicates details specific to this version.

# ShinyCell package
`ShinyCell` is a R package that allows users to create interactive Shiny-based 
web applications to visualise single-cell data via (i) visualising cell 
information and/or gene expression on reduced dimensions e.g. UMAP, (ii) 
visualising the coexpression of two genes on reduced dimensions, (iii) 
visualising the distribution of continuous cell information e.g. nUMI / module 
scores using violin plots / box plots, (iv) visualising the composition of 
different clusters / groups of cells using proportion plots and (v) 
visualising the expression of multiple genes using bubbleplots / heatmap. 
Examples of ShinyCell-generated shiny apps for single and multi datasets can 
be found at http://shinycell1.ddnetbio.com and http://shinycell2.ddnetbio.com 
respectively.

If you are using `ShinyCell`, please cite [Ouyang et al. ShinyCell: Simple and 
sharable visualisation of single-cell gene expression data. Bioinformatics, 
doi:10.1093/bioinformatics/btab209](
http://dx.doi.org/10.1093/bioinformatics/btab209). The manuscript 
is recently accepted and we will update the full citation when it is available.

Key features of `ShinyCell` include: 

1. Written in R and uses the Shiny package, allowing for easy sharing on online 
   platforms e.g. [shinyapps.io](https://www.shinyapps.io/) and Amazon Web 
   Services (AWS) or be hosted via Shiny Server

2. Supports all of the major single-cell data formats (h5ad / loom / Seurat / 
   SingleCellExperiment) and we also include a simple tutorial to process 
   plain-text gene expression matrices
  
3. Web interface have low memory footprint due to the use of hdf5 file system 
   to store the gene expression. Only the expression of genes that are plotted 
   are loaded into memory

4. Inclusion of less common single-cell visualisations, namely coexpression 
   plots and proportion plots that provide additional information on top of 
   low-dimensional embeddings
  
5. Users can export visualisations as PDF or PNG images for presentation or 
   publication use

6. Ability to include multiple single-cell datasets into a single Shiny web app

7. It is easy to use and customise aethetsics e.g. label names and colour 
   palettes. In the simplest form, ShinyCell can convert an input single-cell 
   data into a Shiny app with five lines of code 
   (see [Quick Start Guide](#quick-start-guide))

We also compared ShinyCell with nine other popular scRNA-seq visualisation 
tools, which further highlights the key features of `ShinyCell`. For a more 
detailed description, see the 
[Supplementary Information](docs/OuyangEtAl_Shinycell_SuppInfo.pdf).

![](images/comparison.png)



# Table of Contents and Additional Tutorials
This readme is broken down into the following sections:

- [Installation](#installation) on how to install `ShinyCell`

- [Quick Start Guide](#quick-start-guide) to rapidly deploy a shiny app with 
  a few lines of code

- [Frequently Asked Questions](#frequently-asked-questions)

There are also additional tutorials as follows:

- [Tutorial for customising ShinyCell aesthetics](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)

- [Tutorial for creating a ShinyCell app containing several single-cell datasets](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/2multi.html)

- [Tutorial on other supported file formats (h5ad / loom / SCE)](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/3otherformat.html)

- [Tutorial on processing plain-text gene expression matrices](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/3plaintext.html)

- [Instructions on how to deploy ShinyCell apps online](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/4cloud.html)



# Installation
First, users can run the following code to check if the packages required by 
`ShinyCell` exist and install them if required:
``` r
reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

# If you are using h5ad file as input, run the code below as well
# reticulate::py_install("anndata")
```

Furthermore, on the system where the Shiny app will be deployed, users can run 
the following code to check if the packages required by the Shiny app exist 
and install them if required:
``` r
reqPkg = c("shiny", "shinyhelper", "shinyjs", "data.table", "Matrix", "DT", "hdf5r",
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro", "RCurl", "jsonlite")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}
```

`ShinyCell` can then be installed from GitHub as follows:
``` r
devtools::install_github("SGDDNB/ShinyCell")
```

&#127793; To install this LIS version instead, use
``` r
remove.packages("ShinyCell") # if necessary
devtools::install_github("legumeinfo/ShinyCell")
```



# Quick Start Guide
In short, the `ShinyCell` package takes in an input single-cell object and 
generates a ShinyCell config `scConf` containing labelling and colour palette 
information for the single-cell metadata. The ShinyCell config and single-cell 
object are then used to generate the files and code required for the shiny app. 

In this example, we will use single-cell data (Seurat object) containing 
intermediates collected during the reprogramming of human fibroblast into 
induced pluripotent stem cells using the RSeT media condition, taken from 
[Liu, Ouyang, Rossello et al. Nature (2020)](
https://www.nature.com/articles/s41586-020-2734-6). The Seurat object can be 
[downloaded here](http://files.ddnetbio.com/hrpiFiles/readySeu_rset.rds).

A shiny app can then be quickly generated using the following code:
 
``` r
library(Seurat)
library(ShinyCell)

getExampleData()                       # Download example dataset (~200 MB)
seu = readRDS("readySeu_rset.rds")
scConf = createConfig(seu)
makeShinyApp(seu, scConf, gene.mapping = TRUE,
             shiny.title = "ShinyCell Quick Start") 
```

The generated shiny app can then be found in the `shinyApp/` folder (which is 
the default output folder). To run the app locally, use RStudio to open either 
`server.R` or `ui.R` in the shiny app folder and click on "Run App" in the top 
right corner. The shiny app can also be deployed online via online platforms 
e.g. [shinyapps.io](https://www.shinyapps.io/) and Amazon Web Services (AWS) 
or be hosted via Shiny Server. For further details, refer to 
[Instructions on how to deploy ShinyCell apps online](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/4cloud.html).

&#127793; To generate our application [_Medicago truncatula_ Meliloti vs. Mock Inoculated Root](https://shinycell.legumeinfo.org/medtr.A17.gnm5.ann1_6.expr.Cervantes-Perez_Thibivilliers_2022/), we use
``` r
seurat_medtr <- readRDS("shinycell.rds")
scConf_medtr <- createConfig(seurat_medtr)
title_medtr <- "Medicago truncatula Meliloti vs. Mock Inoculated Root"
dir_medtr <- "medtr.A17.gnm5.ann1_6.expr.Cervantes-Perez_Thibivilliers_2022/"
props_medtr <- list(
  gene_prefix = "medtr.A17.gnm5.ann1_6.",
  gene_lookup_filename = "medicago_gene_id_lookup.txt"
)
saveRDS(props_medtr, paste0(dir_medtr, "properties.rds"))
makeShinyApp(seurat_medtr, scConf_medtr, gene.mapping = TRUE,
  shiny.title = title_medtr, shiny.dir = dir_medtr)
```
which puts the application files in `medtr.A17.gnm5.ann1_6.expr.Cervantes-Perez_Thibivilliers_2022/` instead of the default `shinyApp/`.

The shiny app contains seven tabs (highlighted in blue box), with the opening 
page showing the first tab "CellInfo vs GeneExpr" (see below), plotting both 
cell information and gene expression side-by-side on reduced dimensions e.g. 
UMAP. Users can click on the toggle on the bottom left corner to display the 
cell numbers in each cluster / group and the number of cells expressing a gene.
The next two tabs are similar, showing either two cell information 
side-by-side (second tab: "CellInfo vs CellInfo") or two gene expressions 
side-by-side (third tab: "GeneExpr vs GeneExpr").

![](images/quick-shiny1.png)

The fourth tab "Gene coexpression" blends the gene expression of two genes, 
given by two different colour hues, onto the same reduced dimensions plot. 
Furthermore, the number of cells expressing either or both genes are given. 

![](images/quick-shiny2.png)

The fifth tab "Violinplot / Boxplot" plots the distribution of continuous cell 
information e.g. nUMI or module scores or gene expression across each cluster 
/ group using violin plots or box plots.

![](images/quick-shiny3.png)

The sixth tab "Proportion plot" plots the composition of different clusters / 
groups of cells using proportion plots. Users can also plot the cell numbers 
instead of proportions.

![](images/quick-shiny4.png)

The seventh tab "Bubbleplot / Heatmap" allows users to visualise the 
expression of multiple genes across each cluster / group using bubbleplots / 
heatmap. The genes (rows) and groups (columns) can be furthered clustered 
using hierarchical clustering.

![](images/quick-shiny5.png)

## &#127793; LIS-specific properties

As part of your build process, create an R data file called `properties.rds` in your application folder, containing a list of (optional) fields described in the next two sections. If any are not specified, `ShinyCell` uses their default settings.

For an example, see the **_Medicago truncatula_ Meliloti vs. Mock Inoculated Root** generation process above. It defines `gene_prefix` and `gene_lookup_filename`, but not `tab_index` or `dataset_index`.

### &#127793; Specifying values in the URL

This LIS version allows setting initial values of
* `gene1`: gene to display in 1st gene expression view
* `gene2`: gene to display in 2nd gene expression view
* `tab`: tag for selected tab
* `dataset`: tag for selected dataset (ignored unless multiple datasets exist)

For example, this URL will launch the application in the **Gene coexpression** tab, with genes MtRPG and MtFE selected.<br>
https://shinycell.legumeinfo.org/medtr.A17.gnm5.ann1_6.expr.Cervantes-Perez_Thibivilliers_2022/?gene1=MtRPG&gene2=MtFE&tab=gene-coexpression

To configure the tab and dataset lookup, you may add the following (optional) lists to your `properties.rds`. (Make sure that they match your actual tabs and datasets!)

`tab_index`: List for converting URL `tab` field to the matching tab's index. If missing, `ShinyCell` assumes that your application has the original tab structure,
```
tab_index = list(
  "cell-gene" = 1,
  "cell-cell" = 2,
  "gene-gene" = 3,
  "gene-coexpression" = 4,
  "violin-box" = 5,
  "proportion-numbers" = 6,
  "bubble-heatmap" = 7
)
```

`dataset_index`: List for converting URL `dataset` field to the matching dataset's index, required only when your application has multiple datasets. An example would be
```
dataset_index = list(
  "nodules" = 1,
  "root" = 2
)
```

### &#127793; Gene linkouts

The **Gene expression** panels now display a set of linkouts for the selected gene. We obtain these from the LIS `gene_linkouts` microservice, which takes as argument the full-yuck gene name. Since the gene names in your Seurat object are generally short forms or represented by another symbol, you must add the following (optional) fields to your `properties.rds` to tell `ShinyCell` how to convert them to full-yuck.

`gene_lookup_filename` (optional): Name of a text file containing a lookup table of gene symbol to gene name. This file must contain two tab-separated columns; column 1 is the gene name and column 2 is the gene symbol(s). If multiple symbols point to the same gene, you may list them on the same row in column 2, separated by commas. If missing, `ShinyCell` will assume that all of your gene names really are names and not symbols, and skip the lookup.

`gene_prefix` (optional): Prefix that restores the gene names to full-yuck. If missing, `ShinyCell` will assume that the gene names are already full-yuck, and skip the prefixing step.


# Frequently Asked Questions
- Q: How much memory / storage space does `ShinyCell` and the Shiny app consume?
  - A: The Shiny app itself consumes very little memory and is meant to be a 
       heavy-duty app where multiple users can access the app simultaneously. 
       Unlike typical R objects, the entire gene expression matrix is stored 
       on disk and *not on memory* via the hdf5 file system. Also, the hdf5 
       file system offers superior file compression and takes up less storage 
       space than native R file formats such as rds / Rdata files.
  - A: It should be noted that a large amount of memory is required when 
       *building* the Shiny app. This is because the whole Seurat / SCE object 
       has to be loaded onto memory and additional memory is required to 
       generate the required files. From experience, a typical laptop with 8GB 
       RAM can handle datasets around 30k cells while 16GB RAM machines can 
       handle around 60k-70k cells. 
       
- Q: I have generated additional dimension reductions (e.g. force-directed 
layout / MDS / monocle2 etc.) and would like to include them into the Shiny 
app. How do I do that?
  - A: `ShinyCell` automatically retrieves dimension reduction information from 
       the Seurat or SCE object. Thus, the additional dimension reductions 
       have to be added into the Seurat or SCE object before running `ShinyCell`. 
  - For Seurat objects, users can refer "Storing a custom dimensional reduction 
       calculation" in https://satijalab.org/seurat/v3.0/dim_reduction_vignette.html
  - For SCE objects, users can refer to https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#3_Adding_low-dimensional_representations
  
- Q: I have both RNA and integrated data in my Seurat object. How do I specify 
which gene expression assay to plot in the Shiny app?
  - A: Only one gene expression assay can be visualised per dataset. To 
       specify the assay, use the `gex.assay = "integrated` argument in the 
       `makeShinyApp()` or `makeShinyFiles()` functions. If users want to 
       visualise both gene expression, they have to treat each assay as an 
       individual dataset and include multiple datasets into a single shiny 
       app, following the [Tutorial for creating a ShinyCell app containing several single-cell datasets](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/2multi.html)



<br/><br/>
<br/><br/>
<br/><br/>
<br/><br/>
<br/><br/>