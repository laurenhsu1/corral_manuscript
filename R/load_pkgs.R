library(SingleCellExperiment)
library(corral)
library(ggplot2)
library(Seurat)
library(scran)
library(patchwork)
library(reshape2)
library(scuttle)
library(scater)
library(bluster)
library(stringr)
library(DropletUtils)
library(devtools)

# data packages
try(library(SeuratData))
library(DuoClustering2018)
library(CellBench)


library(biomaRt)
try(library(EnsDb.Hsapiens.v79))