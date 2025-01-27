library(Signac)
library(Seurat)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v86)
#BiocManager::install("EnsDb.Hsapiens.v86")
#library(EnsDb.Hsapiens.v75)
#BiocManager::install("EnsDb.Hsapiens.v75")
#conda install -c bioconda bioconductor-ensdb.hsapiens.v75
library(ggplot2)
library(stringr)
library(patchwork)
library(icesTAF)
library(rtracklayer)
library(Matrix)
library(patchwork)
library(rlang)
library(cowplot)
library(lsa)
library(cisTopic)
set.seed(1234)

##get mtx
##Rscript SB_downloadTOcelltype_mtx_GSE163697.R>& ss.txt &
mkdir("rds")
project <- "GSE163697"
species <- c("dm6","dm6","dm6","dm6","dm6","dm6","dm6")
allsamples <- c("larva","pupa_0h","pupa_3h","pupa_6h","pupa_12h","pupa_72h","adult")
gtfs <- c("dm6.gtf.gz","dm6.gtf.gz","dm6.gtf.gz","dm6.gtf.gz","dm6.gtf.gz","dm6.gtf.gz","dm6.gtf.gz")
ataccounts <- Matrix::readMM("./cntMat_scATAC_240919c_129078r.mtx.gz")
#dim(ataccounts)
#head(rownames(ataccounts))
#cellinfo <- read.table(paste("/data/gts/dbscATAC/project/", project, "/scATAC_cells__cntMatRows.txt.gz", sep=""), skip = 1, header=FALSE)
cellinfo <- read.table("./scATAC_cells__cntMatRows.txt.gz", skip = 1, header=FALSE)
#head(cellinfo)
#length(cellinfo$V1)
colnames(ataccounts) <- cellinfo$V1
#peakinfo <- read.table(paste("/data/gts/dbscATAC/project/", project, "/scATAC_regions__cntMatCols.txt.gz", sep=""), skip = 1, header=FALSE)
peakinfo <- read.table("./scATAC_regions__cntMatCols.txt.gz", skip = 1, header=FALSE)
#head(peakinfo)
#length(peakinfo$V1)
rownames(ataccounts) <- peakinfo$V1

tsv <- read.table("./atacmetadatatrue.txt", sep="\t", header=FALSE)

source("/data/gts/dbscATAC/project/perl/ChangedFeaturePlot.R")
for(t in 2:length(allsamples)){
#t <- 1
tsvtmp <- subset(tsv, V3 == allsamples[t])
mtxPeak <- ataccounts[, tsvtmp$V1]
mtxPeak <- mtxPeak[Matrix::rowSums(mtxPeak) > 0,]
rownames(mtxPeak) <- gsub("\\-|\\:","_",rownames(mtxPeak))
features <- tsvtmp
names(features) <- c("cellId","cellType","stage")
features$cellType <- gsub("\\?","_",features$cellType)
tissues <- unique(features$cellType)
mtxPeaktmp <- mtxPeak

for(i in 1:length(tissues)){
#i<- 2
subcells <-features$cellId[features$cellType==tissues[i]]
tmpcluster <- mtxPeak[,subcells]
tmpcluster <- tmpcluster[Matrix::rowSums(tmpcluster) > 0,]
if(length(subcells)<50){
next
}
if(file.exists(paste0("rds/", sep="", allsamples[t], sep="_", tissues[i], sep=".", project, sep=".", species[t], sep="", ".rds"))=="FALSE"){
saveRDS(tmpcluster, file = paste0("rds/", sep="", allsamples[t], sep="_", tissues[i], sep=".", project, sep=".", species[t], sep="", ".rds"))
}
}

##to here




