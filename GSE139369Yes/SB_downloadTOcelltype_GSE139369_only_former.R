library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicFeatures)
#library(EnsDb.Hsapiens.v75)
library(devtools)
library(leidenbase)
library(remotes)
library(cicero)
library(tidyr)
library(icesTAF)
library(monocle3)
library(dplyr)
library(rtracklayer)
library(Matrix)
library(patchwork)
library(rlang)
library(cowplot)
library(lsa)
library(ggplot2)
set.seed(1234)
#library(agricolae)

#https://stackoverflow.com/questions/58310969/mean-of-the-rows-where-the-row-name-is-same-matrix-in-r
#atac <- readRDS("/data/gts/dbscATAC/project/GSE139369Yes/rds/Healthy-Hematopoiesis_07_GMP.GSE139369.hg38.rds")
#duplicate_rownames <- duplicated(rownames(atac))
#mat_unique <- atac[!duplicate_rownames, ]
#head(rownames(mat_unique))


#aa <- aggregate(atac, list(row.names(atac)), mean)
#aa <- do.call(rbind, by(atac, row.names(atac), FUN = colMeans))
#write.table(rownames(aa), 'example3.txt',  sep = '_', row.names = F, col.names = F, quote = F)
#write.table(rownames(atac), 'example2.txt',  sep = '_', row.names = F, col.names = F, quote = F)

##get mtx
##Rscript SB_downloadTOcelltype_GSE139369.R>& yyy5 &
mkdir("rds")
project <- "GSE139369"
species <- c("hg38","hg38")
allsamples <- c("All-Hematopoiesis-MPAL","Healthy-Hematopoiesis")
gtfs <- c("hg38.gtf.gz","hg38.gtf.gz")
source("/data/gts/dbscATAC/project/perl/ChangedFeaturePlot.R")
groups <- c("Group","BioClassification")

for(t in 1:length(allsamples)){
#t <- 1
atac <- readRDS(paste("/data/gts/dbscATAC/project/", project, "Yes/scATAC-", allsamples[t], "-191120.rds", sep=""))
if(groups[t]=="BioClassification"){
cellinfo <- data.frame(cell=colnames(atac), celltype=colData(atac)$BioClassification)
}else if(groups[t]=="Group"){
cellinfo <- data.frame(cell=colnames(atac), celltype=colData(atac)$Group)
}
#write.table(cellinfo, 'cellinfo.txt',  sep = '\t', row.names = F, col.names = F, quote = F)
peakGRanges <- rowRanges(atac)
peakinfo <- data.frame(seqnames=seqnames(peakGRanges),
  starts=start(peakGRanges)-1,
  ends=end(peakGRanges))
peakinfo$site_name <- paste(peakinfo$seqnames, peakinfo$starts, peakinfo$ends, sep="_")
#write.table(peakinfo, 'peakinfo.txt',  sep = '_', row.names = F, col.names = F, quote = F)
mtxPeak <- assay(atac)
rownames(mtxPeak) <- peakinfo$site_name
colnames(mtxPeak) <- cellinfo$cell

mtxPeak <- mtxPeak[grep("^chr\\d+_\\d+_\\d+$|^chrX_\\d+_\\d+$|^chrY_\\d+_\\d+$", rownames(mtxPeak), value=TRUE),]
write.table(rownames(mtxPeak),"peaks", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
system("perl scATACpeak_hg19TOhg38.pl peaks")
peakinfo <- read.table("./peaks_hg38", header=FALSE)
mtxPeak <- mtxPeak[peakinfo$V1,]
#all(row.names(mtxPeak)==peakinfomm$V1)
rownames(mtxPeak) <- peakinfo$V2
mtxPeak <- mtxPeak[grep("^chr\\d+_\\d+_\\d+$|^chrX_\\d+_\\d+$|^chrY_\\d+_\\d+$", rownames(mtxPeak), value=TRUE),]

#features <- data.frame(cell_type =cellinfo$celltype)
features <- cellinfo
#unique(atac@meta.data$cell_type)
#head(features)
#head(rownames(features))
features$cellType <- gsub("\\.","_",features$celltype)
features$cellId <- features$cell
tissues <- unique(features$cellType)

mtxPeak <- mtxPeak[!duplicated(rownames(mtxPeak)), ]
mtxPeaktmp <- mtxPeak

for(i in 1:length(tissues)){
#i<- 2
subcells <-features$cellId[features$cellType==tissues[i]]
if(length(subcells)<100){
next
}
tmpcluster <- mtxPeak[,subcells]
tmpcluster <- tmpcluster[Matrix::rowSums(tmpcluster) > 0,]
#if(file.exists(paste0("rds/", sep="", allsamples[t], sep="_", tissues[i], sep=".", project, sep=".", species[t], sep="", ".rds"))=="FALSE"){
saveRDS(tmpcluster, file = paste0("rds/", sep="", allsamples[t], sep="_", tissues[i], sep=".", project, sep=".", species[t], sep="", ".rds"))
#}
}

gc()


}

##to here


