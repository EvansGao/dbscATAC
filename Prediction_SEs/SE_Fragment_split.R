library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicFeatures)
#library(EnsDb.Hsapiens.v75)
library(devtools)
library(leidenbase)
library(remotes)
#library(cicero)
library(tidyr)
library(icesTAF)
#library(monocle3)
library(dplyr)
library(rtracklayer)
library(Matrix)
library(patchwork)
library(rlang)
library(cowplot)
library(lsa)
library(ggplot2)
library(argparse)
library(locStra)
library(Matrix)
#library(cicero)
library(tidyr)
library(stringr)

set.seed(1234)

#superenhancer
mkdir("fragment")
project <- "GSE111586"
species <- c("mm10")
allsamples <- c("30078704_celltype")
##celltypes <- c("Testes","Thymus","WholeBrain")
##celltypes <- c("Cerebellum","SmallIntestine","Kidney")
##celltypes <- c("Heart","BoneMarrow","Liver","Lung")
#celltypes <- c("Spleen","PreFrontalCortex")
celltypes <- c("LargeIntestine")
for(m in 1:length(celltypes)){
#m <- 1
rds_file_ls <- list.files(path=paste0("/data/gts/dbscATAC/project/",sep="",project,sep="","Yes/rds"), pattern=paste0("*._",sep="",celltypes[m],sep="","_*."))
fragment <- read.table(paste0("./", sep="", celltypes[m], sep="", "_fragments.sort.bed.gz"))
names(fragment) <- c("chrome","start","end","barcode","overlap")
gc()

#rds_file_ls <- c("30078704_tissue_WholeBrain.GSE149683.mm10.rds")

for(t in 1:length(rds_file_ls)){
#t <- 1
celltype_rds <- readRDS(paste0("rds/", sep="", rds_file_ls[t]))
names <- str_split(rds_file_ls[t], "\\.")
celltypename <- names[[1]][1]
onecelltype <- data.frame(barcode=colnames(celltype_rds), celltype = celltypename)

#for 30078704_tissue_WholeBrain.GSE149683.mm10.rds
#if(dim(celltype_rds)[2] >= 5000){
#onecelltype <- data.frame(barcode=colnames(celltype_rds)[1:5000], celltype = celltypename)
#}else{
#onecelltype <- data.frame(barcode=colnames(celltype_rds), celltype = celltypename)
#}

onecelltypefragment <- merge(fragment, onecelltype, all =FALSE)
#onecelltypefragmentano <- merge(combined_df, onecelltype, all =FALSE)
onecelltypefragment$cluster <- onecelltypefragment$overlap
onecelltypefragment$overlap <- paste0(onecelltypefragment$chrome,sep=":",onecelltypefragment$start,sep=":",onecelltypefragment$end)
onecelltypefragment <- onecelltypefragment[,c(-1)]
write.table(onecelltypefragment, paste0("fragment/", sep="", names[[1]][1], sep=".", names[[1]][2], sep=".", names[[1]][3], sep=".", "bed"), sep = '\t', row.names = F, col.names = F,quote = F)
rm(onecelltypefragment)
gc()
system("sync; echo 3 > /proc/sys/vm/drop_caches")
}

}
