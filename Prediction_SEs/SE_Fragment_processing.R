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
project <- "GSE149683"
species <- c("hg38")
#allsamples <- c("GSM4508928_adrenal_filtered")
#allsamples <- c("GSM4508937_muscle_filtered")
#allsamples <- c("GSM4508932_heart_filtered")
#celltypes <- c("adrenal","muscle")
celltypes <- c("eye","intestine")
#celltypes <- c("placenta","thymus","heart")
#celltypes <- c("kidney")
#celltypes <- c("cerebrum","cerebellum","stomach")
#celltypes <- c("spleen")
for(m in 1:length(celltypes)){
#m <- 1
rds_file_ls <- list.files(path="/data/gts/dbscATAC/project/GSE149683Yes/rds", pattern=paste0("*._",sep="",celltypes[m],sep="","_filtered*."))
fragment <- read.table(file = paste0("./",sep="",celltypes[m],sep="","_fragment.ori.bed.gz"))
names(fragment) <- c("chrome","start","end","barcode","overlap")

for(t in 1:length(rds_file_ls)){
fragfile <- gsub("\\.rds","",rds_file_ls[t])
output_file_c <- paste0("fragment",sep="/",fragfile,sep=".","bed")
if(file.exists(output_file_c)=="FALSE"){
celltype_rds <- readRDS(paste0("rds/", sep="", rds_file_ls[t]))
names <- str_split(rds_file_ls[t], "\\.")
celltypename <- names[[1]][1]
if(dim(celltype_rds)[2] >= 6000){
onecelltype <- data.frame(barcode=colnames(celltype_rds)[1:6000], celltype = celltypename)
}else{
onecelltype <- data.frame(barcode=colnames(celltype_rds), celltype = celltypename)
}

onecelltypefragment <- merge(fragment, onecelltype, all =FALSE)
onecelltypefragment$cluster <- onecelltypefragment$overlap
onecelltypefragment$overlap <- paste0(onecelltypefragment$chrome,sep=":",onecelltypefragment$start,sep=":",onecelltypefragment$end)
onecelltypefragment <- onecelltypefragment[,c(-1)]
write.table(onecelltypefragment, paste0("fragment/", sep="", names[[1]][1], sep=".", names[[1]][2], sep=".", names[[1]][3], sep=".", "bed"), sep = '\t', row.names = F, col.names = F,quote = F)
rm(onecelltypefragment)
gc()
system("sync; echo 3 > /proc/sys/vm/drop_caches")
}
}

}


##to here


#fragment_file_df <- lapply(fragment_file_ls, function(x) {read.table(file = paste0("sample_fragments/", sep="", x))})
#fragment_file_df_name <- lapply(fragment_file_df, function(x) {names(x) <- c("chrome","start","end","barcode","overlap")})
#fragment_file_df_name <- lapply(fragment_file_df, function(x){
#names(x) <- c("chrome","start","end","barcode","overlap")
#return(x)
#})
# Combine them
#fragment <- read.table(file = "./intestine_fragment.ori.bed.gz")
#rm(fragment_file_df)

#gc()
#fragment <- do.call("rbind", lapply(fragment_file_df, as.data.frame))
#names(fragment) <- c("chrome","start","end","barcode","overlap")

#gc()


for(t in 1:length(rds_file_ls)){
#t <- 2
fragfile <- gsub("\\.rds","",rds_file_ls[t])
output_file_c <- paste0("fragment",sep="/",fragfile,sep=".","bed")
if(file.exists(output_file_c)=="FALSE"){
celltype_rds <- readRDS(paste0("rds/", sep="", rds_file_ls[t]))
names <- str_split(rds_file_ls[t], "\\.")
celltypename <- names[[1]][1]
if(dim(celltype_rds)[2] >= 6000){
onecelltype <- data.frame(barcode=colnames(celltype_rds)[1:6000], celltype = celltypename)
}else{
onecelltype <- data.frame(barcode=colnames(celltype_rds), celltype = celltypename)
}
#onecelltype <- data.frame(barcode=colnames(celltype_rds), celltype = celltypename)
#onecelltypefragment <- do.call("rbind", lapply(fragment_file_df, function(x){
#return(merge(x, onecelltype, all =FALSE))
#}))
onecelltypefragment <- merge(fragment, onecelltype, all =FALSE)
onecelltypefragment$cluster <- onecelltypefragment$overlap
onecelltypefragment$overlap <- paste0(onecelltypefragment$chrome,sep=":",onecelltypefragment$start,sep=":",onecelltypefragment$end)
onecelltypefragment <- onecelltypefragment[,c(-1)]
write.table(onecelltypefragment, paste0("fragment/", sep="", names[[1]][1], sep=".", names[[1]][2], sep=".", names[[1]][3], sep=".", "bed"), sep = '\t', row.names = F, col.names = F,quote = F)
rm(onecelltypefragment)
gc()
system("sync; echo 3 > /proc/sys/vm/drop_caches")
}
}

