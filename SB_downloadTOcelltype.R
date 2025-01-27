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
library(icesTAF)
library(cisTopic)
library(patchwork)
set.seed(1234)

##get RDS

atac <-  readRDS("AllTimepoints_cisTopic.Rds")
write.table(atac@cell.data, "./atacmetadata.txt", sep = '\t', row.names = T, col.names = T,quote = F)
#str(atac)
mtxPeak <- atac@assays$data@listData$counts
rownames(mtxPeak) <- atac@rowRanges$peak
mtxPeak <- mtxPeak[Matrix::rowSums(mtxPeak) > 0,]
rownames(mtxPeak) <- gsub("\\-","_",rownames(mtxPeak))
mtxPeak <- mtxPeak[grep("^chr\\d+_\\d+_\\d+$|^chrX_\\d+_\\d+$|^chrY_\\d+_\\d+$", rownames(mtxPeak), value=TRUE),]
features <- data.frame(barcode=colnames(mtxPeak), cell_type = atac@colData@listData$assigned)
features$cell_type <- gsub(" ","_",features$cell_type)
features$cell_type <- gsub("\\.","_",features$cell_type)
features$cell_type <- gsub("\\?","",features$cell_type)
features$cell_type <- gsub("\\/","_",features$cell_type)
tissues <- unique(features$cell_type)
celltypes <- array()
j <- 1
for(i in 1:length(tissues)){
subcells <- features$barcode[features$cell_type==tissues[i]]
#print(length(subcells))
if(length(subcells)<50){
#if(length(subcells)<50 || grepl("Unk",tissues[i],ignore.case = TRUE)=="TRUE"){
next
}
tmpcluster <- mtxPeak[,subcells]
tmpcluster <- tmpcluster[Matrix::rowSums(tmpcluster) > 0,]
tmpcluster <- tmpcluster[,Matrix::colSums(tmpcluster) >= 200]
saveRDS(tmpcluster, file = paste("rds/", tissues[i], "_GSE151230.rds", sep=""))
celltypes[j] <- tissues[i]
j <- j+1
}




#superenhancer
#head(features)
celltype <- data.frame(barcode=rownames(features), celltype = features$cell_type)
fragment <- read.table("/data/xianym/singlecell/signac/atac_v1_pbmc_10k_fragments.tsv.gz")
names(fragment) <- c("chrome","start","end","barcode","overlap")

for(i in 1:length(celltypes)){
onecelltype <- subset(celltype, celltype==celltypes[i])
onecelltypefragment <- merge(fragment, onecelltype, all =FALSE)
onecelltypefragment$cluster <- onecelltypefragment$overlap
onecelltypefragment$overlap <- paste0(onecelltypefragment$chrome,sep=":",onecelltypefragment$start,sep=":",onecelltypefragment$end)
onecelltypefragment <- onecelltypefragment[,c(-1)]
write.table(onecelltypefragment, paste0("fragment", sep="/", celltypes[i],sep="", "_10X_v1_pbmc_10k.bed"), sep = '\t', row.names = F, col.names = F,quote = F)
#print(celltypes[i])
}

write.table(gsub(" ","_",celltypes), paste0("/data/gts/dbscATAC/superenhancer/filelist",sep="/","10X_v1_pbmc_10k_celltypes.txt"), sep = '\t', row.names = F, col.names = F,quote = F)

#bedtools sort -i tmp/cluster_10.bed>tmp/cluster_10_sorted.bed
#bedToBam -i tmp/cluster_10_sorted.bed -g hg19.chrom.sizes > tmp/cluster_10_sorted.bam
#samtools view tmp/cluster_10_sorted.bam | head -5



