library(devtools)
library(leidenbase)
library(remotes)
library(cicero)
library(tidyr)
library(icesTAF)
library(monocle3)
library(dplyr)
library(GenomeInfoDb)
library(Seurat)
library(Signac)
library(patchwork)
#BiocManager::install("GenomeInfoDb")
library(reticulate)
library(GenomicFeatures)
mkdir("rds")
##scATAC-Healthy-Hematopoiesis-191120.rds
ATAC <- readRDS("/data/gts/SingleCell/Cicero/31792411/scATAC-Healthy-Hematopoiesis-191120.rds")
#str(ATAC)
#head(colData(ATAC))
ATACM <- assay(ATAC)
cellinfo <- data.frame(cell=colnames(ATAC), celltype=colData(ATAC)$BioClassification)
#write.table(cellinfo, 'cellinfo.txt',  sep = '\t', row.names = F, col.names = F, quote = F)
peakGRanges <- rowRanges(ATAC)
peakinfo <- data.frame(seqnames=seqnames(peakGRanges),
  starts=start(peakGRanges)-1,
  ends=end(peakGRanges))
peakinfo$site_name <- paste(peakinfo$seqnames, peakinfo$starts, peakinfo$ends, sep="_")
#write.table(peakinfo, 'peakinfo.txt',  sep = '\t', row.names = F, col.names = F, quote = F)

row.names(ATACM) <- peakinfo$site_name
colnames(ATACM) <- cellinfo$cell

#head(assay(ATAC))
#head(metadata(ATAC))
#head(assay(ATAC))
#head(colData(ATAC))
#head(colData(ATAC)$Clusters2)
#head(rowRanges(ATAC)[.1:2])

tissues <- unique(cellinfo$celltype)
for(i in 1:length(tissues)){
subcells <- cellinfo$cell[cellinfo$celltype==tissues[i]]
if(length(subcells)<50 || grepl("Unk",tissues[i],ignore.case = TRUE)=="TRUE"){
next
}
tmpcluster <- ATACM[,subcells]
tmpcluster <- tmpcluster[Matrix::rowSums(tmpcluster) > 0,]
tmpcluster <- tmpcluster[,Matrix::colSums(tmpcluster) >= 200]
#saveRDS(tmpcluster, file = paste("rds/", tissues[i], "_healthy_31792411.rds", sep=""))
if(dim(tmpcluster)[2]<=2000){
saveRDS(tmpcluster, file = paste("rds/", tissues[i], "_healthy_31792411.rds", sep=""))
}else{
samplemtx <- CreateSeuratObject(counts = tmpcluster,project = 'ATAC')
samplemtx <- RunTFIDF(samplemtx)
samplemtx <- FindTopFeatures(samplemtx, min.cutoff = 'q0')
samplemtx <- RunSVD(object = samplemtx, reduction.key = 'LSI_', n = 30, reduction.name = 'lsi')
samplemtx <- RunUMAP(object = samplemtx, reduction = 'lsi', umap.method="uwot", dims = 2:30)
samplemtx <- FindNeighbors(object = samplemtx, reduction = 'lsi', dims = 2:30)
samplemtxnew <- FindClusters(object = samplemtx, verbose = FALSE, resolution = 0.2, algorithm = 3)
#DimPlot(object = samplemtxnew, reduction = "umap", pt.size = 0.8, label = TRUE, label.size = 8) + NoLegend()
features <- samplemtxnew@meta.data
#head(features)
#head(rownames(features))
clusters <- unique(features$seurat_clusters)
	for(j in 1:length(clusters)){
	subclusters <- rownames(features)[features$seurat_clusters==clusters[j]]
	#length(subcells)
	if(length(subclusters)<50){
	next
	}
	tmpclusterdata <- tmpcluster[,subclusters]
	tmpclusterdata <- tmpclusterdata[Matrix::rowSums(tmpclusterdata) > 0,]
	tmpclusterdata <- tmpclusterdata[,Matrix::colSums(tmpclusterdata) >= 200]
	saveRDS(tmpclusterdata, file = paste("rds/", tissues[i], "_", clusters[j], "_healthy_31792411.rds", sep=""))
	#print(head(subcells))
	}
}

}




##scATAC-Healthy-Hematopoiesis-191120.rds
ATAC <- readRDS("/data/gts/SingleCell/Cicero/31792411/scATAC-All-Hematopoiesis-MPAL-191120.rds")
#str(ATAC)
#head(colData(ATAC))
#tissues <- unique(colData(ATAC)$Group)
ATACM <- assay(ATAC)
cellinfo <- data.frame(cell=colnames(ATAC), celltype=colData(ATAC)$Group)
#write.table(cellinfo, 'cellinfo.txt',  sep = '\t', row.names = F, col.names = F, quote = F)
peakGRanges <- rowRanges(ATAC)
peakinfo <- data.frame(seqnames=seqnames(peakGRanges),
  starts=start(peakGRanges)-1,
  ends=end(peakGRanges))
peakinfo$site_name <- paste(peakinfo$seqnames, peakinfo$starts, peakinfo$ends, sep="_")
#write.table(peakinfo, 'peakinfo.txt',  sep = '\t', row.names = F, col.names = F, quote = F)

row.names(ATACM) <- peakinfo$site_name
colnames(ATACM) <- cellinfo$cell

#head(assay(ATAC))
#head(metadata(ATAC))
#head(assay(ATAC))
#head(colData(ATAC))
#head(colData(ATAC)$Clusters2)
#head(rowRanges(ATAC)[.1:2])

tissues <- unique(cellinfo$celltype)
for(i in 1:length(tissues)){
subcells <- cellinfo$cell[cellinfo$celltype==tissues[i]]
if(length(subcells)<50 || grepl("MPAL",tissues[i],ignore.case = TRUE)=="FALSE"){
next
}
tmpcluster <- ATACM[,subcells]
tmpcluster <- tmpcluster[Matrix::rowSums(tmpcluster) > 0,]
tmpcluster <- tmpcluster[,Matrix::colSums(tmpcluster) >= 200]
#saveRDS(tmpcluster, file = paste("rds/", tissues[i], "_both_31792411.rds", sep=""))
if(dim(tmpcluster)[2]<=2000){
saveRDS(tmpcluster, file = paste("rds/", tissues[i], "_31792411.rds", sep=""))
}else{
samplemtx <- CreateSeuratObject(counts = tmpcluster,project = 'ATAC')
samplemtx <- RunTFIDF(samplemtx)
samplemtx <- FindTopFeatures(samplemtx, min.cutoff = 'q0')
samplemtx <- RunSVD(object = samplemtx, reduction.key = 'LSI_', n = 30, reduction.name = 'lsi')
samplemtx <- RunUMAP(object = samplemtx, reduction = 'lsi', umap.method="uwot", dims = 2:30)
samplemtx <- FindNeighbors(object = samplemtx, reduction = 'lsi', dims = 2:30)
samplemtxnew <- FindClusters(object = samplemtx, verbose = FALSE, resolution = 0.2, algorithm = 3)
#DimPlot(object = samplemtxnew, reduction = "umap", pt.size = 0.8, label = TRUE, label.size = 8) + NoLegend()
features <- samplemtxnew@meta.data
#head(features)
#head(rownames(features))
clusters <- unique(features$seurat_clusters)
	for(j in 1:length(clusters)){
	subclusters <- rownames(features)[features$seurat_clusters==clusters[j]]
	#length(subcells)
	if(length(subclusters)<50){
	next
	}
	tmpclusterdata <- tmpcluster[,subclusters]
	tmpclusterdata <- tmpclusterdata[Matrix::rowSums(tmpclusterdata) > 0,]
	tmpclusterdata <- tmpclusterdata[,Matrix::colSums(tmpclusterdata) >= 200]
	saveRDS(tmpclusterdata, file = paste("rds/", tissues[i], "_", clusters[j], "_31792411.rds", sep=""))
	#print(head(subcells))
	}
}

}

