library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicFeatures)
library(EnsDb.Hsapiens.v75)
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

##
#atac <- readRDS("./GSM4508940_spleen_filtered.seurat.RDS")
#str(atac)

##get mtx
##Rscript SB_downloadTOcelltype_GSE149683.R>& ss2.txt &
mkdir("rds")
project <- "GSE149683"
species <- c("hg38","hg38","hg38","hg38","hg38","hg38","hg38","hg38","hg38","hg38","hg38","hg38","hg38","hg38","hg38")
allsamples <- c("GSM4508928_adrenal_filtered","GSM4508929_cerebellum_filtered","GSM4508930_cerebrum_filtered","GSM4508931_eye_filtered","GSM4508932_heart_filtered","GSM4508933_intestine_filtered","GSM4508934_kidney_filtered","GSM4508935_liver_filtered","GSM4508936_lung_filtered","GSM4508937_muscle_filtered","GSM4508938_pancreas_filtered","GSM4508939_placenta_filtered","GSM4508940_spleen_filtered","GSM4508941_stomach_filtered","GSM4508942_thymus_filtered")
gtfs <- c("hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz","hg38.gtf.gz")
source("/data/gts/dbscATAC/project/perl/ChangedFeaturePlot.R")

for(t in 9:length(allsamples)){
#t <- 8
atac <- readRDS(paste("/data/gts/dbscATAC/project/", project, "/", allsamples[t], ".seurat.RDS", sep=""))
mtxPeak <- atac@assays$peaks@counts
mtxPeak <- mtxPeak[Matrix::rowSums(mtxPeak) > 0,]
rownames(mtxPeak) <- gsub("\\-","_",rownames(mtxPeak))
mtxPeak <- mtxPeak[grep("^chr\\d+_\\d+_\\d+$|^chrX_\\d+_\\d+$|^chrY_\\d+_\\d+$", rownames(mtxPeak), value=TRUE),]
write.table(rownames(mtxPeak),"peaks", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
system("perl scATACpeak_hg19TOhg38.pl peaks")
peakinfo <- read.table("./peaks_hg38", header=FALSE)
mtxPeak <- mtxPeak[peakinfo$V1,]
#all(row.names(mtxPeak)==peakinfomm$V1)
row.names(mtxPeak) <- peakinfo$V2
mtxPeak <- mtxPeak[grep("^chr\\d+_\\d+_\\d+$|^chrX_\\d+_\\d+$|^chrY_\\d+_\\d+$", rownames(mtxPeak), value=TRUE),]

features <- atac@meta.data
#unique(atac@meta.data$cell_type)
#head(features)
#head(rownames(features))
features$cell_type <- gsub(" ","_",features$cell_type)
features$cell_type <- gsub("\\/","_",features$cell_type)
features$cell_type <- gsub("\\.","_",features$cell_type)
features$cell_type <- gsub("\\?","",features$cell_type)
features$cellType <- features$cell_type
features$cellId <- rownames(features)
tissues <- unique(features$cellType)
mtxPeaktmp <- mtxPeak

for(i in 1:length(tissues)){
#i<- 2
subcells <-features$cellId[features$cellType==tissues[i]]
tmpcluster <- mtxPeak[,subcells]
tmpcluster <- tmpcluster[Matrix::rowSums(tmpcluster) > 0,]
if(length(subcells)<100){
next
}
if(file.exists(paste0("rds/", sep="", allsamples[t], sep="_", tissues[i], sep=".", project, sep=".", species[t], sep="", ".rds"))=="FALSE"){
saveRDS(tmpcluster, file = paste0("rds/", sep="", allsamples[t], sep="_", tissues[i], sep=".", project, sep=".", species[t], sep="", ".rds"))
}
}

if(file.exists("peaks.txt")=="TRUE"){
file.remove('peaks.txt')
}
if(file.exists("barcodes.txt")=="TRUE"){
file.remove('barcodes.txt')
}
if(file.exists("counts.mm")=="TRUE"){
file.remove('counts.mm')
}
if(file.exists("integrated.fragments.tsv.gz")=="TRUE"){
file.remove('integrated.fragments.tsv.gz')
}
if(file.exists("integrated.fragments.tsv.gz.tbi")=="TRUE"){
file.remove('integrated.fragments.tsv.gz.tbi')
}

write(row.names(mtxPeaktmp),"peaks.txt")
write(colnames(mtxPeaktmp),"barcodes.txt")
writeMM(obj = Matrix(mtxPeaktmp , sparse = T), file="counts.mm")
#please run SSB_rawtoGetFragment.pl to get integrated.fragment.gz
system("perl SSB_rawtoGetFragment.pl")

gc()
#should change species ID
#head(rownames(mtxPeak))
rna_chromassay <- CreateChromatinAssay(
counts = mtxPeak,
sep = c("_", "_"),
#genome = species[t],
fragments = './integrated.fragments.tsv.gz',
)
rna <- CreateSeuratObject(counts = rna_chromassay, assay = "ATAC")
#please change ensembl database
#https://ftp.ensembl.org/pub/release-110/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.110.gtf.gz

gtfGRanges <- import(paste0("/data/gts/dbscATAC/project/gtf/", sep="", gtfs[t]))
gtfGRangestmp<- GRanges(
  seqnames = paste0("chr", seqnames(gtfGRanges)),
  ranges = ranges(gtfGRanges), strand=strand(gtfGRanges)
)
mcols(gtfGRangestmp) <- mcols(gtfGRanges)
names(gtfGRangestmp) <- mcols(gtfGRanges)$gene_id
#danRer10new <- danRer10tmp[seqnames(danRer10tmp) %in% standardChromosomes(danRer10tmp)[1:25]]
##note the genome build
genome(gtfGRangestmp) <- species[t]
#seqlevels(gtfGRangestmp) <- paste0('chr', seqlevels(gtfGRangestmp))
#seqlevelsStyle(gtfGRangestmp) <- "UCSC"
#gtfGRangestmp$gene_biotype <- gtfGRangestmp$type
Annotation(rna) <- gtfGRangestmp

rnanew <- RunTFIDF(rna)
#slotNames(rnanew)
rnanew <- FindTopFeatures(rnanew, min.cutoff = 'q75')
rnanew <- RunSVD(object = rnanew)

rnanew <- RunUMAP(object = rnanew, reduction = 'lsi', dims = 2:50)
rnanew <- FindNeighbors(object = rnanew, reduction = 'lsi', dims = 2:50)
rnanew <- FindClusters(object = rnanew, verbose = FALSE, resolution = 1, algorithm = 3)

gene.activities <- GeneActivity(rnanew)
#gene.activities <- GeneActivity(rnanew, biotypes = "gene")
rnanew[['RNA']] <- CreateAssayObject(counts = gene.activities)
rnanew <- NormalizeData(
  object = rnanew,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(rnanew$nCount_RNA)
)
#fragments[[rnanew]]

DefaultAssay(rnanew) <- 'RNA'
#head(rnanew@meta.data)
#get cell type information
cell_type <- features$cellType
#cell_type <- gsub(" ","_",cell_type)
#cell_type <- gsub("\\.","_",cell_type)
#cell_type <- gsub("\\?","",cell_type)
#cell_type <- gsub("\\/","_",cell_type)
rnanew$cell_type <- cell_type
Idents(rnanew) <- rnanew@meta.data$cell_type
rnanew.markers <- FindAllMarkers(rnanew, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.585, return.thresh=1e-10)
saveRDS(rnanew.markers, file=paste0("./", sep="", "rnanewmarkers", sep="_", allsamples[t], sep=".", project, sep=".", species[t], sep="", ".rds"))

#rnanew.markers <- readRDS("./rnanewmarkers.R")
write.table(rnanew.markers, paste0("./", sep="", "rnanewmarkers", sep="_", allsamples[t], sep=".", project, sep=".", species[t], sep="", ".txt"), sep = '\t', row.names = T, col.names = T, quote = F)
png(filename=paste0("./", sep="", allsamples[t], sep="_", "FeaturePlot", sep=".", project, sep=".", species, sep="", ".png"), width=2000, height=600)
p1 <- DimPlot(rnanew, reduction = "umap", pt.size = 0.8, label=TRUE, repel = TRUE, group.by = "cell_type") + NoLegend()
p2 <- FeaturePlot(rnanew, reduction = "umap", pt.size = 0.8, cols =c("lightgrey", "red"), features = rnanew.markers$gene[1:2], min.cutoff ="q1", ncol=2)
p1 + p2
dev.off()

genes <- unique(subset(rnanew.markers, rnanew.markers$p_val_adj<=1e-10)$gene)
rnanewsub <- rnanew[genes,]
datamtx <- GetAssayData(rnanewsub, slot = "data")
#head(rownames(datamtx))
rnanewsub <- CreateSeuratObject(counts = datamtx, assay = "RNA", meta.data =rnanew@meta.data)
rnanewsub@reductions <- rnanew@reductions
#DefaultDimReduc(object = rnanewsub)
saveRDS(rnanewsub, file=paste0("./", sep="", "rnanewsub", sep="_", allsamples[t], sep=".", project, sep=".", species[t], sep="", ".rds"))
Idents(rnanewsub) <- rnanewsub@meta.data$cell_type
#source("/data/gts/dbscATAC/project/perl/ChangedFeaturePlot.R")
#rnanewsub <- readRDS("rnanewsub_GSM5901075_uterus.GSE196794.macFas6.rds")
png(filename=paste0("./", sep="", allsamples[t], sep="_", "ChangedFeaturePlot", sep=".", project, sep=".", species[t], sep="", ".png"), width=500, height=500)
ChangedFeaturePlot(rnanewsub, pt.size =1.5, features = rnanew.markers$gene[1])
dev.off()
if(file.exists("peaks.txt")=="TRUE"){
file.remove('peaks.txt')
file.remove('barcodes.txt')
file.remove('counts.mm')
file.remove('integrated.fragments.tsv.gz')
file.remove('integrated.fragments.tsv.gz.tbi')
}
gc()

}

##to here


