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


