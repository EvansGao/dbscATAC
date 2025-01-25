#! /usr/bin/env Rscript

#This is for human
#usage example: ./enhancers_test.R -m "Lung_Hematopoietic_progenitors"
#The input file should be .rds scATAC-seq matrix
# four places should be changes input plot_cells

start_time <- Sys.time()

library(argparse)
library(locStra)
library(dplyr)
library(Matrix)
library(cicero)
library(tidyr)
library(stringr)
#install.packages('interp', dependencies=TRUE)
parser <- ArgumentParser()

parser$add_argument("-m", "--matrix", type= "character" , help = "single-cell ATAC-seq matrix")

args <- parser$parse_args()

input <- args$matrix
#input <- "Astrocyte.GSE192772.hg38"
#mkdir("R")
#mkdir("hundred")
#mkdir("peak")
names <- str_split(input, "\\.")
species <- names[[1]][3]

input_file <- paste0("rds",sep="/",input,sep=".","rds")
#input_file <- paste0("rds",sep="/","embryo_2-4_Ectoderm_A_29539636",sep=".","rds")
output_file_j <- paste0("R",sep="/",input,sep="_","jaccard",sep=".","txt")

output_file_c <- paste0("R",sep="/",input,sep="_","cicero",sep=".","txt")

scATAC_matrix <- readRDS(input_file)
#head(rownames(scATAC_matrix))
#head(colnames(scATAC_matrix))
scATAC_matrix <- 1*scATAC_matrix
#any(is.na(scATAC_matrix))
#sum(is.na(scATAC_matrix))
#scATAC_matrix[is.na(scATAC_matrix)] <- 0

#if(dim(scATAC_matrix)[2]>3000){
#T1 <- round(nrow(scATAC_matrix)/5, digits = 0)
#scATAC_matrix[1:T1,][scATAC_matrix[1:T1,] > 0] <- 1
#scATAC_matrix[(T1+1):(2*T1),][scATAC_matrix[(T1+1):(2*T1),] > 0] <- 1
#scATAC_matrix[(2*T1+1):(3*T1),][scATAC_matrix[(2*T1+1):(3*T1),] > 0] <- 1
#scATAC_matrix[(3*T1+1):(4*T1),][scATAC_matrix[(3*T1+1):(4*T1),] > 0] <- 1
#scATAC_matrix[(4*T1+1):nrow(scATAC_matrix),][scATAC_matrix[(4*T1+1):nrow(scATAC_matrix),] > 0] <- 1
#}else{
#scATAC_matrix[scATAC_matrix > 0] <- 1
#}
if(dim(scATAC_matrix)[2]<3000){
scATAC_matrix[scATAC_matrix > 0] <- 1
}

clean <- scATAC_matrix[Matrix::rowSums(scATAC_matrix) > 0,]
rownames(clean) <- gsub("\\:|\\-","_",rownames(clean))

if(species == "galGal6"){
clean <- clean[grep("^chr\\d+_\\d+_\\d+$|^chrZ_\\d+_\\d+$|^chrW_\\d+_\\d+$", rownames(clean), value=TRUE),]
}else if(species == "TAIR10"){
clean <- clean[grep("^chr\\d+_\\d+_\\d+$|^chrMt_\\d+_\\d+$|^chrPt_\\d+_\\d+$", rownames(clean), value=TRUE),]
}else if(species == "dm6"){
clean <- clean[grep("^chr\\w+_\\d+_\\d+$", rownames(clean), value=TRUE),]
}else{
clean <- clean[grep("^chr\\d+_\\d+_\\d+$|^chr2A_\\d+_\\d+$|^chr2B_\\d+_\\d+$|^chrX_\\d+_\\d+$|^chrY_\\d+_\\d+$", rownames(clean), value=TRUE),]
}
#clean <- clean[grep("^chr\\d+_\\d+_\\d+$", rownames(clean), value=TRUE),]

colsums <- Matrix::colSums(clean)
#write.table(colsums,"colsumsori.txt",col.names=FALSE,sep="\t",quote=FALSE)
colsumsnew <- sort(colsums,decreasing = TRUE)
extractcols <- array()
if(length(colsumsnew)>112){
extractcols <- colsumsnew[-(1:as.integer(length(colsumsnew)*0.1))]
extractcols <- extractcols[1:100]
}else{
extractcols <- colsums
}
write.table(as.matrix(clean[,names(extractcols)]),paste0("hundred/hundred_",input,".txt",sep=""),col.names=FALSE,sep="\t",quote=FALSE)
#write.table(Matrix::rowSums(clean),paste0("peak/",dim(clean)[2],"_",input,".txt",sep=""),col.names=FALSE,sep="\t",quote=FALSE)
write.table(Matrix::rowSums(clean),paste0("peak/",dim(clean)[2],"_", round(mean(colsums), digits = 0),"_",input,".txt",sep=""),col.names=FALSE,sep="\t",quote=FALSE)
#saveRDS(clean, file="/data/gts/dbscATAC/project/GSE179705/example2.rds")

if(dim(clean)[2] >= 2500){
clean <- clean[,1:2500]
clean[clean > 0] <- 1
}

if(file.exists(output_file_c)=="FALSE"){
j <- jaccardMatrix(clean)
row.sums <- apply(j, 2, sum) - 1
sum <- sum(row.sums)
weighted <- row.sums/sum
score <- t(t(clean)*weighted)
#col.sums <- apply(score, 1, sum)
rowsums <- Matrix::rowSums(score)
sort <- sort(rowsums)
#write.table(sort,"sort.txt")

random1 <- clean
if(dim(clean)[2] > 2500){
for (col in 1:ncol(clean)){
random1[,col] <- sample(clean[,col])
} 
}else{
random1 <- apply(clean,2,sample)
}
score1 <- t(t(random1)*weighted)
rowsums1 <- Matrix::rowSums(score1)
sort1 <- sort(rowsums1)

random2 <- clean
if(dim(clean)[2] > 2500){
for (col in 1:ncol(clean)){
random2[,col] <- sample(clean[,col])
} 
}else{
random2 <- apply(clean,2,sample)
}
score2 <- t(t(random2)*weighted)
rowsums2 <- Matrix::rowSums(score2)
sort2 <- sort(rowsums2)


random3 <- clean
if(dim(clean)[2] > 2500){
for (col in 1:ncol(clean)){
random3[,col] <- sample(clean[,col])
} 
}else{
random3 <- apply(clean,2,sample)
}
score3 <- t(t(random3)*weighted)
rowsums3 <- Matrix::rowSums(score3)
sort3 <- sort(rowsums3)


value <- length(sort) * 0.95
cutoff <- (sort1[value]+sort2[value]+sort3[value])/3
result <- subset(sort,sort > cutoff)

#write the jaccard index results
write.table(result,output_file_j)

# format cell info
cellist <- colnames(clean)
if(length(cellist) == length(unique(cellist))){
cellist <- as.data.frame(cellist)
rownames(cellist) <- cellist$cellist
}else{
cellist <- as.data.frame(cellist)
rownames(cellist) <- unique(replicate(length(cellist$cellist)+1000, paste0(sample(letters, 10), collapse = "")))[1:length(cellist$cellist)]
colnames(clean) <- rownames(cellist)
}
names(cellist) <- "cells"

# format peak info
peaklist <- rownames(clean)
peaklist <- as.data.frame(peaklist)
peaklist <- tidyr::separate(data = peaklist, col = peaklist, into = c("chr", "bp1", "bp2"))
peaklist$site_name <- paste(peaklist$chr, peaklist$bp1, peaklist$bp2, sep="_")
row.names(peaklist) <- peaklist$site_name



# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(clean,
                                                 cell_metadata = cellist,
                                                 gene_metadata = peaklist))

#saveRDS(input_cds, file="/data/gts/dbscATAC/project/GSE179705/exa2.rds")

input_cds <- monocle3::detect_genes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

set.seed(2023)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
if(dim(clean)[2]>150){
input_cds <- preprocess_cds(input_cds, method = "LSI")
}else{
input_cds <- preprocess_cds(input_cds, num_dim = 30, method = "LSI")
}
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")

#plot_cells(input_cds)

umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds<- NULL
if(dim(clean)[2]>150){
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
}else{
cicero_cds <- make_cicero_cds(input_cds, k = 25, reduced_coordinates = umap_coords)
}


genomepath <- ""
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
if(species == "susScr11"){
#3.0 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/susScr11.chrom.sizes"
dist_const <- 240000
wind1 <- 480000
wind2 <- 100000
}else if(species == "hg19"){
#3.1 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/hg19.chrom.sizes"
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
}else if(species == "hg38"){
#3.3 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/hg38.chrom.sizes"
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
}else if(species == "galGal6"){
#1.2 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/galGal6.chrom.sizes"
dist_const <- 100000
wind1 <- 200000
wind2 <- 40000
}else if(species == "calJac3"){
#2.87 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/calJac3.chrom.sizes"
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
}else if(species == "macFas6"){
#2.91 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/macFas6.chrom.sizes"
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
}else if(species == "rheMac10"){
#2.85 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/rheMac10.chrom.sizes"
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
}else if(species == "panTro5"){
#3 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/panTro5.chrom.sizes"
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
}else if(species == "mm9"){
#2.73 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/mm9.chrom.sizes"
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
}else if(species == "mm10"){
#2.7 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/mm10.chrom.sizes"
dist_const <- 250000
wind1 <- 500000
wind2 <- 100000
}else if(species == "danRer10"){
#1.7 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/danRer10.chrom.sizes"
dist_const <- 140000
wind1 <- 280000
wind2 <- 56000
}else if(species == "dm6"){
#0.144 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/dm6.chrom.sizes"
dist_const <- 12500
wind1 <- 25000
wind2 <- 5000
}else if(species == "TAIR10"){
#0.135 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/TAIR10.chrom.sizes"
dist_const <- 13500
wind1 <- 25000
wind2 <- 5000
}else if(species == "IRGSP1"){
#0.43 ¡Á 109 base-pairs
genomepath <- "/data/gts/dbscATAC/project/chromsizes/IRGSP1.chrom.sizes"
dist_const <- 40000
wind1 <- 70000
wind2 <- 14000
}

distance_parameters <- estimate_distance_parameter(cicero_cds,
                                                   sample_num=100,
                                                   genomic_coords = genomepath,
                                                   s = 0.85,
                                                   distance_constraint = dist_const,
                                                   window = wind1)

model_output <- generate_cicero_models(cicero_cds,
                                       distance_parameter = mean(distance_parameters),
                                       genomic_coords = genomepath,
                                       window = wind2,
                                       s = 0.85)

cicero_cons <- assemble_connections(model_output)


cicero_result <- filter(cicero_cons,coaccess>=0.1)

write.csv(cicero_result,output_file_c)
}
end_time <- Sys.time()
end_time - start_time



