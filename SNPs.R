#script to start SNPs investigation: association with environmental data. You will use this data as input on cluster.
#THE LOGIC IS THIS: you have to perform (on cluster) a mantel test for every gene's SNP with one (of 68) environmental variables.
#What do you need to have to run on cluster: 
   #ForParallel.RData: contains chr_regions, climate_var, GxT_coordinates, K. 
   # All_combinations.txt: it's a table with indexes, gene number and number of env variables. It tells you how many and which combinations you have for every gene with each env variable. The job array index will refer to the index of this table (so for every index you have on gene and one env var)
   # Snps_parallel.txt: number of snps for every gene

#load SNPs and kinship matrix (rhdf5)
library(rhdf5)
K_hdf <- H5Fopen("/Users/sonia.celestini/Desktop/WGCNA_8acn_before/kinship_ibs_binary_mac5.h5py") 
K <- K_hdf$kinship
colnames(K) <- K_hdf$accessions
row.names(K) <- K_hdf$accessions
G_hdf <- h5ls('/Users/sonia.celestini/Desktop/WGCNA_8acn_before/all_chromosomes_binary_Rcompatible.hdf5')

#find coordintes for every gene (only genes GxT)
category_genes<- read.table('/Users/sonia.celestini/Desktop/WGCNA_8acn_after/Category_genes.txt', header = T)
#install package 'ape'
Araport11<- read.gff('/Users/sonia.celestini/Desktop/WGCNA_8acn_before/Araport11_GFF3_genes_transposons.201606.gff')
##clean the table: select only genes and divide the column with the name (stringr and dplyr)
Araport11<- Araport11[grep('gene', Araport11$type),]
Araport11<- cbind(Araport11[,1:8], str_split_fixed(Araport11$attributes, ';', 2))
colnames(Araport11)[9]<- 'ID'
Araport11$ID<- str_remove(Araport11$ID, 'ID=')
GxT_coordinates<- list()
GxT_coordinates$genes<- category_genes$genes[category_genes$category %in% 'GxT']
GxT_coordinates$coordinates<- lapply(GxT_coordinates$genes, function(x) {Araport11[,4:5][Araport11$ID %in% x,]})
names(GxT_coordinates$coordinates)<- GxT_coordinates$genes
GxT_coordinates$chromosome<- lapply(GxT_coordinates$genes, function(x) {as.character(Araport11[,1][Araport11$ID %in% x])})
GxT_coordinates$chromosome<- lapply(GxT_coordinates$chromosome, function(x) {as.numeric(str_remove(x, 'Chr'))})
names(GxT_coordinates$chromosome)<- GxT_coordinates$genes


climate_var<- read.table('/Users/sonia.celestini/Desktop/WGCNA_8acn_before/List_2029_latlon_env.txt', header = T, quote = '')

save(chr_regions, climate_var, GxT_coordinates, K, file = '/Users/sonia.celestini/Desktop/WGCNA_8acn_after/SNPs_Analysis/ForParallel.RData')

#create data input for new_parallel (all combinations gene/env)  
install.packages("rlist")
library(rlist)
array_combinations<- list()
for (gene in 1:length(GxT_coordinates$genes)) {
  for (env in 1:ncol(climate_var[5:71])) {
    array_combinations<-list.append(array_combinations, c(gene, env))
  }
}
#create a dataframe to use in jobscript
all_combinations<-data.frame('N'=numeric(), 'gene'=numeric(), 'env'=numeric())
for (i in 1:length(array_combinations)) {
  gene<- array_combinations[[i]][1]
  env<- array_combinations[[i]][2]
  all_combinations<-rbind(all_combinations, c(i, gene, env))
}
write.table(all_combinations, file = '/Users/sonia.celestini/Desktop/WGCNA_8acn_after/SNPs_Analysis/All_combinations.txt', sep = '\t', row.names = F, col.names = F)

#get number of snps for every gene + promoter regions
chr_regions<- h5read('/Users/sonia.celestini/Desktop/WGCNA_8acn_before/all_chromosomes_binary_Rcompatible.hdf5', '/chr_regions')

Snps_parallel<- c()
library(data.table)
for (gene in 1:length(GxT_coordinates$genes)) {
  #for every gene I need to have the indexes and then positions of its chromosom
  gene<- GxT_coordinates$genes[[gene]]
  Chr<- GxT_coordinates$chromosom[[gene]]
  coordinates_gene<- GxT_coordinates$coordinates[[gene]]
  SNPsmatrix<- gene_SNPmatrix(chr_regions, Chr, coordinates_gene)
  SNPsmatrix<- SNPsmatrix[row.names(SNPsmatrix) %in% climate_var$tg_ecotypeid,]
  SNPsmatrix<- MAF(SNPsmatrix)
  Snps_parallel<- c(Snps_parallel, ncol(SNPsmatrix))
}
write.table(Snps_parallel, file = '/Users/sonia.celestini/Desktop/WGCNA_8acn_after/SNPs_Analysis/Snps_parallel.txt', sep = '\t', row.names = F, col.names = F)

#clean table for acn that you have
climate_var<- climate_var[climate_var$tg_ecotypeid %in% rownames(SNPsmatrix),]

library(ecodist)
#create distance matrixes (order row in the same way)
climate_var<- climate_var[match(rownames(SNPsmatrix), climate_var$tg_ecotypeid),]
#colnames and rownames of K already match with the matrix, the numb of acn can not change
#trial to test if it's working before running on cluster
mantelResults <- data.frame('Name gene'= NA, 'coordinate'= NA, 'environment' = NA, 'mantelr' = NA, 'pval1' = NA, 'pval2' = NA, 'pval3' = NA, 'llim.2.5%' = NA, 'ulim.97.5%' = NA)
for (env in 5:71) {
  climate.d <- dist(climate_var[,env])
  for (snp in 1:ncol(SNPsmatrix)) {
    SNP.d <- dist(as.data.frame(SNPsmatrix)[,snp])
    m <- mantel(climate.d ~ SNP.d + as.dist(K))
    mantelResults <- rbind(mantelResults, c(gene, colnames(SNPsmatrix)[snp], colnames(climate_var)[env], m))
  }    
}
mantelResults <- na.omit(mantelResults)



