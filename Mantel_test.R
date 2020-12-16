slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')                                                                                 
gene_n_env<-as.numeric(slurm_arrayid) 
snp<- as.numeric(commandArgs(trailingOnly = TRUE))

setwd('~/SNPs_Analysis/Scripts_Data/')
all_combinations<- read.table(file='All_combinations.txt', sep='\t') 
gene_n<- all_combinations$V2[all_combinations$V1 %in% gene_n_env]
env<- (all_combinations$V3[all_combinations$V1 %in% gene_n_env])+4
#load environment
library(ecodist)
library(rhdf5)
library(data.table)
load('ForParallel.RData')

#find snps for gene + promoter regions
gene_SNPmatrix<- function(chr_regions, Chr, coordinates_gene){
	  indexes<- chr_regions[, Chr]
  positions<- h5read('all_chromosomes_binary_Rcompatible.hdf5', '/positions', index = list((indexes[1]+1):indexes[2]))
    gene_indexes<- which(positions %in% c((coordinates_gene[[1]]-10000):(coordinates_gene[[2]]+10000)))
    gene_indexes<- gene_indexes + indexes[1]
      SNPsmatrix<- h5read('all_chromosomes_binary_Rcompatible.hdf5', '/snps', index = list(1:1135, min(gene_indexes):max(gene_indexes)))
      rownames(SNPsmatrix)<- h5read('all_chromosomes_binary_Rcompatible.hdf5', '/accessions')
        colnames(SNPsmatrix)<- gene_indexes
        return(SNPsmatrix)}

MAF<- function(SNPsmatrix) {
	  row.n<- c(rownames(SNPsmatrix))
  n0<- apply(SNPsmatrix, 2, function(x) sum(x == 0))
    allele_freq<- unlist(lapply(n0, function(x) {(x * 100)/nrow(SNPsmatrix)}))
    cut_off<- c(which(allele_freq <10), which(allele_freq >90))
      SNPsmatrix<- as.data.table(SNPsmatrix)[, (cut_off) := NULL]
      SNPsmatrix<- as.data.frame(SNPsmatrix)
        rownames(SNPsmatrix)<- row.n
        return(SNPsmatrix)}

##for every gene I need to have the indexes of the promoter (at beginning and end 10 kb - it will be in the function-) and then positions of its chromosom
gene_name<- GxT_coordinates$gene[[gene_n]]
Chr<- GxT_coordinates$chromosom[[gene_n]]
coordinates_gene<- GxT_coordinates$coordinates[[gene_n]]
SNPsmatrix<- gene_SNPmatrix(chr_regions, Chr, coordinates_gene)
SNPsmatrix<- SNPsmatrix[row.names(SNPsmatrix) %in% climate_var$tg_ecotypeid,]

#clean matrix: minor allele frequency cut-off
SNPsmatrix<- MAF(SNPsmatrix)

#clean table for acn that you have
climate_var<- climate_var[climate_var$tg_ecotypeid %in% rownames(SNPsmatrix),]

#execute mantel test(correlation between snps and temperature) 
mantelResults <- data.frame('Name gene'= NA, 'coordinate'= NA, 'environment' = NA, 'mantelr' = NA, 'pval1' = NA, 'pval2' = NA, 'pval3' = NA, 'llim.2.5%' = NA, 'ulim.97.5%' = NA)
climate.d <- dist(climate_var[,env])
SNP.d <- dist(as.data.frame(SNPsmatrix)[,snp])
m <- mantel(climate.d ~ SNP.d + as.dist(K))
mantelResults <- rbind(mantelResults, c(gene_name, colnames(SNPsmatrix)[snp], colnames(climate_var)[env], m))
mantelResults <- na.omit(mantelResults)

write.table(mantelResults, paste('~/SNPs_Analysis/Results/', gene_n, '/MR_', gene_n,'_', env,'_', snp, '.txt', sep=''), quote=F, row.names=F, col.names=T)


