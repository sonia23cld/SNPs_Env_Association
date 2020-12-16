#!/usr/bin/env bash

# SLURM
#SBATCH --mem=5GB
#SBATCH --output=/users/sonia.celestini/SNPs_Analysis/Logs/MantelTest_%A_%a.log
#SBATCH --array=1-3000
#SBATCH --partition=c
#SBATCH --qos=medium
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=20

i=${SLURM_ARRAY_TASK_ID}

while read p; do
   p+=($p) 
done <All_combinations.txt

gene=${p[i*3-1]}

n_snps=( $(cat 'Snps_parallel.txt') )
n_snps=${n_snps[gene-1]}

mkdir ~/SNPs_Analysis/Results/${gene}

module load r/3.6.2-foss-2018b
module load parallel/20190222-gcccore-7.3.0

parallel --delay 0.2 --jobs 20 srun -n1 --exclusive Rscript Mantel_test.R ::: $(seq 1 $n_snps)

n_files=( $(ls ~/SNPs_Analysis/Results/${gene} | wc -l) )
should_be=$(( $n_snps*67 )) 

if [ $n_files -eq $should_be ]
then 
	cat ~/SNPs_Analysis/Results/${gene}/MR_* >  ~/SNPs_Analysis/Results/${gene}.txt
	rm -rf ~/SNPs_Analysis/Results/${gene}
fi
