#!/bin/bash 
#SBATCH --job-name=coverage_regions_bedtools_coverage
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=1-00:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL 

cd /scratch/tpm84680/workDir/genome_assembly/

ml BEDTools/2.30.0-GCC-8.3.0

echo -e "** Script begun on `date` **"

time bedtools coverage -a scaffold_calculation_list.bed -b iLEF_2T1_COMB_WGS.bam > iLEF_2T1_COMB_WGS_coverage_interest_regions_bedtools_coverage.tsv

echo -e "** Script finished on `date` **"