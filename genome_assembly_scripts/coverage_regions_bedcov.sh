#!/bin/bash 
#SBATCH --job-name=coverage_regions_bedcov
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL 

cd /scratch/tpm84680/workDir/genome_assembly/

ml SAMtools/1.9-GCC-8.3.0

echo -e "** Script begun on `date` **"

time samtools bedcov scaffold_calculation_list.bed iLEF_2T1_COMB_WGS.bam > iLEF_2T1_COMB_WGS_coverage_interest_regions_bedcov.tsv

echo -e "** Script finished on `date` **"