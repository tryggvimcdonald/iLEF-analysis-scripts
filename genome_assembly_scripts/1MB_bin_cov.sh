#!/bin/bash 
#SBATCH --job-name=1MB_bin_cov
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL 

cd /scratch/tpm84680/workDir/genome_assembly/

module load SAMtools/1.9-GCC-8.3.0
module load BBMap/38.93-GCC-8.3.0

pileup.sh in=iLEF_2T1_COMB_WGS.bam bincov=iLEF_2T1_COMB_WGS_scaffold_coverage_per_1MB_bin.txt binsize=1000000
