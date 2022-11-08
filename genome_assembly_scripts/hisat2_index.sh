#!/bin/bash 
#SBATCH --job-name=hisat_index_anosag2.1_Y
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

cd $SLURM_SUBMIT_DIR

ml HISAT2/2.2.1-foss-2019b

hisat2-build -f -p 6 AnoSag2.1_PGA_Y.fa AnoSag2.1_Y