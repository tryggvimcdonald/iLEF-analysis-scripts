#!/bin/bash
#SBATCH --job-name=vcf_index
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL

cd /scratch/tpm84680/workDir/genome_assembly

module load  tabix/0.2.6-GCCcore-8.3.0

echo "** Script has begun! **"

time bgzip full_SNPs_filtered.vcf

time tabix -p vcf full_SNPs_filtered.vcf.gz

echo "** Script completed successfully on `date` **"
