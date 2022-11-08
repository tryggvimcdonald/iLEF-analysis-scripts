#!/bin/bash 
#SBATCH --job-name=bowtie2-build_AnoSag2.1_Y
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

module load Bowtie2/2.4.1-GCC-8.3.0

echo -e "\n** Script started on `date` **"

time bowtie2-build -f AnoSag2.1_PGA_Y.fa AnoSag2.1_Y

echo -e "\n** Script ended on `date` **"
echo Done!
