#!/bin/bash 
#SBATCH --job-name=temp-conversion
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

module load SAMtools/1.9-GCC-8.3.0
module load BBMap/38.93-GCC-8.3.0

echo -e "\n** Script started on `date` **"

echo -e "\n** calculate coverage per scaffold **"
time pileup.sh in=iLEF_2T1_COMB_WGS.bam out=iLEF_2T1_COMB_WGS_coverage_per_scaffold.txt bincov=iLEF_2T1_COMB_WGS_scaffold_coverage_per_5k_bin.txt binsize=5000

echo -e "\n** Script ended on `date` **"
echo Done!
