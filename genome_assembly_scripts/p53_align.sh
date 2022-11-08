#!/bin/bash 
#SBATCH --job-name=bowtie2_p53
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL 

cd $SLURM_SUBMIT_DIR

INDEX=/scratch/tpm84680/workDir/genome_assembly/AnoSag2.1_Y
P53=/scratch/tpm84680/workDir/genome_assembly/p53.fa

module load Bowtie2/2.4.1-GCC-8.3.0
module load SAMtools/1.9-GCC-8.3.0
module load BBMap/38.93-GCC-8.3.0


bowtie2 -p 6 \
--local \
--very-sensitive \
-x ${INDEX} \
-f ${P53} \
-S iLEF_2T1_PCR_WGS_p53_alignment.sam 

samtools view -Sb iLEF_2T1_PCR_WGS_p53_alignment.sam | samtools sort - -T iLEF_2T1_PCR_WGS_p53_alignment.sorted -o iLEF_2T1_PCR_WGS_p53_alignment.bam 