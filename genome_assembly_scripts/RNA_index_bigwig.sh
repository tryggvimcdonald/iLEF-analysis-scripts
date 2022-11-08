#!/bin/bash 
#SBATCH --job-name=RNA-bw
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL 

cd $SLURM_SUBMIT_DIR

ml purge
ml SAMtools/1.6-foss-2019b
ml deepTools/3.3.1-intel-2019b-Python-3.7.4

echo -e "\n** Script started on `date` **"
echo -e "\n** Build index for bam files **"

time samtools index /scratch/tpm84680/workDir/genome_assembly/2T1_TOT_RNA.bam 2T1_TOT_RNA.bam.bai

echo -e "\n** Convert bam to bw **"

time bamCoverage -p 6 -b /lustre2/scratch/tpm84680/workDir/genome_assembly/2T1_TOT_RNA.bam -o 2T1_TOT_RNA_50bp.bw \
--binSize 50

time bamCoverage -p 6 -b /lustre2/scratch/tpm84680/workDir/genome_assembly/2T1_TOT_RNA.bam -o 2T1_TOT_RNA_500bp.bw \
--binSize 500

time bamCoverage -p 6 -b /lustre2/scratch/tpm84680/workDir/genome_assembly/2T1_TOT_RNA.bam -o 2T1_TOT_RNA_5000bp.bw \
--binSize 5000

echo -e "\n** Script ended on `date` **"
echo Done!
