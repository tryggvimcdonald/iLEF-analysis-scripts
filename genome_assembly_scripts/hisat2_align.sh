#!/bin/bash 
#SBATCH --job-name=hisat2_align_total
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL 

cd $SLURM_SUBMIT_DIR

ml HISAT2/2.2.1-foss-2019b
ml SAMtools/1.9-GCC-8.3.0
ml BBMap/38.93-GCC-8.3.0

INDEX=/scratch/tpm84680/workDir/genome_assembly/AnoSag2.1_Y

F1=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_DMSO_R1_F.fastq.gz
R1=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_DMSO_R1_R.fastq.gz
F2=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_DMSO_R2_F.fastq.gz
R2=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_DMSO_R2_R.fastq.gz
F3=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_DMSO_R3_F.fastq.gz
R3=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_DMSO_R3_R.fastq.gz
F4=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_DMSO_R4_F.fastq.gz
R4=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_DMSO_R4_R.fastq.gz

F5=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_SAG_R1_F.fastq.gz
R5=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_SAG_R1_R.fastq.gz
F6=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_SAG_R2_F.fastq.gz
R6=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_SAG_R2_R.fastq.gz
F7=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_SAG_R3_F.fastq.gz
R7=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_SAG_R3_R.fastq.gz
F8=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_SAG_R4_F.fastq.gz
R8=/scratch/tpm84680/workDir/genome_assembly/raw_fastq_files/2T1_SAG_R4_R.fastq.gz

F9=/scratch/tpm84680/workDir/genome_assembly/iLEF_WGS_PCR_CKDL220019162-1A_H7V3VDSX5_L4_1.fq.gz
R9=/scratch/tpm84680/workDir/genome_assembly/iLEF_WGS_PCR_CKDL220019162-1A_H7V3VDSX5_L4_2.fq.gz
F10=/scratch/tpm84680/workDir/genome_assembly/iLEF_WGS_CKDL220016286-1A_H37VJDSX5_L2_1.fq.gz
R10=/scratch/tpm84680/workDir/genome_assembly/iLEF_WGS_CKDL220016286-1A_H37VJDSX5_L2_2.fq.gz


hisat2 -p 8 \
-x ${INDEX} \
-1 ${F1},${F2},${F3},${F4},${F5},${F6},${F7},${F8},${F9},${F10} -2 ${R1},${R2},${R3},${R4},${R5},${R6},${R7},${R8},${R9},${R10} \
-S 2T1_TOT_RNA.sam

samtools view -Sb 2T1_TOT_RNA.sam | samtools sort - -T 2T1_TOT_RNA.sorted -o 2T1_TOT_RNA.bam

time pileup.sh in=2T1_TOT_RNA.bam out=2T1_TOT_RNA_coverage_per_scaffold.txt bincov=2T1_TOT_RNA_scaffold_coverage_per_5k_bin.txt binsize=5000

