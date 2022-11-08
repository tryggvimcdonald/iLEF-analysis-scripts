#!/bin/bash 
#SBATCH --job-name=bowtie2_2T1_WGS
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

INDEX=/scratch/tpm84680/workDir/genome_assembly/AnoSag2.1_Y
F1=/scratch/tpm84680/workDir/genome_assembly/iLEF_WGS_PCR_CKDL220019162-1A_H7V3VDSX5_L4_1.fq.gz
R1=/scratch/tpm84680/workDir/genome_assembly/iLEF_WGS_PCR_CKDL220019162-1A_H7V3VDSX5_L4_2.fq.gz
F2=/scratch/tpm84680/workDir/genome_assembly/iLEF_WGS_CKDL220016286-1A_H37VJDSX5_L2_1.fq.gz
R2=/scratch/tpm84680/workDir/genome_assembly/iLEF_WGS_CKDL220016286-1A_H37VJDSX5_L2_2.fq.gz

module load Bowtie2/2.4.1-GCC-8.3.0
module load SAMtools/1.9-GCC-8.3.0
module load BBMap/38.93-GCC-8.3.0

echo -e "\n** Script started on `date` **"

time bowtie2 -p 6 \
--time \
--end-to-end \
--very-sensitive \
-x ${INDEX} \
-1 ${F1},${F2} -2 ${R1},${R2} \
-S iLEF_2T1_COMB_WGS.sam 

echo -e "\n** Convert sam to bam using samtools **"
time samtools view -Sb iLEF_2T1_COMB_WGS.sam | samtools sort - -T iLEF_2T1_COMB_WGS.sorted -o iLEF_2T1_COMB_WGS.bam 

echo -e "\n** calculate coverage per scaffold **"
time pileup.sh in=iLEF_2T1_COMB_WGS.bam out=iLEF_2T1_COMB_WGS_coverage_per_scaffold.txt bincov=iLEF_2T1_COMB_WGS_scaffold_coverage_per_5k_bin.txt binsize=5000

echo -e "\n** Script ended on `date` **"
echo Done!
