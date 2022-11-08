#!/bin/bash
#SBATCH --job-name=snp_calling
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=120:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL

cd /scratch/tpm84680/workDir/genome_assembly

module load SAMtools/1.9-GCC-8.3.0
module load picard/2.16.0-Java-1.8.0_144
module load BCFtools/1.10.2-GCC-8.3.0

echo -e "\n** Script started on `date` **"

# Detect and remove duplicates
time java -Xmx20g -classpath "/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar /apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
	METRICS_FILE=out.metrics \
	REMOVE_DUPLICATES=true \
	ASSUME_SORTED=true  \
	VALIDATION_STRINGENCY=LENIENT \
	INPUT=2T1_TOT_RNA.bam \
	OUTPUT=2T1_TOT_RNA_no_dups.bam

# Index deduplicated sorted bam file
time samtools index 2T1_TOT_RNA_no_dups.bam

# Create the pileup 
time samtools mpileup -uf AnoSag2.1_PGA_Y.fa 2T1_TOT_RNA_no_dups.bam > full_raw.bcf

# Call SNPs (Make sure to load BCFtools/1.9-foss-2016b)
time bcftools call --keep-alts --multiallelic-caller --variants-only --output-type v full_raw.bcf -o full_SNPs.bcf

# Filter SNPs
time bcftools view full_SNPs.bcf | vcfutils.pl varFilter - > full_SNPs_filtered.vcf

echo -e "\n** Script ended on `date` **"
echo Done!