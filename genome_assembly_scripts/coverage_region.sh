#!/bin/bash 
#SBATCH --job-name=coverage_region
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --export=NONE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=tpm84680@uga.edu 
#SBATCH --mail-type=ALL 
#SBATCH --array=1-9

cd /scratch/tpm84680/workDir/genome_assembly/

ml SAMtools/1.9-GCC-8.3.0

i=$SLURM_ARRAY_TASK_ID

regions=("scaffold_3:1-10000000" "scaffold_3:20000000-62340000" "scaffold_3:62340001-281461234" "scaffold_6:1-115420000" "scaffold_6:115420001-135963048" "scaffold_7:1-20850000" "scaffold_7:20850001-91780800" "scaffold_7:91780801-96876404" "scaffold_7:96876405-110948685")

echo -e "** Script begun on `date` **"

time samtools depth -a -r "${regions[$i-1]}" iLEF_2T1_COMB_WGS.bam | awk -v region="${regions[$i-1]}" '{sum+=$3; sumsq+=$3*$3} END { print "Average for ",region," = ",sum/NR; print "Stdev for ",region," = ",sqrt(sumsq/NR - (sum/NR)**2)}' >> iLEF_2T1_COMB_WGS_coverage_interest_regions_summed.tsv

echo -e "** Script finished on `date` **"
