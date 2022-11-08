#!/bin/bash
sbatch --dependency=afterok:$(squeue --noheader --format %i --name bowtie2_2T1_WGS --me | cut -d_ -f1 | head -n 1) /scratch/tpm84680/workDir/genome_assembly/generate_bigwig.sh