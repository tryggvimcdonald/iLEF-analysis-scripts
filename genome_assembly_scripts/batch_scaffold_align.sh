#!/bin/bash
sbatch --dependency=afterok:$(squeue --noheader --format %i --name bowtie2-build_AnoSag2.1_Y --me | cut -d_ -f1 | head -n 1) /scratch/tpm84680/workDir/genome_assembly/scaffold_align.sh