#!/bin/bash

tail -n +4 iLEF_2T1_COMB_WGS_scaffold_coverage_per_1MB_bin.txt | awk -F'\t' 'BEGIN {OFS=FS} {print $1, $3-1000000, $3, $2}' > iLEF_2T1_COMB_WGS_scaffold_coverage_per_1MB_bin_headless_rearranged.txt