#!/bin/bash
#
# Run visualization of Cutadapt adapter trimming
#

RES_DIR="results"

conda activate teraseq

mkdir -p $RES_DIR

####################################################################################################

./src/R/trimming.R data/trimming.tsv $RES_DIR/trimming.pdf

echo ">>> ALL DONE <<<"