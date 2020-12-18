#!/bin/bash
#
# Run visualization of Cutadapt adapter trimming
#

RES_DIR="results"

mkdir -p $RES_DIR

####################################################################################################

conda activate teraseq

./src/R/trimming.R data/trimming.tsv $RES_DIR/trimming.pdf

echo ">>> ALL DONE <<<"
