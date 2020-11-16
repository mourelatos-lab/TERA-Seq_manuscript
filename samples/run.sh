#!/bin/bash
#
# Preprocess all libraries - prepare fastq, align to genome and transcriptome, prepare database
#
################################################################################

# If you have fast5 files you can re-run basecalling
#echo ">>> GUPPY BASECALLING <<<"
#./run_Guppy.sh

echo ">>> PREPARE SAMPLES <<<"

echo ">> RUNNING 5TERA-SHORT <<"
./run_5TERA-short.sh
echo ">> RUNNING 5TERA <<"
./run_5TERA.sh
echo ">> RUNNING TERA3-SHORT <<"
./run_TERA3.sh
echo ">> RUNNING 5TERA3 <<"
./run_5TERA3.sh
echo ">> MERGING 5TERA3 <<"
./run_5TERA3-merge.sh
echo ">> RUNNING 5TERA-SIRV <<"
./run_5TERA-SIRV.sh
echo ">> RUNNING MOUSE SIRV <<"
./run_mouse_SIRV.sh
echo ">> RUNNING AKRON5 <<"
./run_Akron5Seq.sh
echo ">> RUNNING RIBOSEQ <<"
./run_RiboSeq.sh
echo ">> RUNNING RNASEQ <<"
./run_RNASeq.sh

echo ">>> ALL DONE <<<"
