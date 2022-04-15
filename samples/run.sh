#!/bin/bash
#
# Preprocess all libraries - prepare fastq, align to genome and transcriptome, prepare database
#
################################################################################

# If you have fast5 files you can re-basecall
#echo ">>> GUPPY BASECALLING <<<"
#./run_Guppy.sh

# If you want to run poly(A) tail estimates, you must download fast5 files and/or re-basecall them.

echo ">>> SAMPLES PROCESSING <<<"

echo ">> RUNNING 5TERA-SHORT <<"
./run_5TERA-short.sh &> run_5TERA-short.out
echo ">> RUNNING 5TERA <<"
./run_5TERA.sh &> run_5TERA.out
echo ">> RUNNING TERA3-SHORT <<"
./run_TERA3.sh &> run_TERA3.out
echo ">> RUNNING 5TERA3 <<"
./run_5TERA3.sh &> run_5TERA3.out
echo ">> MERGING 5TERA3 <<"
./run_5TERA3-merge.sh &> run_5TERA3-merge.out
echo ">> RUN dRNA (NO ADAPT) <<"
./run_dRNASeq.sh &> run_dRNASeq.out
echo ">> RUNNING MOUSE SIRV <<"
./run_mouse_SIRV.sh &> run_mouse_SIRV.out
echo ">> RUNNING AKRON5 <<"
./run_Akron5Seq.sh &> run_Akron5Seq.out
echo ">> RUNNING RIBOSEQ <<"
./run_RiboSeq.sh &> run_RiboSeq.out
echo ">> RUNNING RNASEQ <<"
./run_RNASeq.sh &> run_RNASeq.out

echo ">>> ALL DONE <<<"
