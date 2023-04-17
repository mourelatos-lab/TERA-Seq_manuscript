#!/bin/bash
#
# Run script to install all dependencies, prepare all references and run all analyses
#
# You still have to make sure the sample directories are populated!
#

main_dir=$(dirname "$0") # Save location of this script to a variable

echo ">>> INSTALLING ADDITIONAL SOFTWARE <<<"
cd tools/
./run.sh &> log.out

echo ">>> COMPILING REFERENCES <<<"
cd ${main_dir}/tools
./run.sh &> log.out

#echo ">>> BASECALL FAST5 <<<"
# Note: You have to populate fast5 directories first if you want to re-basecall
#cd ${main_dir}/samples/
#./run_Guppy.sh

echo ">>> PREPARING SAMPLES <<<"
cd ${main_dir}/samples/
./run.sh &> log.out

echo ">>> RUNNING ADAPTER LENGTH VISUALIZATION <<<"
cd ${main_dir}/adapter/
./run.sh &> log.out

echo ">>> RUNNING ALIGNMENT STATISTICS <<<"
cd ${main_dir}/align-stats/
./run.sh &> log.out

echo ">>> RUNNING ANNOTATION CORRECTION AND META-COORDINATES <<<"
cd ${main_dir}/metacoord-correction/
./run.sh &> log.out

echo ">>> RUNNING CHANGES AFTER CORRECTION <<<"
cd ${main_dir}/reannot-change/
./run.sh &> log.out

echo ">>> RUNNING TRANSCRIPT COVERAGE <<<"
cd ${main_dir}/trans-coverage/
./run.sh &> log.out

echo ">>> RUNNING POLY(A) LENGTH <<<"
cd ${main_dir}/polya/
./run.sh &> log.out

echo ">>> RUNNING CAGE AND APA <<<"
cd ${main_dir}/cage_apa/
./run.sh &> log.out

echo ">>> RUNNING PROMOTER HEATMAP <<<"
cd ${main_dir}/promoter-heatmap/
./run.sh &> log.out

echo ">>> RUNNING RELATIVE POSITION DISTRIBUTION <<<"
cd ${main_dir}/relative-pos-distro/
./run.sh &> log.out

echo ">>> RUNNING CONSERVATION <<<"
cd ${main_dir}/conservation/
./run.sh &> log.out

echo ">>> RUNNING SIRV <<<"
cd ${main_dir}/sirv/
./run.sh &> log.out

echo ">>> RUNNING META-COORDINATES AND POLY(A) CORRELATION <<<"
cd ${main_dir}/metacoord-vs-polya/
./run.sh &> log.out

echo ">>> RUNNING EXPRESSION CORRELATION <<<"
cd ${main_dir}/expression/
./run.sh &> log.out

echo ">>> ALL DONE <<<"
