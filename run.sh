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
cd ../data/
./run.sh &> log.out

#echo ">>> POPULATING SAMPLES DIRECTORY <<<"
#cd ../samples/
#./run-populate.sh

echo ">>> PREPARING SAMPLES <<<"
cd ../samples/
./run.sh &> log.out

echo ">>> RUNNING ADAPTER LENGTH VISUALIZATION <<<"
cd ../adapter/
./run.sh &> log.out

echo ">>> RUNNING ALIGNMENT STATISTICS <<<"
cd ../align-stats/
./run.sh &> log.out

echo ">>> RUNNING ANNOTATION CORRECTION AND META-COORDINATES <<<"
cd ../metacoord_correction/
./run.sh &> log.out

echo ">>> RUNNING CHANGES AFTER CORRECTION <<<"
cd ../reannot-change/
./run.sh &> log.out

echo ">>> RUNNING TRANSCRIPT COVERAGE <<<"
cd ../trans-coverage/
./run.sh &> log.out

echo ">>> RUNNING POLY(A) LENGTH <<<"
cd ../polya/
./run.sh &> log.out

echo ">>> RUNNING CAGE AND APA <<<"
cd ../cage_apa/
./run.sh &> log.out

echo ">>> RUNNING PROMOTER HEATMAP <<<"
cd ../promoter-heatmap/
./run.sh &> log.out

echo ">>> RUNNING RELATIVE POSITION DISTRIBUTION <<<"
cd ../relative-pos-distro/
./run.sh &> log.out

echo ">>> RUNNING CONSERVATION <<<"
cd ../conservation/
./run.sh &> log.out

echo ">>> RUNNING SIRV <<<"
cd ../sirv/
./run.sh &> log.out

echo ">>> RUNNING META-COORDINATES AND POLY(A) CORRELATION <<<"
cd ../metacoord-vs-polya/
./run.sh &> log.out

echo ">>> RUNNING EXPRESSION CORRELATION <<<"
cd ../expression/
./run.sh &> log.out

echo ">>> ALL DONE <<<"
