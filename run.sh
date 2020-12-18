#!/bin/bash
#
# Run script to install all dependencies, prepare all references and run all analyses
#
# You still have to make sure the sample directories are populated!
#

main_dir = $(dirname "$0") # Save location of this script to a variable

cd $main_dir/

echo ">>> INSTALLING ADDITIONAL SOFTWARE <<<"
cd tools/
./run.sh

echo ">>> COMPILING REFERENCES <<<"
cd ../data/
./run.sh

#echo ">>> POPULATING SAMPLES DIRECTORY <<<"
#cd ../samples/
#./run-populate.sh

echo ">>> PREPARING SAMPLES <<<"
cd ../samples/
./run.sh

echo ">>> RUNNING ADAPTER LENGTH VISUALIZATION <<<"
cd ../adapter/
./run.sh

echo ">>> RUNNING ALIGNMENT STATISTICS <<<"
cd ../align-stats/
./run.sh

echo ">>> RUNNING ANNOTATION CORRECTION AND META-COORDINATES <<<"
cd ../metacoord_correction/
./run.sh

echo ">>> RUNNING CHANGES AFTER CORRECTION <<<"
cd ../reannot-change/
./run.sh

echo ">>> RUNNING TRANSCRIPT COVERAGE <<<"
cd ../trans-coverage/
./run.sh

echo ">>> RUNNING POLY(A) LENGTH <<<"
cd ../polya/
./run.sh

echo ">>> RUNNING CAGE AND APA <<<"
cd ../cage_apa/
./run.sh

echo ">>> RUNNING RELATIVE POSITION DISTRIBUTION <<<"
cd ../relative-pos-distro/
./run.sh

echo ">>> RUNNING SIRV <<<"
cd ../sirv/
./run.sh

echo ">>> RUNNING META-COORDINATES AND POLY(A) CORRELATION <<<"
cd ../metacoord-vs-polya/
./run.sh

echo ">>> RUNNING EXPRESSION CORRELATION <<<"
cd ../expression/
./run.sh

echo ">>> ALL DONE <<<"
