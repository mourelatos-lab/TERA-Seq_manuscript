#!/bin/bash
#
# Comparison of polyA and degradation
#
# IMPORTANT: This analysis can be run only after the metacoord_correction and polya analyses are finished
#

source ../PARAMS.sh

threads=10
assembly="hg38"

distro_res_dir="$DIR/metacoord_correction/results" # Directory with metacoordinates & correction results
polya_resdir="$DIR/polya/results" # Directory with poly(A) length results

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> GET DEGRADATION DATA - LINKS <<<"

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
)

for i in "${samples[@]}"; do
	echo $i
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	for a in $distro_res_dir/$i/coords-on-corrected-meta-mrnas*tab
	do
		ln -s $a $sdir/$(basename $a)
	done
done

echo ">>> APPEND POLYA LENGTH LIBRARIES - REL5 <<<"

source $INSTALL/perl-virtualenv/teraseq/bin/activate

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
)

for i in "${samples[@]}"; do
	echo $i
	sdir=$RES_DIR/$i

	echo ">> ON CORRECTED ANNOTATIONS <<"
	table-join --left --table1 $sdir/coords-on-corrected-meta-mrnas.unambig.rel5-inclLen.tab --key1 qname \
		--table2 $polya_resdir/$i/rna_tails.rel5.tsv --key2 readname --key-as qname \
		> $sdir/coords-on-corrected-meta-mrnas.unambig.rel5-inclLen-inclPolya.tab &
	table-join --left --table1 $sdir/coords-on-corrected-meta-mrnas.unambig.norel5-inclLen.tab --key1 qname \
		--table2 $polya_resdir/$i/rna_tails.norel5.tsv --key2 readname --key-as qname \
		> $sdir/coords-on-corrected-meta-mrnas.unambig.norel5-inclLen-inclPolya.tab &
done
wait

echo ">>> COMPARE POLYA AND DEGRADATION <<<"

deactivate
conda activate teraseq

echo ">> VIOLIN PLOT AND CORRELATION <<"

for i in "${samples[@]}"; do
	echo $i
	sdir=$RES_DIR/$i

	src/R/degradation-polya-simple.R $sdir/coords-on-corrected-meta-mrnas.unambig.rel5-inclLen-inclPolya.tab $sdir &
done
wait

echo ">> INDIVIDUAL MOLECULES WITH POLYA <<"
# Slightly different visualization and only for selected transcripts (transcripts of genes which have only one protein-coding transcript "child")

for i in "${samples[@]}"; do
	echo $i
	sdir=$RES_DIR/$i

	src/R/degradation-vs-polya.R --ifile $sdir/coords-on-corrected-meta-mrnas.unambig.rel5-inclLen-inclPolya.tab \
		--tx_list data/HeLa_transcripts.all-perGene-one2one.transcript-gene.txt \
		--outdir $sdir > ${sdir}/degradation-vs-polya.log &
done
wait

echo ">>> ALL DONE <<<"
