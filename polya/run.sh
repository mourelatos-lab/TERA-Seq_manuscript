#!/bin/bash
#
# Analyza and plot poly(A) tail 
#

source ../PARAMS.sh

threads=6
assembly="hg38"

RES_DIR="results"
mkdir -p $RES_DIR

##########################################################################

echo ">>> PREPARE FILES WITH POLY-A LENGTH - NANOPOLISH <<<"

conda activate teraseq
source $INSTALL/perl-virtualenv/teraseq/bin/activate

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for i in "${samples[@]}"; do
	echo " Working for" $i
	sdir=$RES_DIR/$i

	mkdir -p $sdir

	cat $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab \
		| paste-table-col --ifile - --col-name sample --col-val $i \
		> $sdir/polya_length.tab &
done
wait

echo ">>> MAKE POLYA PLOTS - 5TERA <<<"

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for sample in ${samples[@]}; do
	echo $sample
	sdir=$RES_DIR/$sample
	mkdir -p $sdir

	## Get separately w/wo adapter and only the coding transcripts
	# wo adapter
	sqlite3 $SAMPLE_DIR/$sample/db/sqlite.db 'SELECT qname FROM transcr WHERE (rel5 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' > $sdir/reads.norel5.txt &
	# w adapter
	sqlite3 $SAMPLE_DIR/$sample/db/sqlite.db 'SELECT qname FROM transcr WHERE (rel5 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' > $sdir/reads.rel5.txt &
done
wait	

for sample in ${samples[@]}; do
	echo $sample
	sdir=$RES_DIR/$sample
	mkdir -p $sdir

	subset-table \
		--ifile $sdir/polya_length.tab \
		--column "readname" \
		--ilist $sdir/reads.norel5.txt \
		--ofile $sdir/rna_tails.norel5.tsv 
	
	sed -i "s/$sample\t/${sample}.norel5\t/g" $sdir/rna_tails.norel5.tsv

	subset-table \
		--ifile $sdir/polya_length.tab \
		--column "readname" \
		--ilist $sdir/reads.rel5.txt \
		--ofile $sdir/rna_tails.rel5.tsv 

	sed -i "s/$sample\t/${sample}.rel5\t/g" $sdir/rna_tails.rel5.tsv

	table-cat \
		$sdir/rna_tails.rel5.tsv \
		$sdir/rna_tails.norel5.tsv \
		| src/R/polya-distro.R \
			--ifile stdin \
			--software "nanopolish" \
			--maxL 300 \
			--removezero \
			--figfile $sdir/polya_length_nozero.pdf &

	table-cat \
		$sdir/rna_tails.rel5.tsv \
		$sdir/rna_tails.norel5.tsv \
		| src/R/polya-distro.R \
			--ifile stdin \
			--software "nanopolish" \
			--maxL 300 \
			--figfile $sdir/polya_length_all.pdf &
done
wait

echo ">>> MAKE POLYA PLOTS - TERA3 <<<"

samples=(
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for sample in ${samples[@]}; do
	echo $sample
	sdir=$RES_DIR/$sample
	mkdir -p $sdir

	## Get separately w/wo adapter and only the coding transcripts
	# wo adapter
	sqlite3 $SAMPLE_DIR/$sample/db/sqlite.db 'SELECT qname FROM transcr WHERE (rel3 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' > $sdir/reads.norel3.txt &
	# w adapter
	sqlite3 $SAMPLE_DIR/$sample/db/sqlite.db 'SELECT qname FROM transcr WHERE (rel3 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));' > $sdir/reads.rel3.txt &
done
wait

for sample in ${samples[@]}; do
	echo $sample
	sdir=$RES_DIR/$sample
	mkdir -p $sdir

	subset-table \
		--ifile $sdir/polya_length.tab \
		--column "readname" \
		--ilist $sdir/reads.norel3.txt \
		--ofile $sdir/rna_tails.norel3.tsv 
	
	sed -i "s/$sample\t/${sample}.norel3\t/g" $sdir/rna_tails.norel3.tsv

	subset-table \
		--ifile $sdir/polya_length.tab \
		--column "readname" \
		--ilist $sdir/reads.rel3.txt \
		--ofile $sdir/rna_tails.rel3.tsv 

	sed -i "s/$sample\t/${sample}.rel3\t/g" $sdir/rna_tails.rel3.tsv

	table-cat \
		$sdir/rna_tails.rel3.tsv \
		$sdir/rna_tails.norel3.tsv \
		| src/R/polya-distro.R \
			--ifile stdin \
			--software "nanopolish" \
			--maxL 300 \
			--removezero \
			--figfile $sdir/polya_length_nozero.pdf &

	table-cat \
		$sdir/rna_tails.rel3.tsv \
		$sdir/rna_tails.norel3.tsv \
		| src/R/polya-distro.R \
			--ifile stdin \
			--software "nanopolish" \
			--maxL 300 \
			--figfile $sdir/polya_length_all.pdf &
done
wait

echo ">>> COMPARE LIBRARIES <<<"

table-cat \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/rna_tails.rel5.tsv \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1/rna_tails.rel5.tsv \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.REL5.long.1/rna_tails.rel5.tsv \
    $RES_DIR/hsa.dRNASeq.HeLa.total.REL3.2/rna_tails.rel3.tsv \
	| src/R/polya-distro.R \
		--ifile stdin \
		--software "nanopolish" \
		--maxL 300 \
		--removezero \
		--figfile $RES_DIR/polya_length_nozero.pdf &

table-cat \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/rna_tails.rel5.tsv \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1/rna_tails.rel5.tsv \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.REL5.long.1/rna_tails.rel5.tsv \
    $RES_DIR/hsa.dRNASeq.HeLa.total.REL3.2/rna_tails.rel3.tsv \
	| src/R/polya-distro.R \
		--ifile stdin \
		--software "nanopolish" \
		--maxL 300 \
		--figfile $RES_DIR/polya_length_all.pdf &
wait

echo ">>> ALL DONE <<<"