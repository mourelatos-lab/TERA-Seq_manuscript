#!/bin/env bash
#
#
#

source ../PARAMS.sh

assembly="hg38"
threads=10

SPANUP=0 # First position to plot from the mRNA start
SPANDOWN=150 # Last position to plot from the mRNA start

dist_meta_dir="../metacoord_correction/results" # metacoord_correction directory with annotation correction

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> DENSITY OF READ 5p AROUND MRNA 5p FOR READS WITH REL5 (TRANSCRIPTOME) - DEFAULT ANNOT <<<"

conda activate teraseq
source $INSTALL/perl-virtualenv/teraseq/bin/activate

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
)

for i in "${samples[@]}"; do
	echo " Working for" $i;

	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Uncorrected primary
	bin/dist-from-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--bed-pos 5p \
		--pos 5p \
		--span-up $SPANUP \
		--span-down $SPANDOWN \
		--where "$MATE1 AND $PRIMARY AND rel5 IS NOT NULL AND $CODINGTRANSCRIPT" \
	| paste-table-col --ifile - --col-name sample --col-val $i \
	| src/R/rel-pos-zscore-heatmap.R \
		--ifile stdin \
		--figfile $sdir/rel-pos-5p-5p-zscore-heatmap.mrna.rel5.uncorr-primary.pdf \
		--posMin "-$SPANUP" --posMax "$SPANDOWN" \
		--subset 0.3 \
		--name $i-5p-5p.mrna.rel5 &

	# Uncorrected unambiguous
	bin/dist-from-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--bed-pos 5p \
		--pos 5p \
		--span-up $SPANUP \
		--span-down $SPANDOWN \
		--where "$MATE1 AND $UNAMBIG AND rel5 IS NOT NULL AND $CODINGTRANSCRIPT" \
	| paste-table-col --ifile - --col-name sample --col-val $i \
	| src/R/rel-pos-zscore-heatmap.R \
		--ifile stdin \
		--figfile $sdir/rel-pos-5p-5p-zscore-heatmap.mrna.rel5.uncorr-unambig.pdf \
		--posMin "-$SPANUP" --posMax "$SPANDOWN" \
		--subset 0.3 \
		--name $i-5p-5p.mrna.rel5 &
done
wait

echo ">>> DENSITY OF READ 5p AROUND MRNA 5p FOR READS WITH REL5 (TRANSCRIPTOME) - CORRECTED ANNOT <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;

	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Corrected by primary
	bin/dist-from-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $dist_meta_dir/coords-corrected/genic_elements.mrna.primary.tab \
		--bed-pos 5p \
		--pos 5p \
		--span-up $SPANUP \
		--span-down $SPANDOWN \
		--where "$MATE1 AND $PRIMARY AND rel5 IS NOT NULL AND $CODINGTRANSCRIPT" \
	| paste-table-col --ifile - --col-name sample --col-val $i \
	| src/R/rel-pos-zscore-heatmap.R \
		--ifile stdin \
		--figfile $sdir/rel-pos-5p-5p-zscore-heatmap.mrna.rel5.corr-primary.pdf \
		--posMin "-$SPANUP" --posMax "$SPANDOWN" \
		--subset 0.3 \
		--name $i-5p-5p.mrna.rel5 &

	# Corrected by unambig
	bin/dist-from-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $dist_meta_dir/coords-corrected/genic_elements.mrna.unambig.tab \
		--bed-pos 5p \
		--pos 5p \
		--span-up $SPANUP \
		--span-down $SPANDOWN \
		--where "$MATE1 AND $UNAMBIG AND rel5 IS NOT NULL AND $CODINGTRANSCRIPT" \
	| paste-table-col --ifile - --col-name sample --col-val $i \
	| src/R/rel-pos-zscore-heatmap.R \
		--ifile stdin \
		--figfile $sdir/rel-pos-5p-5p-zscore-heatmap.mrna.rel5.corr-unambig.pdf \
		--posMin "-$SPANUP" --posMax "$SPANDOWN" \
		--subset 0.3 \
		--name $i-5p-5p.mrna.rel5 &
done
wait

echo ">>> COUNT TOTAL AND FULL-CDS READS - DEFAULT <<<"

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
)

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Primary
	bin/count-reads-by-decay-level \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--cds $DATA_DIR/$assembly/genic_elements.cds.bed \
		--where "$MATE1 AND $PRIMARY AND rel5 IS NOT NULL AND $CODINGTRANSCRIPT" \
	| table-paste-col --table - --col-name sample --col-val $i \
		> $sdir/counts.uncorrected-mrna.primary.rel5.tab

	src/R/compare-counts-for-decay-levels.R \
		--ifile $sdir/counts.uncorrected-mrna.primary.rel5.tab \
		--figfile $sdir/counts.uncorrected-mrna.primary.rel5.pdf \
		> $sdir/metrics.uncorrected-mrna.primary.rel5.tab

	# Unambiguous
	bin/count-reads-by-decay-level \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--cds $DATA_DIR/$assembly/genic_elements.cds.bed \
		--where "$MATE1 AND $UNAMBIG AND rel5 IS NOT NULL AND $CODINGTRANSCRIPT" \
	| table-paste-col --table - --col-name sample --col-val $i \
		> $sdir/counts.uncorrected-mrna.unambig.rel5.tab

	src/R/compare-counts-for-decay-levels.R \
		--ifile $sdir/counts.uncorrected-mrna.unambig.rel5.tab \
		--figfile $sdir/counts.uncorrected-mrna.unambig.rel5.pdf \
		> $sdir/metrics.uncorrected-mrna.unambig.rel5.tab	

done
wait

echo ">>> COUNT TOTAL AND FULL-CDS READS - CORRECTED <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Primary
	bin/count-reads-by-decay-level \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $dist_meta_dir/coords-corrected/genic_elements.mrna.primary.tab \
		--cds $DATA_DIR/$assembly/genic_elements.cds.bed \
		--where "$MATE1 AND $PRIMARY AND rel5 IS NOT NULL AND $CODINGTRANSCRIPT" \
	| table-paste-col --table - --col-name sample --col-val $i \
		> $sdir/counts.corrected-mrna.primary.rel5.tab

	src/R/compare-counts-for-decay-levels.R \
		--ifile $sdir/counts.corrected-mrna.primary.rel5.tab \
		--figfile $sdir/counts.corrected-mrna.primary.rel5.pdf \
		> $sdir/metrics.corrected-mrna.primary.rel5.tab

	# Unambiguous
	bin/count-reads-by-decay-level \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $dist_meta_dir/coords-corrected/genic_elements.mrna.unambig.tab \
		--cds $DATA_DIR/$assembly/genic_elements.cds.bed \
		--where "$MATE1 AND $PRIMARY AND rel5 IS NOT NULL AND $CODINGTRANSCRIPT" \
	| table-paste-col --table - --col-name sample --col-val $i \
		> $sdir/counts.corrected-mrna.unambig.rel5.tab

	src/R/compare-counts-for-decay-levels.R \
		--ifile $sdir/counts.corrected-mrna.unambig.rel5.tab \
		--figfile $sdir/counts.corrected-mrna.unambig.rel5.pdf \
		> $sdir/metrics.corrected-mrna.unambig.rel5.tab	

done
wait

echo ">>> ALL DONE <<<"