#!/bin/bash
#
# Analyze conservation around 5' end of the mapped reads
#

source ../PARAMS.sh

SPANUP=100
SPANDOWN=100
LENLIM=300
QAREA=100000

assembly="hg38"

RES_DIR="results"
mkdir -p $RES_DIR

################################################################################

echo ">>> DOWNLOAD PHASTCONS DATA <<<"

mkdir data

rsync -avz --progress \
 	rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/ data/

echo ">> PREPARE JSON WITH CONSERVATION <<"

src/Perl/gff-to-cons-json.pl \
	--gtf $DATA_DIR/$assembly/genes.gtf \
	--cons-dir data/hg38.100way.phastCons/ \
	--chr-pref \
	> data/transcripts.cons.json

echo ">>> CONSERVATION AROUND 5p (W ADAPTER) <<<"

samples=(
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
)


for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	bin/htsdb-var-distro-in-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--json data/transcripts.cons.json \
		--bed $DATA_DIR/$assembly/genic_elements.cds.bed \
		--pos 5p \
		--span-up $SPANUP \
		--span-down $SPANDOWN \
		--len-lim $LENLIM \
		--q-area $QAREA \
		--where "rel5 IS NOT NULL" \
	| paste-table-col --ifile - --col-name sample --col-val $i \
	| paste-table-col --ifile - --col-name cond --col-val rel5 \
		> $sdir/cons.5p.cds.rel5.tab &
done
wait

echo ">>> PLOT CONSERVATION AROUND 5p (W ADAPTER) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	src/R/plot-var-distro.R \
		--ifile $sdir/cons.5p.cds.rel5.tab \
		--figfile $sdir/cons.5p.cds.rel5.pdf \
		--name $i-5p-cds-rel5
done
wait

echo ">>> CONSERVATION AROUND 5p (WO ADAPTER) <<<"

for i in "${samples[@]}"; do
    echo " Working for" $i;
    sdir=$RES_DIR/$i
    mkdir -p $sdir

    bin/htsdb-var-distro-in-feats \
        --db $SAMPLE_DIR/$i/db/sqlite.db \
        --table transcr \
        --json data/transcripts.cons.json \
        --bed $DATA_DIR/$assembly/genic_elements.cds.bed \
        --pos 5p \
        --span-up $SPANUP \
        --span-down $SPANDOWN \
        --len-lim $LENLIM \
        --q-area $QAREA \
        --where "rel5 IS NULL" \
    | paste-table-col --ifile - --col-name sample --col-val $i \
    | paste-table-col --ifile - --col-name cond --col-val norel5 \
        > $sdir/cons.5p.cds.norel5.tab &
done
wait

echo ">>> PLOT CONSERVATION AROUND 5p (WO ADAPTER) <<<"

for i in "${samples[@]}"; do
    echo " Working for" $i;
    sdir=$RES_DIR/$i
    mkdir -p $sdir

    src/R/plot-var-distro.R \
        --ifile $sdir/cons.5p.cds.norel5.tab \
        --figfile $sdir/cons.5p.cds.norel5.pdf \
        --name $i-5p-cds-norel5
done
wait

echo ">>> CONSERVATION AROUND 5p (RANDOM) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	bin/htsdb-var-distro-in-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--json data/transcripts.cons.json \
		--bed $DATA_DIR/$assembly/genic_elements.cds.bed \
		--pos 5p \
		--span-up $SPANUP \
		--span-down $SPANDOWN \
		--len-lim $LENLIM \
		--q-area $QAREA \
		--random \
	| paste-table-col --ifile - --col-name sample --col-val $i \
	| paste-table-col --ifile - --col-name cond --col-val random \
		> $sdir/cons.5p.cds.random.tab &
done
wait

echo ">>> PLOT CONSERVATION AROUND 5p (RANDOM) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	src/R/plot-var-distro.R \
		--ifile $sdir/cons.5p.cds.random.tab \
		--figfile $sdir/cons.5p.cds.random.pdf \
		--name $i-5p-cds-random
done
wait

echo ">>> CONSERVATION AROUND 5p (RANDOM, KEEP CONTENT) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	bin/htsdb-var-distro-in-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--json data/transcripts.cons.json \
		--bed $DATA_DIR/$assembly/genic_elements.cds.bed \
		--pos 5p \
		--span-up $SPANUP \
		--span-down $SPANDOWN \
		--len-lim $LENLIM \
		--q-area $QAREA \
		--random \
		--fasta $DATA_DIR/$assembly/transcripts.fa \
		--win-up 1 \
		--win-down 0 \
	| paste-table-col --ifile - --col-name sample --col-val $i \
	| paste-table-col --ifile - --col-name cond --col-val keepcontent \
		> $sdir/cons.5p.cds.keepcontent.tab &
done
wait

echo ">>> PLOT CONSERVATION AROUND 5p (RANDOM, KEEP CONTENT) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	src/R/plot-var-distro.R \
		--ifile $sdir/cons.5p.cds.keepcontent.tab \
		--figfile $sdir/cons.5p.cds.keepcontent.pdf \
		--name $i-5p-cds-keepcontent
done
wait

echo ">>> AGGREGATE CONSERVATION AROUND 5p SAMPLES <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i
	rdir="$RES_DIR/groups/$i"
	mkdir -p $rdir

	table-cat \
		$RES_DIR/$i/cons.5p.cds.rel5.tab \
        $RES_DIR/$i/cons.5p.cds.norel5.tab \
		$RES_DIR/$i/cons.5p.cds.random.tab \
		$RES_DIR/$i/cons.5p.cds.keepcontent.tab \
		> $rdir/cons.5p.cds.tab

	src/R/plot-var-distro.R \
		--ifile $rdir/cons.5p.cds.tab \
		--figfile $rdir/cons.5p.cds.pdf \
		--name $i-5p
done

echo ">>> ALL DONE <<<"