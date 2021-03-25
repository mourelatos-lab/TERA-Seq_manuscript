#!/bin/bash
#
# Calculate overlap with annotation (CAGE/APA) and plot the ouputs
#

source ../PARAMS.sh

threads=10
assembly="hg38"

dist_cage=10 # Max distance of read len to signal
dist_apa=$dist_cage

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> EXTRACT BAM COORDINATES <<<"

conda activate teraseq

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
)

for i in "${samples[@]}"; do
	echo "Working on $i"
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Make BAM to BED and overlap with protein-coding transcripts only
	# 	Has to be the best hit (XP:i:1) and have only one best hit (XN:i:1) (unambiguous)
	samtools view -@ $threads -h -F 4 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.toGenome.sorted.bam \
		| grep -e '^@HD\|^@PG\|^@SQ\|[[:space:]]XP:i:1[[:space:]]\|[[:space:]]XP:i:1$' \
		| grep -e '^@HD\|^@PG\|^@SQ\|[[:space:]]XN:i:1[[:space:]]\|[[:space:]]XN:i:1$' \
		| samtools view -@ $threads -bh - \
		| bedtools bamtobed -ed -i stdin | sort -k1,1 -k2,2n \
		| bedtools intersect -wa -u -s \
		  -a stdin -b $DATA_DIR/$assembly/ensembl_transcripts_protein_coding.gtf > $sdir/reads.1.sanitize.toGenome.sorted.bed
done

echo ">>> GET OVERLAPS - CAGE <<<"

source $INSTALL/perl-virtualenv/teraseq/bin/activate

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
)

for i in "${samples[@]}"; do
	echo "Working on $i"
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Extract only 5' ends - sample
	cat $sdir/reads.1.sanitize.toGenome.sorted.bed | grep -P "\t\+$" | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$2+1,$4,$5,$6}' \
        > $sdir/tmp.plus.bed &
	cat $sdir/reads.1.sanitize.toGenome.sorted.bed | grep -P "\t\-$" | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$3-1,$3,$4,$5,$6}' \
        > $sdir/tmp.minus.bed
	wait

	cat $sdir/tmp.plus.bed $sdir/tmp.minus.bed | sort -k1,1 -k2,2n > $sdir/reads.1.sanitize.toGenome.sorted.5p.bed
	rm $sdir/tmp.plus.bed $sdir/tmp.minus.bed

	# This is the most likely overlap - reporting the first upstream CAGE signal - sample
	bedtools closest -s -D a -id -t first -mdb all \
		-names fantom5.rep1 fantom5.rep2 fantom5.rep3 \
		-a $sdir/reads.1.sanitize.toGenome.sorted.5p.bed \
		-b $DATA_DIR/fantom5/HeLa.rep1.hg38.ctss.bed $DATA_DIR/fantom5/HeLa.rep2.hg38.ctss.bed $DATA_DIR/fantom5/HeLa.rep3.hg38.ctss.bed \
		| paste-table-col --ifile - \
		--col-val $i.5p \
		> $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.bed &
	wait
done

echo ">>> GENERATE CONTROLS - RANDOM <<<"

sdir=$RES_DIR/common
mkdir -p $sdir

cat $DATA_DIR/$assembly/genes.gtf | grep -P "\texon\t" \
    | grep -w -f data/HeLa_transcripts.all.txt > $sdir/genes.exon.sub.gtf

lines=10000000
bedtools random -l 1 -n $lines -seed $RANDOM -g $DATA_DIR/$assembly/genome/genome.fa.fai \
    | bedtools intersect -wa -u -s -f 1 \
    -a stdin -b $sdir/genes.exon.sub.gtf \
    | sort -k1,1 -k2,2n > $sdir/genes.exon.sub.random.bed

echo ">> GET OVERLAPS - RANDOM <<"

bedtools closest -s -D a -id -t first -mdb all \
    -names fantom5.rep1 fantom5.rep2 fantom5.rep3 \
    -a $sdir/genes.exon.sub.random.bed \
    -b $DATA_DIR/$assembly/fantom5/HeLa.rep1.hg38.ctss.bed $DATA_DIR/$assembly/fantom5/HeLa.rep2.hg38.ctss.bed $DATA_DIR/$assembly/fantom5/HeLa.rep3.hg38.ctss.bed \
    | paste-table-col --ifile - \
    --col-val control-transcripts \
    > $sdir/genes.exon.sub.random.5p-cage-upstream.bed

echo ">>> PLOT OVERLAPS AND STATS - CAGE <<<"

for i in "${samples[@]}"; do
	echo "Working on $i"
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	subset-table --ifile $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.bed \
		--ilist $SAMPLE_DIR/$i/fastq/reads.1.sanitize.w_rel5.names.txt \
		--nonames \
		--column 5 \
		--ofile $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.wAdapter.bed
		sed -i "s/${i}/${i}.wAd/g" $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.wAdapter.bed

	subset-table --ifile $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.bed \
		--ilist $SAMPLE_DIR/$i/fastq/reads.1.sanitize.wo_rel5.names.txt \
		--nonames \
		--column 5 \
		--ofile $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.woAdapter.bed
		sed -i "s/${i}/${i}.woAd/g" $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.woAdapter.bed

	cat $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.wAdapter.bed $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.woAdapter.bed \
		| src/R/plot-distance.R --ifile stdin --direction up --distance $dist_cage --ofile $sdir/5p-cage-upstream.w_vs_woAdapter.pdf \
        > $sdir/5p-cage-upstream.w_vs_woAdapter.txt &

    cat $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.wAdapter.bed $RES_DIR/common/genes.exon.sub.random.5p-cage-upstream.be \
        | src/R/plot-distance.R --ifile stdin --direction up --distance $dist_cage --ofile $sdir/5p-cage-upstream.w_vs_controlAdapter.pdf \
        > $sdir/5p-cage-upstream.w_vs_controlAdapter.txt &

    cat $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.woAdapter.bed $RES_DIR/common/genes.exon.sub.random.5p-cage-upstream.be \
        | src/R/plot-distance.R --ifile stdin --direction up --distance $dist_cage --ofile $sdir/5p-cage-upstream.wo_vs_controlAdapter.pdf \
        > $sdir/5p-cage-upstream.wo_vs_controlAdapter.txt &

    # wider span
    cat $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.wAdapter.bed $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.woAdapter.bed \
        | src/R/plot-distance.R --ifile stdin --direction up --distance 500 --ofile $sdir/5p-cage-upstream.w_vs_woAdapter.500.pdf \
        > $sdir/5p-cage-upstream.w_vs_woAdapter.500.txt &

    cat $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.wAdapter.bed $RES_DIR/common/genes.exon.sub.random.5p-cage-upstream.be \
        | src/R/plot-distance.R --ifile stdin --direction up --distance 500 --ofile $sdir/5p-cage-upstream.w_vs_controlAdapter.500.pdf \
        > $sdir/5p-cage-upstream.w_vs_controlAdapter.500.txt &

    cat $sdir/reads.1.sanitize.toGenome.sorted.5p-cage-upstream.woAdapter.bed $RES_DIR/common/genes.exon.sub.random.5p-cage-upstream.be \
        | src/R/plot-distance.R --ifile stdin --direction up --distance 500 --ofile $sdir/5p-cage-upstream.wo_vs_controlAdapter.500.pdf \
        > $sdir/5p-cage-upstream.wo_vs_controlAdapter.500.txt &
    wait
done
wait

echo ">>> GET OVERLAPS - APA (ALT. POLYA SITES) <<<"
# In APA or polyA we don't want direct overlap but we assume the signal is a bit upstream from the end of the read BUT
#	the way polyasite APA are annotated requires us to consider overlaps as well. It annotates clusters which are quite
#	wide and we have to include the overlap as well

# Make a temporary subset of the polyasite db since it has more columns than cage and we want to keep the consisency
mkdir $RES_DIR/common
cat $DATA_DIR/polyasite-2.0/atlas.clusters.hg38.2-0.bed | cut -f1-6 > $RES_DIR/common/polyasite-2.0.bed

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
)

for i in "${samples[@]}"; do
	echo "Working on $i"
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Extract only 3' ends - sample
	cat $sdir/reads.1.sanitize.toGenome.sorted.bed | grep -P "\t\+$" | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$3-1,$3,$4,$5,$6}' \
        > $sdir/tmp.plus.bed &
	cat $sdir/reads.1.sanitize.toGenome.sorted.bed | grep -P "\t\-$" | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$2+1,$4,$5,$6}' \
        > $sdir/tmp.minus.bed
	wait

	cat $sdir/tmp.plus.bed $sdir/tmp.minus.bed | sort -k1,1 -k2,2n > $sdir/reads.1.sanitize.toGenome.sorted.3p.bed
	rm $sdir/tmp.plus.bed $sdir/tmp.minus.bed

	# This is the most likely overlap - reporting the first upstream apa signal
	# Note: -names x doesn't work if there is only one database
	bedtools closest -s -D a -id -t first -mdb all \
    	-a $sdir/reads.1.sanitize.toGenome.sorted.3p.bed \
    	-b $RES_DIR/common/polyasite-2.0.bed \
	| awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$6,"polyasite-2.0",$7,$8,$9,$10,$11,$12,$13}' \
	| paste-table-col --ifile - \
	--col-val $i.3p \
	> $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.bed &
	wait
done

echo ">> GET OVERLAPS - RANDOM TRANSCRIPTS <<"
# We reuse the same random positions we had for CAGE since these are random positions in transcripts

sdir=$RES_DIR/common

bedtools closest -s -D a -id -t first -mdb all \
    -a $sdir/genes.exon.sub.random.bed \
    -b $RES_DIR/common/polyasite-2.0.bed \
    | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$6,"polyasite-2.0",$7,$8,$9,$10,$11,$12,$13}' \
    | paste-table-col --ifile - \
    --col-val control-transcripts \
    > $sdir/genes.exon.sub.random.3p-apa-upstream.bed

echo ">>> PLOT OVERLAPS AND STATS - APA (ALT. POLYA SITES) <<<"

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
	"hsa.dRNASeq.HeLa.total.REL3.1"
	"hsa.dRNASeq.HeLa.total.REL3.2"
	"hsa.dRNASeq.HeLa.total.REL3.3"
)

for i in "${samples[@]}"; do
	echo "Working on $i"
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Upstream
	subset-table --ifile $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.bed \
	   --ilist $SAMPLE_DIR/$i/fastq/reads.1.sanitize.w_rel[3,5].names.txt \
	   --nonames \
	   --column 5 \
	   --ofile $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.wAdapter.bed
	sed -i "s/${i}/${i}.wAd/g" $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.wAdapter.bed

	subset-table --ifile $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.bed \
	   --ilist $SAMPLE_DIR/$i/fastq/reads.1.sanitize.wo_rel[3,5].names.txt \
	   --nonames \
	   --column 5 \
	   --ofile $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.woAdapter.bed
	sed -i "s/${i}/${i}.woAd/g" $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.woAdapter.bed

	cat $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.wAdapter.bed $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.woAdapter.bed \
		| src/R/plot-distance.R --ifile stdin --direction up --distance $dist_apa --ofile $sdir/3p-apa-upstream.w_vs_woAdapter.pdf \
        > $sdir/3p-apa-upstream.w_vs_woAdapter.txt &

    cat $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.wAdapter.bed $RES_DIR/common/genes.exon.sub.random.3p-apa-upstream.bed \
        | src/R/plot-distance.R --ifile stdin --direction up --distance $dist_apa --ofile $sdir/3p-apa-upstream.w_vs_controlAdapter.pdf \
        > $sdir/3p-apa-upstream.w_vs_controlAdapter.txt &

    cat $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.woAdapter.bed $RES_DIR/common/genes.exon.sub.random.3p-apa-upstream.bed \
        | src/R/plot-distance.R --ifile stdin --direction up --distance $dist_apa --ofile $sdir/3p-apa-upstream.wo_vs_controlAdapter.pdf \
        > $sdir/3p-apa-upstream.wo_vs_controlAdapter.txt &

    # wider span
    cat $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.wAdapter.bed $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.woAdapter.bed \
        | src/R/plot-distance.R --ifile stdin --direction up --distance 500 --ofile $sdir/3p-apa-upstream.w_vs_woAdapter.500.pdf \
        > $sdir/3p-apa-upstream.w_vs_woAdapter.500.txt &

    cat $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.wAdapter.bed $RES_DIR/common/genes.exon.sub.random.3p-apa-upstream.bed \
        | src/R/plot-distance.R --ifile stdin --direction up --distance 500 --ofile $sdir/3p-apa-upstream.w_vs_controlAdapter.500.pdf \
        > $sdir/3p-apa-upstream.w_vs_controlAdapter.500.txt &

    cat $sdir/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.woAdapter.bed $RES_DIR/common/genes.exon.sub.random.3p-apa-upstream.bed \
        | src/R/plot-distance.R --ifile stdin --direction up --distance 500 --ofile $sdir/3p-apa-upstream.wo_vs_controlAdapter.500.pdf \
        > $sdir/3p-apa-upstream.wo_vs_controlAdapter.500.txt &
	wait
done

echo ">>> ALL DONE <<<"
