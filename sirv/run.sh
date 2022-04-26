#!/bin/bash
#
# Analyze spike-in samples - Lexogen SIRV E2
#

source ../PARAMS.sh

threads=10
assembly="mm10"

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> SIRV MOUSE ANALYSIS <<<"
# Getting only forward/sense mapping reads

if [ -z ${CONDA_PREFIX} ]; then
    echo "Variable \$CONDA_PREFIX is not set. Please make sure you specified if in PARAMS.sh."
    exit
fi

source $CONDA_PREFIX/bin/activate # Source Conda base
conda activate teraseq

source $INSTALL/perl-virtualenv/teraseq/bin/activate

samples=(
    "mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680379.1"
    "mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363659.1"
    "mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680375.1"
    "mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363657.1"
)

workdir=$(pwd)

for i in "${samples[@]}"; do
    echo "Working on $i"

    sdir=$RES_DIR/$i
    mkdir -p $sdir
    mkdir $sdir/primary
    mkdir $sdir/unambig

    samtools view -H $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
        > $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam
    cat $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam > $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam

    samtools view -@ $threads -F 4 -F 16 -F 256 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
        | grep -P "\tSIRV" >> $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam
    samtools view -@ $threads -F 4 -F 16 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
        | grep -w "XN:i:1" | grep -w "XP:i:1" \
        | grep -P "\tSIRV" >> $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam

    samtools view -@ $threads -bh $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam \
        > $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.bam && rm $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam
    samtools index $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.bam

    samtools idxstats $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.bam | grep -i "^SIRV" \
        | sort -k1,1 | cut -f1-3 | sed '1i transcript_id\tlength_idx\tcount' \
        | paste-table-col --ifile - --col-name library --col-val $i \
        > $sdir/primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.idxstats &

    cd $sdir/primary/
    $workdir/src/R/compare_sirv.R $workdir/data/SIRV-annot.tsv
    cd $workdir/

    samtools view -@ $threads -bh $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam \
        > $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.bam && rm $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.sam
    samtools index $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.bam

    bedtools bamtobed -ed -i $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.bam \
    | sort -k 1,1 -k2,2n | sed  '1i chr\tstart\tend\tname\tscore\tstrand' > $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.bed &

    samtools idxstats $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.bam | grep -i "^SIRV" \
        | sort -k1,1 | cut -f1-3 | sed '1i transcript_id\tlength_idx\tcount' \
        | paste-table-col --ifile - --col-name library --col-val $i \
        > $sdir/unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.idxstats &

    cd $sdir/unambig/
    $workdir/src/R/compare_sirv.R $workdir/data/SIRV-annot.tsv
    cd $workdir/
done
wait

mkdir $RES_DIR/SIRV-primary
mkdir $RES_DIR/SIRV-unambig

table-cat $(ls $RES_DIR/*/primary/*-SIRV.sorted.idxstats) > $RES_DIR/SIRV-primary/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.idxstats
table-cat $(ls $RES_DIR/*/unambig/*-SIRV.sorted.idxstats) > $RES_DIR/SIRV-unambig/reads.1.sanitize.noribo.toTranscriptome-SIRV.sorted.idxstats

cd $RES_DIR/SIRV-primary
$workdir/src/R/compare_sirv.R $workdir/data/SIRV-annot.tsv
cd $workdir/

cd $RES_DIR/SIRV-unambig
$workdir/src/R/compare_sirv.R $workdir/data/SIRV-annot.tsv

echo ">>> ALL DONE <<<"
