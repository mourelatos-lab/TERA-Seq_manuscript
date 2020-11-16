#!/bin/bash
#
# Analyze spike-in samples - Lexogen SIRV E2
#

source ../PARAMS.sh

threads=10
assembly="mm10"

RES_DIR="results"
mkdir $RES_DIR

####################################################################################################

echo ">>> SIRV MOUSE ANALYSIS <<<"
# Getting only forward/sense mapping reads

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
cd $workdir/

echo ">>> ADJUST SIRV ANNOTATION <<<"

mkdir $RES_DIR/common

# Make exon bed from gtf
awk 'BEGIN {FS="\t";OFS="\t"} {if ($3=="exon") {print $1, $4-1,$5, $8, $7}}' $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf \
    > $RES_DIR/common/SIRVome_isoforms_C_170612a.E2.exon.bed.tmp
awk 'BEGIN {FS="\t";OFS="\t"} {if ($3=="exon") {print $0}}' $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.gtf \
    | cut -f9 | sed 's/; /;/g' | awk 'BEGIN {FS=";"}{print $3}' | sed 's/"//g' | sed 's/exon_assignment //g' \
    > $RES_DIR/common/SIRVome_isoforms_C_170612a.E2.exon.bed.tmp.names
paste $RES_DIR/common/SIRVome_isoforms_C_170612a.E2.exon.bed.tmp $RES_DIR/common/SIRVome_isoforms_C_170612a.E2.exon.bed.tmp.names \
    | awk 'BEGIN {FS="\t";OFS="\t"} {print $1,$2,$3,$6,$4,$5}' > $RES_DIR/common/SIRVome_isoforms_C_170612a.E2.exon.bed \
    && rm $RES_DIR/common/SIRVome_isoforms_C_170612a.E2.exon.bed.tmp $RES_DIR/common/SIRVome_isoforms_C_170612a.E2.exon.bed.tmp.names

echo ">>> ANALYZE 5TERA SIRV <<<"

samples=(
   "hsa.dRNASeq.SIRV.polyA.REL5.long.2"
)

for i in "${samples[@]}"; do
	echo "Working on $i"

    sdir=$RES_DIR/$i
	mkdir -p $sdir

    mkdir $sdir/unambig

    # Trimmed reads
    samtools view -H $SAMPLE_DIR/$i/align/reads.1.sanitize.toGenome.sorted.bam > $sdir/unambig/reads.1.sanitize.toGenome.sorted.sam
    samtools view -@ $threads -F 4 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.toGenome.sorted.bam \
        | grep -w "XN:i:1" | grep -w "XP:i:1" \
        >> $sdir/unambig/reads.1.sanitize.toGenome.sorted.sam
    samtools view -bh $sdir/unambig/reads.1.sanitize.toGenome.sorted.sam | samtools sort -@ $threads - \
        > $sdir/unambig/reads.1.sanitize.toGenome.sorted.bam && rm $sdir/unambig/reads.1.sanitize.toGenome.sorted.sam
    samtools index $sdir/unambig/reads.1.sanitize.toGenome.sorted.bam

    # Untrimmed reads for visualization
    samtools view -H $SAMPLE_DIR/$i/align/reads.1.sanitize.untrim.toGenome.sorted.bam > $sdir/unambig/reads.1.sanitize.untrim.toGenome.sorted.sam
    samtools view -@ $threads -F 4 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.untrim.toGenome.sorted.bam \
        | grep -w "XN:i:1" | grep -w "XP:i:1" \
        >> $sdir/unambig/reads.1.sanitize.untrim.toGenome.sorted.sam
    samtools view -bh $sdir/unambig/reads.1.sanitize.untrim.toGenome.sorted.sam | samtools sort -@ $threads - \
        > $sdir/unambig/reads.1.sanitize.untrim.toGenome.sorted.bam && rm $sdir/unambig/reads.1.sanitize.untrim.toGenome.sorted.sam
    samtools index $sdir/unambig/reads.1.sanitize.untrim.toGenome.sorted.bam

    # Genomic mapping
    bedtools bamtobed -split -i $sdir/unambig/reads.1.sanitize.toGenome.sorted.bam \
        | sort -k 1,1 -k2,2n | sed  '1i chr\tstart\tend\tname\tscore\tstrand' > $sdir/unambig/reads.1.sanitize.toGenome.sorted.split.bed &
    bedtools bamtobed -ed -i $sdir/unambig/reads.1.sanitize.toGenome.sorted.bam \
        | sort -k 1,1 -k2,2n | sed  '1i chr\tstart\tend\tname\tscore\tstrand' > $sdir/unambig/reads.1.sanitize.toGenome.sorted.bed &
    wait

    # Add overlaps between exons and split read alignments
    # we could problably use bam as a direct input instead bamtobed
    bedtools intersect -a $sdir/unambig/reads.1.sanitize.toGenome.sorted.split.bed -b $RES_DIR/common/SIRVome_isoforms_C_170612a.E2.exon.bed -wao -s \
        | sed  '1i chr\tstart\tend\tname\tscore\tstrand\tchr_ex\tstart_ex\tend_ex\tname_ex\tscore_ex\tstrand_ex\toverlap' \
            > $sdir/unambig/reads.1.sanitize.toGenome.sorted.split.overlaps.bed

    ./src/R/assign-reads-transcripts.R \
        --ibed $sdir/unambig/reads.1.sanitize.toGenome.sorted.bed \
        --ibed_split $sdir/unambig/reads.1.sanitize.toGenome.sorted.split.overlaps.bed \
        --outdir $sdir/unambig/reassign/genome \
        --molarity data/SIRV-molarity-E2.tsv \
        --tolerance_end 20 \
        --tolerance_splice 10 \
        --overlap 3

    # Subset BAM for each SIRV separately based on the assignment
    mkdir $sdir/unambig/reassign/genome/tmpdir
    mkdir $sdir/unambig/reassign/genome/w_wo_rel5

    tail -n+2 $sdir/unambig/reassign/genome/assign_read-to-trans.tsv | cut -f1 | sort | uniq | while read trans; do
        tail -n+2 $sdir/unambig/reassign/genome/assign_read-to-trans.tsv | grep -w $trans | cut -f2 | sort | uniq \
        > $sdir/unambig/reassign/genome/assign_read-to-trans.$trans.reads.txt

        echo "picard FilterSamReads INPUT=$sdir/unambig/reads.1.sanitize.untrim.toGenome.sorted.bam FILTER=includeReadList READ_LIST_FILE=$sdir/unambig/reassign/genome/assign_read-to-trans.$trans.reads.txt OUTPUT=$sdir/unambig/reassign/genome/assign_read-to-trans.$trans.untrim.bam VALIDATION_STRINGENCY=LENIENT"
    done > $sdir/unambig/reassign/genome/cmds.txt
    wait

    cat $sdir/unambig/reassign/genome/cmds.txt | parallel --tmpdir $sdir/unambig/reassign/genome/tmpdir --eta -j $threads --load 95% --noswap '{}'
    rm $sdir/unambig/reassign/genome/cmds.txt

    for bam in $sdir/unambig/reassign/genome/assign_read-to-trans.*.bam; do
        samtools index -@ $threads $bam
    done

    # Subset bam to get w/wo REL5 reads
    for bam in $sdir/unambig/reassign/genome/*bam; do
        echo "picard FilterSamReads INPUT=$bam FILTER=includeReadList READ_LIST_FILE=$SAMPLE_DIR/$i/fastq/reads.1.sanitize.w_rel5.names.txt OUTPUT=$sdir/unambig/reassign/genome/w_wo_rel5/$(basename $bam .bam).w_rel5.bam VALIDATION_STRINGENCY=LENIENT"
        echo "picard FilterSamReads INPUT=$bam FILTER=includeReadList READ_LIST_FILE=$SAMPLE_DIR/$i/fastq/reads.1.sanitize.wo_rel5.names.txt OUTPUT=$sdir/unambig/reassign/genome/w_wo_rel5/$(basename $bam .bam).wo_rel5.bam VALIDATION_STRINGENCY=LENIENT"
    done > $sdir/unambig/reassign/genome/cmds.txt

    cat $sdir/unambig/reassign/genome/cmds.txt | parallel --tmpdir $sdir/unambig/reassign/genome/tmpdir --eta -j $threads --load 95% --noswap '{}'
    rm $sdir/unambig/reassign/genome/cmds.txt
    rm -r $sdir/unambig/reassign/genome/tmpdir

    # Get coverage for a quick look
    # You can also see the individual BAM files in alignment viewers such as IGV with soft-clipping on to see the adapters and poly(A)
    for bam in $sdir/unambig/reassign/genome/w_wo_rel5/assign_read-to-trans.*rel5.bam; do
        trans=`echo $(basename $bam) | sed 's/assign_read-to-trans.//g' | sed 's/.untrim.*//g'` # get transcript name
        lib=`echo $(basename $bam) | sed 's/.*untrim.//' | sed 's/.bam//'` # Get w/wo REL5

        if [[ $trans != "SIRVunassigned" ]]; then
            grep -w $trans $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.E2.bed \
            | cut -f1-6 > ${bam}.bed
            samtools depth -a -b ${bam}.bed -d 10000 $bam | awk -v var="$lib" 'BEGIN {OFS="\t"}{print $0, var}' \
            | sed  '1i chr\tcoord\tcoverage\tLibrary' > ${bam%.bam}.coverage.tab
            rm ${bam}.bed
        fi

        samtools index -@ $threads $bam
    done
    wait

    for tab in $sdir/unambig/reassign/genome/w_wo_rel5/assign_read-to-trans.SIRV*.w_rel5.coverage.tab; do
        echo $tab
        table-cat $tab ${tab%w_rel5.coverage.tab}wo_rel5.coverage.tab > ${tab%w_rel5.coverage.tab}w_wo_rel5.coverage.tab
        ./src/R/coverage_simple.R ${tab%w_rel5.coverage.tab}w_wo_rel5.coverage.tab ${tab%w_rel5.coverage.tab}w_wo_rel5.coverage.pdf
    done

    gs -dFirstPage=1 -dLastPage=1 -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
        -sOutputFile=$sdir/unambig/reassign/genome/w_wo_rel5/assign_read-to-trans.SIRVall.w_wo_rel5.coverage.pdf \
        $sdir/unambig/reassign/genome/w_wo_rel5/assign_read-to-trans.SIRV*.w_wo_rel5.coverage.pdf &
done
wait

echo ">>> ALL DONE <<<"
