#!/bin/bash
#
# Plot coverage of transcripts
#
# IMPORTANT: This analysis can be run only after the metacoord_correction and polya analyses are finished
#

source ../PARAMS.sh

threads=10
assembly="hg38"

gtf="$DATA_DIR/$assembly/Homo_sapiens.GRCh38.91.gtf.gz"
transcripts="data/transcripts_select.txt" # Selected transcripts (you can make your own list or use all)

polya_resdir="$DIR/polya/results" # Directory with poly(A) length results

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> SPLICE GTF <<<"

conda activate teraseq

mkdir $RES_DIR/common

# Selected transcripts
grep -f $transcripts $DATA_DIR/$assembly/ensembl_genes.gtf > $RES_DIR/common/$(basename $transcripts .txt).gtf
src/R/splice-gtf.R --gtf $RES_DIR/common/$(basename $transcripts .txt).gtf --ofile $RES_DIR/common/$(basename $transcripts .txt).spliced.bed

echo ">>>> TRANSCRIPTOMIC COVERAGE <<<<"

echo ">>> SELECTED TRANSCRIPTS COVERAGE <<<"

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
	"hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for i in "${samples[@]}"; do
	echo "Working on $i"

	sdir=$RES_DIR/$i
	mkdir -p $sdir/trans

 	# Manually selected transcripts
	samtools view -H $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam > $sdir/trans/coverage.sam
	samtools view -@ $threads -F 4 -F 16 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
        | grep -f $transcripts \
		| grep -w "XN:i:1" | grep -w "XP:i:1" >> $sdir/trans/coverage.sam
	samtools view -@ $threads -bh $sdir/trans/coverage.sam > $sdir/trans/coverage.bam
	samtools index $sdir/trans/coverage.bam
	bedtools bamtobed -ed -i $sdir/trans/coverage.bam | sort -k 1,1 -k2,2n | sed '1i chr\tstart\tend\tname\tscore\tstrand' \
        > $sdir/trans/coverage.bed &
done
wait

echo ">>> ADD CAGE AND POLYA SIGNALS TO THE PLOTS <<<"

echo ">> CONVERT CAGE AND TSS GENOMIC BED TO TRANSCRIPTOMIC <<"

mkdir $RES_DIR/common

# fantom5 cage
# Note: this is identical to metacoord_correction conversion
./src/R/genomic-to-transcriptomic-coord.R --ifile $DATA_DIR/fantom5/HeLa.rep1.hg38.ctss.bed \
	--ofile $RES_DIR/common/HeLa.rep1.hg38.ctss.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &
./src/R/genomic-to-transcriptomic-coord.R --ifile $DATA_DIR/fantom5/HeLa.rep2.hg38.ctss.bed \
	--ofile $RES_DIR/common/HeLa.rep2.hg38.ctss.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &
./src/R/genomic-to-transcriptomic-coord.R --ifile $DATA_DIR/fantom5/HeLa.rep3.hg38.ctss.bed \
	--ofile $RES_DIR/common/HeLa.rep3.hg38.ctss.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &

# net-cage
./src/R/genomic-to-transcriptomic-coord.R --ifile $DATA_DIR/NET-CAGE/Rep1-HeLaS3-NETCAGE-1M_AGT.hg38.ctss_chr.bed \
	--ofile $RES_DIR/common/Rep1-HeLaS3-NETCAGE-1M.hg38.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &
./src/R/genomic-to-transcriptomic-coord.R --ifile $DATA_DIR/NET-CAGE/Rep2-HeLaS3-NETCAGE-1M_ACG.hg38.ctss_chr.bed \
	--ofile $RES_DIR/common/Rep2-HeLaS3-NETCAGE-1M.hg38.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &

# polya
cat $DATA_DIR/polyasite-2.0/atlas.clusters.hg38.2-0.bed | cut -f1-6 > $RES_DIR/common/polyasite-2.0.bed
cat $RES_DIR/common/polyasite-2.0.bed | \
	./src/R/genomic-to-transcriptomic-coord.R --ifile stdin \
	--ofile $RES_DIR/common/atlas.clusters.hg38.2-0.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &
wait

# fantom5 cage
# Note: this is identical to metacoord_correction conversion
cat \
	$RES_DIR/common/HeLa.rep1.hg38.ctss.genom-to-trans-ByExon.bed \
	$RES_DIR/common/HeLa.rep2.hg38.ctss.genom-to-trans-ByExon.bed \
	$RES_DIR/common/HeLa.rep3.hg38.ctss.genom-to-trans-ByExon.bed \
	| sort --parallel=$threads -T $RES_DIR/common -k 1,1 -k2,2n \
	| uniq \
	> $RES_DIR/common/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.bed &

# net-cage
cat \
	$RES_DIR/common/Rep1-HeLaS3-NETCAGE-1M.hg38.genom-to-trans-ByExon.bed \
	$RES_DIR/common/Rep2-HeLaS3-NETCAGE-1M.hg38.genom-to-trans-ByExon.bed \
	| sort --parallel=$threads -T $RES_DIR/common -k 1,1 -k2,2n \
	| uniq \
	> $RES_DIR/common/RepAll-HeLaS3-NETCAGE-1M.hg38.genom-to-trans-ByExon.bed &
wait

# Merge genomic CAGE/NET-CAGE signals
cat \
    $DATA_DIR/fantom5/HeLa.rep1.hg38.ctss.bed \
    $DATA_DIR/fantom5/HeLa.rep2.hg38.ctss.bed \
    $DATA_DIR/fantom5/HeLa.rep3.hg38.ctss.bed \
    | sort --parallel=$threads -T $RES_DIR/common -k 1,1 -k2,2n \
    | uniq \
    > $RES_DIR/common/HeLa.repAll.hg38.ctss.bed &

cat \
    $DATA_DIR/NET-CAGE/Rep1-HeLaS3-NETCAGE-1M_AGT.hg38.ctss_chr.bed \
    $DATA_DIR/NET-CAGE/Rep2-HeLaS3-NETCAGE-1M_ACG.hg38.ctss_chr.bed \
    | sort --parallel=$threads -T $RES_DIR/common -k 1,1 -k2,2n \
    | uniq \
    > $RES_DIR/common/RepAll-HeLaS3-NETCAGE-1M.hg38.bed
wait

echo ">>> SELECTED TRANSCRIPTS COVERAGE <<<"

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
	mkdir $sdir

	# Manually selected
	subset-table -i $sdir/trans/coverage.bed \
		-c name \
		-l $SAMPLE_DIR/$i/fastq/reads.1.sanitize.w_rel5.names.txt | \
		paste-table-col --ifile - \
		--col-name Library \
		--col-val ${i}.w_rel5 \
		> $sdir/trans/coverage.w_rel5.bed
	subset-table -i $sdir/trans/coverage.bed \
		-c name \
		-l $SAMPLE_DIR/$i/fastq/reads.1.sanitize.wo_rel5.names.txt | \
		paste-table-col --ifile - \
		--col-name Library \
		--col-val ${i}.wo_rel5 \
		> $sdir/trans/coverage.wo_rel5.bed

	table-cat $sdir/trans/coverage.w_rel5.bed $sdir/trans/coverage.wo_rel5.bed \
        > $sdir/trans/coverage.w_wo_rel5.bed
done
wait

echo ">> PLOT TRANSCRIPTOMIC COVERAGE <<"

for i in "${samples[@]}"; do
	echo "Working on $i"

	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# fantom5 cage vs netcage
	src/R/transcript-coverage.R --ifile $sdir/trans/coverage.w_wo_rel5.bed \
        --bin 5 --subset 9999 \
		--addlines $RES_DIR/common/$(basename $transcripts .txt).spliced.bed \
		--suffix fantom5-netcage \
		--addpoints $RES_DIR/common/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.bed \
		--addpoints2 $RES_DIR/common/RepAll-HeLaS3-NETCAGE-1M.hg38.genom-to-trans-ByExon.bed \
		--addregions $RES_DIR/common/atlas.clusters.hg38.2-0.genom-to-trans-ByExon.bed > $sdir/trans/coverage.w_wo_rel5.stats.log &
done
wait

echo "> MERGE OUTPUTS <"

for i in "${samples[@]}"; do
	echo "Working on $i"

	sdir=$RES_DIR/$i

	gs -dFirstPage=1 -dLastPage=2 -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
		-sOutputFile=$sdir/trans/coverage-all.fantom5-netcage.pdf $sdir/trans/coverage-*-fantom5-netcage.pdf &
done
wait

echo "> CIP.DECAP VS 5P W ADAPTER <"

mkdir -p $RES_DIR/compare/trans

table-cat \
	$RES_DIR/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/trans/coverage.w_rel5.bed \
	$RES_DIR/hsa.dRNASeq.HeLa.polyA.REL5.long.1/trans/coverage.w_rel5.bed \
	| src/R/transcript-coverage.R --ifile stdin --bin 5 --subset 9999 \
	--addlines $RES_DIR/common/$(basename $transcripts .txt).spliced.bed \
	--suffix fantom5-netcage \
	--addpoints $RES_DIR/common/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.bed \
	--addpoints2 $RES_DIR/common/RepAll-HeLaS3-NETCAGE-1M.hg38.genom-to-trans-ByExon.bed \
	--addregions $RES_DIR/common/atlas.clusters.hg38.2-0.genom-to-trans-ByExon.bed > $RES_DIR/compare/trans/coverage.w_wo_rel5.coverage-stats.log
mv $RES_DIR/compare/*pdf $RES_DIR/compare/trans/

gs -dFirstPage=1 -dLastPage=2 -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
	-sOutputFile=$RES_DIR/compare/trans/coverage-all-fantom5-netcage.pdf $RES_DIR/compare/trans/coverage-*-fantom5-netcage.pdf

echo ">>>> GENOMIC COVERAGE <<<<"

echo ">> GET ONLY PRIMARY MAPPINGS <<"
# Note: Gviz uses GRanges and they cannot handle Minimap2 secondary mappings with no SEQ

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
  	"hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for i in "${samples[@]}"; do
    echo "Working on $i"

    sdir=$RES_DIR/$i
    mkdir -p $sdir/genome

    samtools view -@ $threads -F4 -F256 -bh $SAMPLE_DIR/$i/align/reads.1.sanitize.toGenome.sorted.bam \
        > $sdir/genome/reads.1.sanitize.toGenome.sorted.primary.bam
    samtools index $sdir/genome/reads.1.sanitize.toGenome.sorted.primary.bam &
done
wait

echo ">> COMPARE LIBRARIES <<"

mkdir -p $RES_DIR/compare/genome

./src/R/gviz-CIPdecapvs5P-dirty-wCage.R \
    --nanopore "$RES_DIR/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/genome/reads.1.sanitize.toGenome.sorted.primary.bam" \
    --nanopore2 "$RES_DIR/hsa.dRNASeq.HeLa.polyA.REL5.long.1/genome/reads.1.sanitize.toGenome.sorted.primary.bam" \
    --assembly "hg38" \
    --gtf $gtf \
    --bedcage "$RES_DIR/common/HeLa.repAll.hg38.ctss.bed" \
    --bedcage2 "$RES_DIR/common/RepAll-HeLaS3-NETCAGE-1M.hg38.bed" \
    --subset "$SAMPLE_DIR/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/fastq/reads.1.sanitize.w_rel5.names.txt" \
    --subset2 "$SAMPLE_DIR/hsa.dRNASeq.HeLa.polyA.REL5.long.1/fastq/reads.1.sanitize.w_rel5.names.txt" \
    --list "data/transcripts_select.txt" \
    --transcriptid TRUE \
    --geneid FALSE \
    --genename FALSE \
    --extend 50 \
    --outdir "$RES_DIR/compare/genome"

gs -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
    -sOutputFile=$RES_DIR/compare/genome/all-alignments.pdf $RES_DIR/compare/genome/*-alignments.pdf

echo ">>> COVERAGE PLOTS WITH ILLUMINA <<<"
### Get alignment+coverage plots of Illumina and Nanopore

for i in "${samples[@]}"; do
	echo "Working on $i"

    sdir=$RES_DIR/$i
    mkdir -p $sdir/genome/wIllumina

    src/R/gviz.R \
        --nanopore $sdir/genome/reads.1.sanitize.toGenome.sorted.primary.bam \
        --illumina $SAMPLE_DIR/hsa.RNASeq.HeLa.xxx.polyA.ENCSR000CPR.1/align/reads.1.Aligned.sortedByCoord.out.bam --strand "reverse" \
        --assembly "hg38" --gtf $gtf \
        --bed "$RES_DIR/common/polyasite-2.0.bed" --extend 50 \
        --list data/transcripts_select.txt --transcriptid \
        --subset $SAMPLE_DIR/$i/fastq/reads.1.sanitize.w_rel5.names.txt \
        --outdir $sdir/genome/wIllumina

    gs -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
        -sOutputFile=$sdir/genome/wIllumina/all-alignments.pdf $sdir/genome/wIllumina/*-alignments.pdf &
done
wait

echo ">> VISUALIZE ALIGNMENTS WITH POLYA TAIL - ADD POLYA <<"

for i in "${samples[@]}"; do
    echo "Working on $i"

    sdir=$RES_DIR/$i
    mkdir -p $sdir

    # Remove softclipping to make appending polyA tail easier
    # Alternative here https://ashokragavendran.wordpress.com/2014/10/20/remove-soft-clipped-reads-from-bam-files/ and here https://www.biostars.org/p/84452/
    java -jar $CONDA_PREFIX/share/jvarkit/remove-softlip.jar --samoutputformat BAM $sdir/genome/reads.1.sanitize.toGenome.sorted.primary.bam \
    > $sdir/genome/reads.1.sanitize.toGenome.sorted.primary.wo-softclip.bam

#   For nanopolish: qc_tag is an additional flag used to indicate the validity of the estimate. Generally speaking, you should only use rows of the output file with this value set to PASS;
#       all other rows with (e.g.) the qc_tag set to SUFFCLIP, ADAPTER, etc. display signs of irregularity that indicate that we believe the estimate to be unreliable.
    tail -n+2 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab \
    | egrep -w "PASS|NOREGION" | cut -f1,9 > $sdir/genome/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.short.tab # get simple polya table

    python3 src/Python/append-polya-bam.py $sdir/genome/reads.1.sanitize.toGenome.sorted.primary.wo-softclip.bam \
        $sdir/genome/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.short.tab \
        $sdir/genome/reads.1.sanitize.toGenome.sorted.primary.wo-softclip.w-polya.bam
done

for i in "${samples[@]}"; do
    echo "Working on $i"

    sdir=$RES_DIR/$i
    mkdir -p $sdir/genome

	# Note: This script CANNOT do read subset, you have to subset the input BAM
    ./src/R/gviz-polyA-dirty-wCage.R \
        --nanopore "$sdir/genome/reads.1.sanitize.toGenome.sorted.primary.wo-softclip.w-polya.bam" \
        --name $i \
        --assembly "hg38" \
        --gtf $gtf \
        --fasta "$DATA_DIR/$assembly/genome/genome.fa" \
        --bedcage "$RES_DIR/common/HeLa.repAll.hg38.ctss.bed" \
        --bedcage2 "$RES_DIR/common/RepAll-HeLaS3-NETCAGE-1M.hg38.bed" \
        --list "data/transcripts_select.txt" \
        --transcriptid TRUE \
        --geneid FALSE \
        --genename FALSE \
        --extend 300 \
        --outdir "$sdir/genome"

    gs -dNOPAUSE -dBATCH -dAutoRotatePages=/None -sDEVICE=pdfwrite \
        -sOutputFile=$sdir/genome/all-alignments.pdf $sdir/genome/*-alignments.pdf &
done
wait

echo ">>> ALL DONE <<<"
