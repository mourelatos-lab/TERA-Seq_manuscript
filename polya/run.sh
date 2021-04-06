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
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
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
   	"hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
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
	$RES_DIR/hsa.dRNASeq.HeLa.polyA.REL5OH.long.1/rna_tails.rel5.tsv \
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
	$RES_DIR/hsa.dRNASeq.HeLa.polyA.REL5OH.long.1/rna_tails.rel5.tsv \
    $RES_DIR/hsa.dRNASeq.HeLa.total.REL3.2/rna_tails.rel3.tsv \
	| src/R/polya-distro.R \
		--ifile stdin \
		--software "nanopolish" \
		--maxL 300 \
		--figfile $RES_DIR/polya_length_all.pdf &
wait

echo ">>> CHECK HISTONE GENES POLYA <<<"

echo ">> GET HISTONE GENES <<"

mkdir $RES_DIR/common

cat $DATA_DIR/$assembly/ensembl_genes.gtf | grep "gene_name \"HIST" | cut -f 9 | tr ';' '\n' \
    | grep "transcript_id" | cut -d ' ' -f 3 | sed 's/"//g' | sort | uniq > $RES_DIR/common/histones-transcripts.txt

samples=(
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
)

for sample in ${samples[@]}; do

    echo $sample
    sdir=$RES_DIR/$sample
    mkdir -p $sdir

    # get only histone transcripts
    head -1 $sdir/polya_length.tab > $sdir/polya_length.histones.tab  
    grep -wf $RES_DIR/common/histones-transcripts.txt $sdir/polya_length.tab >> $sdir/polya_length.histones.tab
 
    # get all but histone transcripts   
    grep -v -wf $RES_DIR/common/histones-transcripts.txt $sdir/polya_length.tab > $sdir/polya_length.wo_histones.tab

    cat $sdir/polya_length.histones.tab \
        | src/R/polya-distro.R \
            --ifile stdin \
            --software "nanopolish" \
            --maxL 300 \
            --figfile $sdir/polya_length.histones.all.pdf &

    cat $sdir/polya_length.wo_histones.tab \
        | src/R/polya-distro.R \
            --ifile stdin \
            --software "nanopolish" \
            --maxL 300 \
            --figfile $sdir/polya_length.wo_histones.all.pdf &           
done

echo "This section will fail unless you ran Nanopolish poly(A) estimate and you have the fast5 files. If have both please comment out the next line"
echo ">>> ALL DONE <<<" && exit

echo ">>> CHECK NANOPOLISH POLYA SEGMENTATION <<<"
#https://github.com/jts/nanopolish/issues/882
#https://github.com/jts/nanopolish/blob/master/scripts/polya_training/workflow.md

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
)

conda deactivate    

for i in ${samples[@]}; do
    echo $i
    sdir="$RES_DIR/$i/nanopolish"

    mkdir -p $sdir/PASS
    mkdir -p $sdir/NOREGION

    # Get examples of w/wo polya tail (PASS or NOREGION) - 50 random rows
    head -1 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab > $sdir/PASS/list-polya.tab
    grep PASS $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab \
        | shuf -n 50 >> $sdir/PASS/list-polya.tab

    head -1 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab > $sdir/NOREGION/list-polya.tab
    grep NOREGION $SAMPLE_DIR/$i//align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab \
        | shuf -n 50 >> $sdir/NOREGION/list-polya.tab

    # Convert multi fast5 to single fast5
    source $INSTALL/ont-fast5-api/venv/bin/activate

    mkdir $SAMPLE_DIR/$i/fast5_single
    multi_to_single_fast5 -i $SAMPLE_DIR/$i/fast5 -s $SAMPLE_DIR/$i/fast5_single -t $threads

    # Make new fast5_single indexes (and backup the original ones)
    deactivate

    conda activate teraseq

    mkdir $SAMPLE_DIR/$i/fastq/tmp
    mv $SAMPLE_DIR/$i/fastq/*index* $SAMPLE_DIR/$i/fastq/tmp/

    nanopolish index \
        --directory $SAMPLE_DIR/$i/fast5_single/ \
        $SAMPLE_DIR/$i/fastq/reads.1.fastq.gz

    mkdir $SAMPLE_DIR/$i/fastq/fast5_single
    mv $SAMPLE_DIR/$i/fastq/*index* $SAMPLE_DIR/$i/fastq/fast5_single
    mv $SAMPLE_DIR/$i/fastq/tmp/* $SAMPLE_DIR/$i/fastq/
    rmdir $SAMPLE_DIR/$i/fastq/tmp

    conda deactivate
    source $INSTALL/ont-fast5-api/bin/activate    

    # Plot the chosen reads
    cd $CUR_DIR/$sdir/PASS/
    tail -n+2 list-polya.tab | cut -f1 | while read name; do
        python3 $INSTALL/nanopolish-ab9722b/scripts/polya_training/hmmplot.py --read $name \
            $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab \
            $SAMPLE_DIR/$i/fastq/fast5_single/reads.1.fastq.gz.index.readdb # list of reads --out .
        if [ -x "$(command -v convert)" ]; then
            convert *.png all.pdf
        fi
    done

    cd $CUR_DIR/$sdir/NOREGION/
    tail -n+2 list-polya.tab | cut -f1 | while read name; do
        # modified hnnplot.py from filtring for PASS to NOREGION
        python3 $CUR_DIR/src/Python/hmmplot-noregion.py --read $name \
            $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab \
            $SAMPLE_DIR/$i/fastq/fast5_single/reads.1.fastq.gz.index.readdb # list of reads --out .
        if [ -x "$(command -v convert)" ]; then
            convert *.png all.pdf
        fi
    done

    deactivate

    cd $CUR_DIR/
done


echo ">>> ALL DONE <<<"