#!/bin/bash
#
# Run 5TERA preprocessing, alignment, and postprocessing for SIRV E2 spike-ins
#

source ../PARAMS.sh

threads=6
assembly="hg38"

####################################################################################################

samples=(
	"hsa.dRNASeq.SIRV.polyA.REL5.long.2"
)

echo ">>> MAKE DIRECTORY STRUCTURE <<<"

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	mkdir -p $sdir/logfiles
	mkdir $sdir/align
	mkdir $sdir/db
done

echo ">>> SANITIZE FASTQ HEADERS <<<"

source $INSTALL/perl-virtualenv/teraseq/bin/activate

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	zcat $sdir/fastq/reads.1.fastq.gz \
	| fastq-sanitize-header --input - --delim : --keep 0 \
	| gzip \
	> $sdir/fastq/reads.1.sanitize.fastq.gz &
done
wait

deactivate

echo ">>> REMOVE REL5 ADAPTOR <<<"

conda activate teraseq
source $INSTALL/cutadapt-2.5/venv/bin/activate

# Adapter trimming settings are adjusted for SIVR - they might be different from regular 5TERA/TERA3 trimming!!!

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	cutadapt \
		-g XAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
		--overlap 31 \
		--minimum-length 25 \
		--error-rate 0.30 \
		--output $sdir/fastq/reads.1.sanitize.w_rel5.fastq.gz \
		--untrimmed-output $sdir/fastq/reads.1.sanitize.wo_rel5.fastq.gz \
		$sdir/fastq/reads.1.sanitize.fastq.gz \
		&> $sdir/logfiles/cutadapt.rel5.log &
done
wait

deactivate

echo ">>> MERGE READS WITH AND WITHOUT REL5 ADAPTOR <<<"

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	zcat $sdir/fastq/reads.1.sanitize.w_rel5.fastq.gz | paste - - - - | cut -f1 | sed 's/^@//g' \
		> $sdir/fastq/reads.1.sanitize.w_rel5.names.txt &
	zcat $sdir/fastq/reads.1.sanitize.wo_rel5.fastq.gz | paste - - - - | cut -f1 | sed 's/^@//g' \
		> $sdir/fastq/reads.1.sanitize.wo_rel5.names.txt &

	cat $sdir/fastq/reads.1.sanitize.w_rel5.fastq.gz $sdir/fastq/reads.1.sanitize.wo_rel5.fastq.gz \
		> $sdir/fastq/reads.1.sanitize.rel5_trim.fastq.gz &
done
wait

echo ">>> MAKE A LINK OF READS FOR TRANSCRIPTOME MAPPING <<<"

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	ln -s $sdir/fastq/reads.1.sanitize.rel5_trim.fastq.gz $sdir/fastq/reads.1.sanitize.noribo.fastq.gz
done

echo ">>> ALIGN READS TO TRANSCRIPTOME (WITH SECONDARY) <<<"

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	minimap2 \
		-a \
		-x map-ont \
		-k 12 \
		-p 1 \
		-u n \
		--for-only \
		-t $threads \
		--secondary=yes \
		$DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/minimap2.17/transcripts_sirv1.E2.k12.mmi \
		$sdir/fastq/reads.1.sanitize.noribo.fastq.gz \
  	| samtools view -h -F 16 - \
	| add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
	| grep -v "SA:Z:" \
	| sam-count-secondary --tag X0 \
	| samtools view -b - \
	| samtools sort - \
	> $sdir/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam

	ln -s $sdir/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam \
	$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam
done
wait

echo ">>> ALIGN READS TO GENOME (WITH SECONDARY) <<<"
# added --splice-flank=no for SIRV
# longest intron in SIRV E2 is 38218 - set -G to 40k; will disable chimeric alignemnts, if present

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	minimap2 \
		-a \
		-x splice \
		-k 12 \
		-p 1 \
		-u f \
		-t $threads \
		--secondary=yes \
		--splice-flank=no \
		--junc-bed $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf.E2.tmp.bed12 \
		--junc-bonus 10 \
		-G 40000 \
		$DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/minimap2.17/genome.k12.mmi \
		$sdir/fastq/reads.1.sanitize.rel5_trim.fastq.gz \
	| add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
	| grep -v "SA:Z:" \
	| sam-count-secondary --tag X0 \
	| samtools view -b - \
	| samtools sort - \
	> $sdir/align/reads.1.sanitize.toGenome.sorted.bam
done
wait

echo ">>> ALIGN READS TO GENOME (WITH SECONDARY) WITH UNTRIMMED READS (FOR ADAPTOR VISUALIZATION) <<<"
# added --splice-flank=no for SIRV
# longest intron in SIRV E2 is 38218 - set -G to 40k; will disable chimeric alignemnts, if present

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	minimap2 \
		-a \
		-x splice \
		-k 12 \
		-p 1 \
		-u f \
		-t $threads \
		--secondary=yes \
		--splice-flank=no \
		--junc-bed $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_C_170612a.gtf.E2.tmp.bed12 \
		--junc-bonus 10 \
		-G 40000 \
		$DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/minimap2.17/genome.k12.mmi \
		$sdir/fastq/reads.1.sanitize.fastq.gz \
	| add-tag-max-sam --tag ms --newtag XP --newtag2 XN \
	| grep -v "SA:Z:" \
	| sam-count-secondary --tag X0 \
	| samtools view -b - \
	| samtools sort - \
	> $sdir/align/reads.1.sanitize.untrim.toGenome.sorted.bam
done
wait

echo ">>> INDEX ALIGNMENTS <<<"

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	samtools index \
		$sdir/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam &
	ln -s $sdir/align/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.bam.bai \
		$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam.bai
	samtools index \
		$sdir/align/reads.1.sanitize.toGenome.sorted.bam &
	samtools index \
		$sdir/align/reads.1.sanitize.untrim.toGenome.sorted.bam &
	wait
done

CONDA_PATH=$CONDA_PREFIX # Temporary store path to the Conda environment
conda deactivate

echo ">>> SAM TO SQLITE (TRANSCRIPTOME) <<<"

source $INSTALL/perl-virtualenv/teraseq/bin/activate

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	cat $sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
	| $CONDA_PATH/bin/samtools view -h -F 4 -F 16 -F 2048 - \
	| sam_to_sqlite \
		--database $sdir/db/sqlite.db \
		--table transcr \
		--records_class GenOOx::Data::File::SAMminimap2::Record \
		--drop &
done
wait

echo ">>> SAM TO SQLITE (GENOME) <<<"

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	cat $sdir/align/reads.1.sanitize.toGenome.sorted.bam \
	| $CONDA_PATH/bin/samtools view -h -F 4 -F 2048 - \
	| sam_to_sqlite \
		--database $sdir/db/sqlite.db \
		--table genome \
		--records_class GenOOx::Data::File::SAMminimap2::Record \
		--drop &
done
wait

deactivate

echo ">>> ANNOTATE WITH REL5 <<<"

conda activate teraseq

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	echo ">> ANNOTATE WITH REL5 (TRANSCRIPTOME) <<"

	annotate-sqlite-with-fastq \
		--database $sdir/db/sqlite.db \
		--db_col_bind "qname" \
		--db_col_add "rel5" \
		--db_tables "transcr" \
		--ifile $sdir/fastq/reads.1.sanitize.w_rel5.fastq.gz

	echo ">> ANNOTATE W/WO REL5 (GENOME) <<"

	annotate-sqlite-with-fastq \
		--database $sdir/db/sqlite.db \
		--db_col_bind "qname" \
		--db_col_add "rel5" \
		--db_tables "genome" \
		--ifile $sdir/fastq/reads.1.sanitize.w_rel5.fastq.gz
done

# Note: If you have the fast5 files you can remove the "exit" and continue with the analysis
echo "This section will fail unless you ran Nanopolish poly(A) estimate. If you have fast5 files please comment out line with \"exit"\"
echo ">>> ALL DONE<<<"
exit

echo ">>> NANOPOLISH INDEX <<<"

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	nanopolish index \
		--directory $sdir/fast5/ \
		$sdir/fastq/reads.1.fastq.gz &
done
wait

echo ">>> NANOPOLISH POLY(A) <<<"

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo " Working for" $i

	nanopolish polya \
		--reads $sdir/fastq/reads.1.fastq.gz \
		--bam $sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
		--genome $DATA_DIR/spikein/sirv/SIRV_Set1_Sequences_170612a/SIRVome_isoforms_170612a.transcripts_sirv1.E2.fa \
		--threads $threads \
		> $sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab
done
wait

echo ">>> ALL DONE <<<"
