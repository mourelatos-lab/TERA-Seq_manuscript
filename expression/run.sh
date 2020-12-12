#
# Compare gene/transcript expressions between the libraries
#

source ../PARAMS.sh

threads=10
assembly="hg38"

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> GET TPM EXPRESSION <<<"
# We don't care about having or not having adapters as we want to compare only the expression
# We take only primary alignments
# ONT-bam libraries

conda activate teraseq
source $INSTALL/perl-virtualenv/teraseq/bin/activate

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.1"
    "hsa.dRNASeq.HeLa.polyA.PNK.REL5.1"
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.4"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.5"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.6"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i

	mkdir -p $sdir

	sqlite3 $SAMPLE_DIR/$i/db/sqlite.db 'SELECT rname FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL) AND number_of_mappings == 1);' | \
	sort | uniq -c > $sdir/quant.short.sf.raw # Extract uniquely mapped protein-codin reads

	lib_size=`cat $sdir/quant.short.sf.raw | awk -F' ' '{sum+=$1;} END{print sum/1000000;}'` # Get library size

	echo $lib_size

    lib_name=`echo $i | sed 's/hsa.dRNASeq.HeLa.//g'`
	awk -v var=$lib_size 'BEGIN{FS=" ";OFS="\t";}{print $2,$1/var;}' $sdir/quant.short.sf.raw | \
		sed '1i transcript_id\tTPM' | \
		table-paste-col \
		--table - \
		--col-name lib_name \
		--col-val $lib_name > $sdir/quant.short.sf && rm $sdir/quant.short.sf.raw # Make TPM
done
wait

#echo ">>> ADD TRANSCRIPT LENGTH <<<"
#get-cdna-length transcript_id $DATA_DIR/$assembly/ensembl_transcripts_protein_coding.gtf | sed '1i transcript_id\tcdna_len' > data/ensembl_transcripts.protein-coding.cdna_len.txt
#get-cdna-length transcript_id $DATA_DIR/$assembly/ensembl_genes.gtf | sed '1i transcript_id\tcdna_len' > data/ensembl_transcripts.cdna_len.txt

# Merge all tables and add transcript lengths
for i in $RES_DIR/*/quant*short.sf; do
    echo $i

    table-join --left \
        --table1 $i --key1 "transcript_id" \
        --table2 data/ensembl_transcripts.protein-coding.cdna_len.txt --key1 "transcript_id" \
        | cut -f2-5 > ${i%.sf}.w_len.sf
done

table-cat \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/quant.short.w_len.sf \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1/quant.short.w_len.sf \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.PNK.REL5.1/quant.short.w_len.sf \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.REL5.1/quant.short.w_len.sf \
    $RES_DIR/hsa.dRNASeq.HeLa.polyA.REL5.long.1/quant.short.w_len.sf \
    $RES_DIR/hsa.dRNASeq.HeLa.total.REL3.1/quant.short.w_len.sf \
    $RES_DIR/hsa.dRNASeq.HeLa.total.REL3.2/quant.short.w_len.sf \
    $RES_DIR/hsa.dRNASeq.HeLa.total.REL3.3/quant.short.w_len.sf \
    > $RES_DIR/quant.w_len.tsv

echo ">>> PLOT CORRELATION <<<"

src/R/correlate-exp.R --ifile $RES_DIR/quant.w_len.tsv \
	--figfile $RES_DIR/quant.allExp.pdf # colnames to input to R script: transcript_id\tlib_name\tTPM\tcdna_len

echo ">>> ALL DONE <<<"
