#!/bin/bash
#
# Merge all 5TERA3 together to one sample for future analysis
#

source ../PARAMS.sh

threads=6

####################################################################################################

conda activate teraseq

i="hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
sdir="$SAMPLE_DIR/$i"

mkdir -p $sdir/align
mkdir $sdir/fastq
mkdir $sdir/db

echo ">>> MERGE GENOME BAMS <<<"

inbam="reads.1.sanitize.toGenome.sorted.bam"

samtools merge -@ $threads $sdir/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/align/$inbam
samtools index $sdir/align/$inbam &

echo ">>> MERGE TRANSCRIPTOME BAMS <<<"

inbam="reads.1.sanitize.noribo.toTranscriptome.sorted.bam"

samtools merge -@ $threads $sdir/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/align/$inbam \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/align/$inbam
samtools index $sdir/align/$inbam &
wait

echo ">>> MERGE W/O ADAPTER LISTS - REL5 <<<"

inlist="reads.1.sanitize.w_rel5.names.txt"

cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/$inlist \
    > $sdir/fastq/$inlist

inlist="reads.1.sanitize.wo_rel5.names.txt"

cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/$inlist \
    > $sdir/fastq/$inlist

echo ">>> MERGE W/O ADAPTER LISTS - REL3 <<<"

inlist="reads.1.sanitize.w_rel3.names.txt"

cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/$inlist \
    > $sdir/fastq/$inlist

inlist="reads.1.sanitize.wo_rel3.names.txt"

cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/$inlist \
    > $sdir/fastq/$inlist

echo ">>> MERGE DB FILES <<<"
# Since we have just a few reads we will merge them. They are more biological than technical replicates

samples=(
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.4"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.5"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.6"
)

## Ttranscript table
for i in "${samples[@]}"; do
    echo " Working for" $i

    sqlite3 $i/db/sqlite.db ".dump transcr" > $sdir/db/$i.sqlite.transcr.sql & #  make sqldump
done
wait

# Merge all sql dumps
# Use one of the samples to get the header
line_num=`head -100 $sdir/db/${samples[0]}.sqlite.transcr.sql | grep -nP "^CREATE TABLE transcr" | cut -d":" -f1` # Get line number/size of header
head -${line_num} $sdir/db/${samples[0]}.sqlite.transcr.sql > $sdir/db/sqlite.transcr.sql # Get the header from ultimate 6

# Check the INDEX & COMMIT of the db and in case we are missing it append it manually
end_line=`cat $sdir/db/${samples[0]}.sqlite.transcr.sql | grep -nP "^CREATE INDEX transcr" | cut -d":" -f1` # Get line number/size of tail
if [ -z "$end_line" ]; then
    echo -e "CREATE INDEX transcr_loc ON transcr (rname, start);\nCOMMIT;" > $sdir/db/sqlite.transcr.sql.tail.tmp # We can just add it manually
else
    tail -n +${end_line} $sdir/db/${samples[0]}.sqlite.transcr.sql > $sdir/db/sqlite.transcr.sql.tail.tmp # Temporarily get the tail to append to the main db from ultimate 6
fi

for i in "${samples[@]}"; do
    echo $i
    add_num=${RANDOM}0000000 # add $RANDOM ten-milion number to the start; assume we don't have more that 10M reads per library

    grep -P "^INSERT" $sdir/db/$i.sqlite.transcr.sql | sed "s/INSERT INTO transcr VALUES(/INSERT INTO transcr VALUES($add_num/g" \
        >> $sdir/db/sqlite.transcr.sql # make unique id by adding ten-milions otherwise we get an error about not-unique id; keep it number makes it easier

    rm $sdir/db/$i.sqlite.transcr.sql
done

cat $sdir/db/sqlite.transcr.sql.tail.tmp >> $sdir/db/sqlite.transcr.sql && rm $sdir/db/sqlite.transcr.sql.tail.tmp # Append the tail and remove the temp

## genome table
for i in "${samples[@]}"; do
    echo " Working for" $i

    sqlite3 $i/db/sqlite.db ".dump genome" > $sdir/db/$i.sqlite.genome.sql & #  make sqldump
done
wait

# Merge all sql dumps
# use one of the samples to get the header
line_num=`head -100 $sdir/db/${samples[0]}.sqlite.genome.sql | grep -nP "^CREATE TABLE genome" | cut -d":" -f1` # Get line number/size of header
head -${line_num} $sdir/db/${samples[0]}.sqlite.genome.sql > $sdir/db/sqlite.genome.sql # Get the header from ultimate 6

# Check the INDEX & COMMIT of the db and in case we are missing it append it manually
end_line=`cat $sdir/db/${samples[0]}.sqlite.genome.sql | grep -nP "^CREATE INDEX genome" | cut -d":" -f1` # Get line number/size of tail
if [ -z "$end_line" ]; then
    echo -e "CREATE INDEX genome_loc ON genome (rname, start);\nCOMMIT;" > $sdir/db/sqlite.genome.sql.tail.tmp # We can just add it manually
else
    tail -n +${end_line} $sdir/db/${samples[0]}.sqlite.genome.sql > $sdir/db/sqlite.genome.sql.tail.tmp # Temporarily get the tail to append to the main db from ultimate 6
fi

for i in "${samples[@]}"; do
    echo $i
    add_num=${RANDOM}0000000 # add $RANDOM ten-milion number to the start; assume we don't have more that 10M reads per library

    db=$i.sqlite.genome.sql

    grep -P "^INSERT" $sdir/db/$i.sqlite.genome.sql | sed "s/INSERT INTO genome VALUES(/INSERT INTO genome VALUES($add_num/g" \
        >> $sdir/db/sqlite.genome.sql # make unique id by adding ten-milions otherwise we get an error about not-unique id; keep it number makes it easier
    rm $sdir/db/$i.sqlite.genome.sql
done

cat $sdir/db/sqlite.genome.sql.tail.tmp >> $sdir/db/sqlite.genome.sql && rm $sdir/db/sqlite.genome.sql.tail.tmp # Append the tail and remove the temp

# Make the "final" db
cat $sdir/db/sqlite.genome.sql $sdir/db/sqlite.transcr.sql | sqlite3 $sdir/db/sqlite.db && rm $sdir/db/sqlite.genome.sql $sdir/db/sqlite.transcr.sql

# Note: If you have the fast5 files you can remove the "exit" and continue with the analysis
echo "This section will fail unless you ran Nanopolish poly(A) estimate. If you have please comment out line with \"exit"\"
echo ">>> ALL DONE<<<" && exit

echo ">>> MERGE NANOPOLISH POLYA <<<"

source $INSTALL/perl-virtualenv/teraseq/bin/activate

inlist="reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab"

table-cat \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/align/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/align/$inlist \
    hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/align/$inlist \
    > $sdir/align/$inlist

echo ">>> ANNOTATE WITH POLYA <<<"
### Make list of polyA reads (>=0) to append them to db
# From Nanopolish:
#	qc_tag is an additional flag used to indicate the validity of the estimate. Generally speaking, you should only use rows of the output file with
#	this value set to PASS; all other rows with (e.g.) the qc_tag set to SUFFCLIP, ADAPTER, etc. display signs of irregularity that indicate that we believe the estimate to be unreliable.
inlist="reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab"

tail -n+2 $sdir/align/$inlist | egrep -w "PASS|NOREGION" | cut -f1,9 > $sdir/align/${inlist%.tab}.filt.tab

annotate-sqlite-with-file --db_col_add "polya" --db_col_bind "qname" \
    --db_tables "genome" --database "$sdir/db/sqlite.db" --ifile "$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab" \
    --round

annotate-sqlite-with-file --db_col_add "polya" --db_col_bind "qname" \
    --db_tables "transcr" --database "$sdir/db/sqlite.db" --ifile "$sdir/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.tab" \
    --round

echo ">>> ALL DONE <<<"
