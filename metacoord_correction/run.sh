#!/bin/bash
#
# Re-annotata and plot distribution of alignments along the meta-mRNA
#

source ../PARAMS.sh

assembly="hg38"
threads=10

gtf=$(ls $DATA_DIR/$assembly/Homo_sapiens.GRCh38.91.gtf.gz) # If we don't have the gtf named as the original name file downloaded from Ensembl we would have to specify the database version in the R script

RES_DIR="results"
mkdir -p $RES_DIR

##########################################################################################

# Make one annotation correction for all the samples

echo ">>> MAKE A SINGLE CORRECTION FOR ALL SAMPLES <<<"
# CIP.decap w REL5 for 5' end, all polyA libraries for 3' end correction
# We go from sqlite db which we already have

conda activate teraseq

samples=(
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
)

mkdir $RES_DIR/coords-corrected

# We are not using REL5 reads from any library but CIP.decap so we have to delete the REL5 annotation before merging samples together
for i in "${samples[@]}"; do
	echo " Working for" $i

	cp $SAMPLE_DIR/$i/db/sqlite.db $RES_DIR/coords-corrected/$i.sqlite.db # Make copy of db

	sqlite3 $RES_DIR/coords-corrected/$i.sqlite.db "UPDATE transcr SET rel5=NULL;" # Remove all annotations of rel5
	sqlite3 $RES_DIR/coords-corrected/$i.sqlite.db ".dump transcr" > $RES_DIR/coords-corrected/$i.sqlite.transcr.sql #  make sqldump
	rm $RES_DIR/coords-corrected/$i.sqlite.db
done

# We keep rel5 annotation from CIP.decap
sqlite3 $SAMPLE_DIR/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/db/sqlite.db ".dump transcr" \
	> $RES_DIR/coords-corrected/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1.sqlite.transcr.sql

# Merge all sql dumps
line_num=`head -100 $RES_DIR/coords-corrected/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1.sqlite.transcr.sql \
	| grep -nP "^CREATE TABLE transcr" | cut -d":" -f1` # Get line number/size of header
head -${line_num} $RES_DIR/coords-corrected/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1.sqlite.transcr.sql \
	> $RES_DIR/coords-corrected/sqlite.transcr.coords.sql # Get the header from CIP.decap

# Check the INDEX & COMMIT of the db and in case we are missing it append it manually
end_line=`cat $RES_DIR/coords-corrected/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1.sqlite.transcr.sql \
	| grep -nP "^CREATE INDEX transcr" | cut -d":" -f1` # Get line number/size of tail
if [ -z "$end_line" ]; then
    echo -e "CREATE INDEX transcr_loc ON transcr (rname, start);\nCOMMIT;" > $RES_DIR/coords-corrected/sqlite.transcr.coords.sql.tail.tmp # We can just add it manually
else
    tail -n +${end_line} $RES_DIR/coords-corrected/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1.sqlite.transcr.sql > $RES_DIR/coords-corrected/sqlite.transcr.coords.sql.tail.tmp # Temporarily get the tail to append to the main db from ultimate 6
fi

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
)

for i in "${samples[@]}"; do
	echo $i
	add_num=${RANDOM}0000000 # add $RANDOM ten-milion number to the start; assume we don't have more that 10M reads per library

	db=$RES_DIR/coords-corrected/$i.sqlite.transcr.sql

	grep -P "^INSERT" $db | sed "s/INSERT INTO transcr VALUES(/INSERT INTO transcr VALUES($add_num/g" \
		>> $RES_DIR/coords-corrected/sqlite.transcr.coords.sql # make unique id by adding ten-milions otherwise we get an error about not-unique id; keep it number makes it easier

	rm $db
done

cat $RES_DIR/coords-corrected/sqlite.transcr.coords.sql.tail.tmp >> $RES_DIR/coords-corrected/sqlite.transcr.coords.sql \
	&& rm $RES_DIR/coords-corrected/sqlite.transcr.coords.sql.tail.tmp # Append the tail and remove the temp

# Make the "final" db
cat $RES_DIR/coords-corrected/sqlite.transcr.coords.sql | sqlite3 $RES_DIR/coords-corrected/sqlite.transcr.coords.db # load sqldump back to database
rm $RES_DIR/coords-corrected/sqlite.transcr.coords.sql

# Redefine the coordinates
# 5' end
bin/redefine-mrna-coords \
	--db $RES_DIR/coords-corrected/sqlite.transcr.coords.db \
	--table transcr \
	--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
	--bed-cds $DATA_DIR/$assembly/genic_elements.cds.bed \
	--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel5 IS NOT NULL" \
	--thres 2 \
	> $RES_DIR/coords-corrected/new-mrna-coords.5p.primary.tab &
bin/redefine-mrna-coords \
	--db $RES_DIR/coords-corrected/sqlite.transcr.coords.db \
	--table transcr \
	--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
	--bed-cds $DATA_DIR/$assembly/genic_elements.cds.bed \
	--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel5 IS NOT NULL" \
	--thres 2 \
	> $RES_DIR/coords-corrected/new-mrna-coords.5p.unambig.tab &
# 3' end
bin/redefine-mrna-coords \
	--db $RES_DIR/coords-corrected/sqlite.transcr.coords.db \
	--table transcr \
	--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
	--bed-cds $DATA_DIR/$assembly/genic_elements.cds.bed \
	--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY" \
	--thres 2 \
	> $RES_DIR/coords-corrected/new-mrna-coords.3p.primary.tab &
bin/redefine-mrna-coords \
	--db $RES_DIR/coords-corrected/sqlite.transcr.coords.db \
	--table transcr \
	--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
	--bed-cds $DATA_DIR/$assembly/genic_elements.cds.bed \
	--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG" \
	--thres 2 \
	> $RES_DIR/coords-corrected/new-mrna-coords.3p.unambig.tab &
wait

echo "> REINTEGRATE THE ANNOTATION CORRECTION <"

src/R/integrate-coord-correction.R \
	--ifile5p $RES_DIR/coords-corrected/new-mrna-coords.5p.primary.tab \
	--ifile3p $RES_DIR/coords-corrected/new-mrna-coords.3p.primary.tab \
	--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
	--column_ifile 1 \
	--column_bed 1 \
	--ofile $RES_DIR/coords-corrected/genic_elements.mrna.primary.tab &

src/R/integrate-coord-correction.R \
	--ifile5p $RES_DIR/coords-corrected/new-mrna-coords.5p.unambig.tab \
	--ifile3p $RES_DIR/coords-corrected/new-mrna-coords.3p.unambig.tab \
	--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
	--column_ifile 1 \
	--column_bed 1 \
	--ofile $RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab &
wait

echo ">>> COORDINATES OF READS W/WO REL5 TO META-COORDINATES ON MRNAS (WITH CORRECTED) <<<"

source $INSTALL/perl-virtualenv/teraseq/bin/activate

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
	"hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

# For the corrected valuse we use ONLY correction from the REL5 libraries to be sure we are capturing the end
for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Primary
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $RES_DIR/coords-corrected/genic_elements.mrna.primary.tab \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel5 IS NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-norel5 \
		> $sdir/coords-on-corrected-meta-mrnas.primary.norel5.tab &
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $RES_DIR/coords-corrected/genic_elements.mrna.primary.tab \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel5 IS NOT NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-rel5 \
		> $sdir/coords-on-corrected-meta-mrnas.primary.rel5.tab &

	# Unambiguous
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel5 IS NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-norel5 \
		> $sdir/coords-on-corrected-meta-mrnas.unambig.norel5.tab &
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel5 IS NOT NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-rel5 \
		> $sdir/coords-on-corrected-meta-mrnas.unambig.rel5.tab &
        wait
done

echo ">>> PLOT META-COORDINATES OF READS W/WO REL5 ON MRNAS (WITH CORRECTED) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir/corrected

	# Primary
	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.primary.norel5.tab \
		$sdir/coords-on-corrected-meta-mrnas.primary.rel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.primary.norel5_vs_rel5.pdf &
	cat \
		$sdir/coords-on-corrected-meta-mrnas.primary.norel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.primary.norel5.pdf &
	cat \
		$sdir/coords-on-corrected-meta-mrnas.primary.rel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.primary.rel5.pdf &

	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-corrected-meta-mrnas.primary.norel5.tab \
		--sub 50000 \
		--figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.primary.norel5.pdf &
	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-corrected-meta-mrnas.primary.rel5.tab \
		--sub 50000 \
		--figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.primary.rel5.pdf &

	# Unambiguous
	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.unambig.norel5.tab \
		$sdir/coords-on-corrected-meta-mrnas.unambig.rel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.norel5_vs_rel5.pdf &
	cat \
		$sdir/coords-on-corrected-meta-mrnas.unambig.norel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.norel5.pdf &
	cat \
		$sdir/coords-on-corrected-meta-mrnas.unambig.rel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.rel5.pdf &

	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-corrected-meta-mrnas.unambig.norel5.tab \
		--sub 50000 \
		--figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.unambig.norel5.pdf &
	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-corrected-meta-mrnas.unambig.rel5.tab \
		--sub 50000 \
		--figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.unambig.rel5.pdf &
    wait
done

echo ">>> COORDINATES OF READS W/WO REL5 TO META-COORDINATES ON MRNAS (ON DEFAULT) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Primary
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel5 IS NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-norel5 \
		> $sdir/coords-on-meta-mrnas.primary.norel5.tab &
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel5 IS NOT NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-rel5 \
		> $sdir/coords-on-meta-mrnas.primary.rel5.tab &

	# Unambiguous
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel5 IS NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-norel5 \
		> $sdir/coords-on-meta-mrnas.unambig.norel5.tab &
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel5 IS NOT NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-rel5 \
		> $sdir/coords-on-meta-mrnas.unambig.rel5.tab &
    wait
done

echo "> PLOT META-COORDINATES OF READS W/WO REL5 ON MRNAS (ON DEFAULT) <"
for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir/uncorrected

	# Primary
	table-cat \
		$sdir/coords-on-meta-mrnas.primary.norel5.tab \
		$sdir/coords-on-meta-mrnas.primary.rel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.primary.norel5_vs_rel5.pdf &
	cat \
		$sdir/coords-on-meta-mrnas.primary.norel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.primary.norel5.pdf &
	cat \
		$sdir/coords-on-meta-mrnas.primary.rel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.primary.rel5.pdf &


	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-meta-mrnas.primary.norel5.tab \
		--sub 50000 \
		--figfile $sdir/uncorrected/coords-on-meta-mrnas-heatmap.primary.norel5.pdf &
	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-meta-mrnas.primary.rel5.tab \
		--sub 50000 \
		--figfile $sdir/uncorrected/coords-on-meta-mrnas-heatmap.primary.rel5.pdf &

	# Unambiguous
	table-cat \
		$sdir/coords-on-meta-mrnas.unambig.norel5.tab \
		$sdir/coords-on-meta-mrnas.unambig.rel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.unambig.norel5_vs_rel5.pdf &
	cat \
		$sdir/coords-on-meta-mrnas.unambig.norel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.unambig.norel5.pdf &
	cat \
		$sdir/coords-on-meta-mrnas.unambig.rel5.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.unambig.rel5.pdf &

	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-meta-mrnas.unambig.norel5.tab \
		--sub 50000 \
		--figfile $sdir/uncorrected/coords-on-meta-mrnas-heatmap.unambig.norel5.pdf &
	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-meta-mrnas.unambig.rel5.tab \
		--sub 50000 \
		--figfile $sdir/uncorrected/coords-on-meta-mrnas-heatmap.unambig.rel5.pdf &
    wait
done

echo ">>> COORDINATES OF READS REL3 TO META-COORDINATES ON MRNAS (WITH CORRECTED) <<<"
# For the corrected valuse we use ONLY correction from the REL5 libraries to be sure we are capturing the end

samples=(
	"hsa.dRNASeq.HeLa.total.REL3.1"
	"hsa.dRNASeq.HeLa.total.REL3.2"
	"hsa.dRNASeq.HeLa.total.REL3.3"
	"hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Primary
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $RES_DIR/coords-corrected/genic_elements.mrna.primary.tab \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel3 IS NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-norel3 \
		> $sdir/coords-on-corrected-meta-mrnas.primary.norel3.tab &
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $RES_DIR/coords-corrected/genic_elements.mrna.primary.tab \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel3 IS NOT NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-rel3 \
		> $sdir/coords-on-corrected-meta-mrnas.primary.rel3.tab &

	# Unambiguous
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel3 IS NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-norel3 \
		> $sdir/coords-on-corrected-meta-mrnas.unambig.norel3.tab &
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel3 IS NOT NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-rel3 \
		> $sdir/coords-on-corrected-meta-mrnas.unambig.rel3.tab &
    wait
done

echo ">>> PLOT META-COORDINATES OF READS W/WO REL3 ON MRNAS (WITH CORRECTED) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir/corrected

	# Primary
	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.primary.norel3.tab \
		$sdir/coords-on-corrected-meta-mrnas.primary.rel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.primary.norel3_vs_rel3.pdf &
	cat \
		$sdir/coords-on-corrected-meta-mrnas.primary.norel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.primary.norel3.pdf &
	cat \
		$sdir/coords-on-corrected-meta-mrnas.primary.rel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.primary.rel3.pdf &

	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-corrected-meta-mrnas.primary.norel3.tab \
		--sub 50000 \
		--figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.primary.norel3.pdf &
	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-corrected-meta-mrnas.primary.rel3.tab \
		--sub 50000 \
		--figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.primary.rel3.pdf &

	# Unambiguous
	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.unambig.norel3.tab \
		$sdir/coords-on-corrected-meta-mrnas.unambig.rel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.norel3_vs_rel3.pdf &
	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.unambig.norel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.norel3.pdf &
	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.unambig.rel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.rel3.pdf &

	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-corrected-meta-mrnas.unambig.norel3.tab \
		--sub 50000 \
		--figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.unambig.norel3.pdf &
	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-corrected-meta-mrnas.unambig.rel3.tab \
		--sub 50000 \
		--figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.unambig.rel3.pdf &
    wait
done

echo ">>> COORDINATES OF READS W/WO REL3 TO META-COORDINATES ON MRNAS (ON DEFAULT) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir

	# Primary
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel3 IS NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-norel3 \
		> $sdir/coords-on-meta-mrnas.primary.norel3.tab &
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $PRIMARY AND rel3 IS NOT NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-rel3 \
		> $sdir/coords-on-meta-mrnas.primary.rel3.tab &

	# Unambiguous
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel3 IS NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-norel3 \
		> $sdir/coords-on-meta-mrnas.unambig.norel3.tab &
	bin/coords-on-meta-feats \
		--db $SAMPLE_DIR/$i/db/sqlite.db \
		--table transcr \
		--bed $DATA_DIR/$assembly/genic_elements.mrna.bed \
		--where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND rel3 IS NOT NULL" \
		--bins 20 \
		--min-len 200 \
		| table-paste-col --table - --col-name group1 --col-val $i-rel3 \
		> $sdir/coords-on-meta-mrnas.unambig.rel3.tab &
    wait
done

echo ">>> PLOT META-COORDINATES OF READS W/WO REL3 ON MRNAS (ON DEFAULT) <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir/uncorrected

	# Primary
	table-cat \
		$sdir/coords-on-meta-mrnas.primary.norel3.tab \
		$sdir/coords-on-meta-mrnas.primary.rel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.primary.norel3_vs_rel3.pdf &
	cat \
		$sdir/coords-on-meta-mrnas.primary.norel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.primary.norel3.pdf &
	cat \
		$sdir/coords-on-meta-mrnas.primary.rel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.primary.rel3.pdf &

	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-meta-mrnas.primary.norel3.tab \
		--sub 50000 \
		--figfile $sdir/uncorrected/coords-on-meta-mrnas-heatmap.primary.norel3.pdf &
	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-meta-mrnas.primary.rel3.tab \
		--sub 50000 \
		--figfile $sdir/uncorrected/coords-on-meta-mrnas-heatmap.primary.rel3.pdf &

	# Unambiguous
	table-cat \
		$sdir/coords-on-meta-mrnas.unambig.norel3.tab \
		$sdir/coords-on-meta-mrnas.unambig.rel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.unambig.norel3_vs_rel3.pdf &
	cat \
		$sdir/coords-on-meta-mrnas.unambig.norel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.unambig.norel3.pdf &
	cat \
		$sdir/coords-on-meta-mrnas.unambig.rel3.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
			--subset 0.3 \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.unambig.rel3.pdf &

	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-meta-mrnas.unambig.norel3.tab \
		--sub 50000 \
		--figfile $sdir/uncorrected/coords-on-meta-mrnas-heatmap.unambig.norel3.pdf &
	src/R/plot-meta-coords-heatmap.R \
		--ifile $sdir/coords-on-meta-mrnas.unambig.rel3.tab \
		--sub 50000 \
		--figfile $sdir/uncorrected/coords-on-meta-mrnas-heatmap.unambig.rel3.pdf &
    wait
done

echo ">>> COORDINATES OF READS W/WO POLYA TO META-COORDINATES ON MRNAS <<<"
echo "Note: This section will succeed only if you ran Nanopolish section in the sample processing section."

samples=(
	"hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for i in "${samples[@]}"; do
    echo " Working for" $i;
    sdir=$RES_DIR/$i
    mkdir -p $sdir

    # Unambiguous
    bin/coords-on-meta-feats \
        --db $SAMPLE_DIR/$i/db/sqlite.db \
        --table transcr \
        --bed $RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab \
        --where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND (rel5 IS NOT NULL) AND polya > 0" \
        --bins 20 \
        --min-len 200 \
        | table-paste-col --table - --col-name group1 --col-val $i-rel5 \
        > $sdir/coords-on-corrected-meta-mrnas.unambig.rel5.coding.w_polya.tab &

    bin/coords-on-meta-feats \
        --db $SAMPLE_DIR/$i/db/sqlite.db \
        --table transcr \
        --bed $RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab \
        --where "$MATE1 AND $CODINGTRANSCRIPT AND $UNAMBIG AND (rel5 IS NOT NULL) AND polya == 0" \
        --bins 20 \
        --min-len 200 \
        | table-paste-col --table - --col-name group1 --col-val $i-rel5 \
        > $sdir/coords-on-corrected-meta-mrnas.unambig.rel5.coding.wo_polya.tab &
done
wait

echo ">>> PLOT META-COORDINATES OF READS W/WO POLYA ON MRNAS <<<"

for i in "${samples[@]}"; do
	echo " Working for" $i;
	sdir=$RES_DIR/$i
	mkdir -p $sdir/corrected

	# Unambiguous
    # Coding
	cat \
        $sdir/coords-on-corrected-meta-mrnas.unambig.rel5.coding.w_polya.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
            --subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.rel5.coding.w_polya.pdf &
	cat \
        $sdir/coords-on-corrected-meta-mrnas.unambig.rel5.coding.wo_polya.tab \
		| src/R/plot-meta-coords.R \
			--ifile stdin \
            --subset 0.3 \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.rel5.coding.wo_polya.pdf &
done
wait

echo ">>> GET DETECTED GENES, ADD LENGTHS AND DIVIDE TO QUARTILES <<<"

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

	for a in $sdir/coords-on-corrected-meta-mrnas*tab
	do
		outname=$(basename $a)

		src/R/add-gene-length-quart.R \
			--ifile $sdir/$(basename $a) \
			--ofile $sdir/${outname%.*}-inclLen.tab \
			--gtf $gtf
	done

	for a in $sdir/coords-on-meta-mrnas*tab
	do
		outname=$(basename $a)

		src/R/add-gene-length-quart.R \
			--ifile $sdir/$(basename $a) \
			--ofile $sdir/${outname%.*}-inclLen.tab \
			--gtf $gtf
	done
done

echo ">>> PLOT META MRNA DISTRIB OVER QUARTILES - REL5 <<<"

samples=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
	"hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for i in "${samples[@]}"; do
	echo " Working for" $i
	sdir=$RES_DIR/$i

	# Corrected
	for a in $sdir/coords-on-corrected-meta-mrnas.*.*-inclLen.tab; do
		name=$(basename $a)
		cat $a \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/corrected/${name%.*}-sameQuart.pdf \
			--recalquart \
			--addtransnum &
	done
    wait

	# With recalculated quartiles
	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.primary.norel5-inclLen.tab \
		$sdir/coords-on-corrected-meta-mrnas.primary.rel5-inclLen.tab \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.primary.rel5_vs_norel5-inclLen-sameQuart.pdf \
			--recalquart &

	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.unambig.norel5-inclLen.tab \
		$sdir/coords-on-corrected-meta-mrnas.unambig.rel5-inclLen.tab \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.rel5_vs_norel5-inclLen-sameQuart.pdf \
			--recalquart &
    wait
	# Uncorrected
	for a in $sdir/coords-on-meta-mrnas.*.*-inclLen.tab; do
		name=$(basename $a)
		cat $a \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/uncorrected/${name%.*}-sameQuart.pdf \
			--recalquart \
			--addtransnum &
	done
    wait

	table-cat \
		$sdir/coords-on-meta-mrnas.primary.norel5-inclLen.tab \
		$sdir/coords-on-meta-mrnas.primary.rel5-inclLen.tab \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.primary.rel5_vs_norel5-inclLen-sameQuart.pdf \
			--recalquart &

	table-cat \
		$sdir/coords-on-meta-mrnas.unambig.norel5-inclLen.tab \
		$sdir/coords-on-meta-mrnas.unambig.rel5-inclLen.tab \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.unambig.rel5_vs_norel5-inclLen-sameQuart.pdf \
			--recalquart &
    wait
done

echo ">>> PLOT META MRNA DISTRIB OVER QUARTILES - REL3 <<<"

samples=(
	"hsa.dRNASeq.HeLa.total.REL3.1"
	"hsa.dRNASeq.HeLa.total.REL3.2"
	"hsa.dRNASeq.HeLa.total.REL3.3"
	"hsa.dRNASeq.HeLa.total.REL5.long.REL3.X"
)

for i in "${samples[@]}"; do
	echo " Working for" $i
	sdir=$RES_DIR/$i

	# Corrected
	for a in $sdir/coords-on-corrected-meta-mrnas.*.*-inclLen.tab; do
		name=$(basename $a)
		cat $a \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/corrected/${name%.*}-sameQuart.pdf \
			--recalquart \
			--addtransnum &
	done
    wait

	# With recalculated quartiles
	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.primary.norel3-inclLen.tab \
		$sdir/coords-on-corrected-meta-mrnas.primary.rel3-inclLen.tab \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.primary.rel3_vs_norel3-inclLen-sameQuart.pdf \
			--recalquart &

	table-cat \
		$sdir/coords-on-corrected-meta-mrnas.unambig.norel3-inclLen.tab \
		$sdir/coords-on-corrected-meta-mrnas.unambig.rel3-inclLen.tab \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/corrected/coords-on-corrected-meta-mrnas.unambig.rel3_vs_norel3-inclLen-sameQuart.pdf \
			--recalquart &
    wait

	# Uncorrected
	for a in $sdir/coords-on-meta-mrnas.*.*-inclLen.tab; do
		name=$(basename $a)
		cat $a \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/uncorrected/${name%.*}-sameQuart.pdf \
			--recalquart \
			--addtransnum &
	done
    wait

	table-cat \
		$sdir/coords-on-meta-mrnas.primary.norel3-inclLen.tab \
		$sdir/coords-on-meta-mrnas.primary.rel3-inclLen.tab \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.primary.rel3_vs_norel3-inclLen-sameQuart.pdf \
			--recalquart &

	table-cat \
		$sdir/coords-on-meta-mrnas.unambig.norel3-inclLen.tab \
		$sdir/coords-on-meta-mrnas.unambig.rel3-inclLen.tab \
		| src/R/plot-meta-coords-inclLen.R \
			--ifile stdin \
			--figfile $sdir/uncorrected/coords-on-meta-mrnas.unambig.rel3_vs_norel3-inclLen-sameQuart.pdf \
			--recalquart &
    wait
done

echo ">>>> ANALYZE META-COORDINATES - CAGE <<<<"

sdir="$RES_DIR/CAGE"
mkdir -p $sdir

echo ">> CONVERT COORDS FROM GENOMIC TO TRANSCRIPTOMIC <<"

mkdir $RES_DIR/common

# fantom5 cage
./src/R/genomic-to-transcriptomic-coord.R --ifile $DATA_DIR/$assembly/fantom5/HeLa.rep1.hg38.ctss.bed \
    --ofile $RES_DIR/common/HeLa.rep1.hg38.ctss.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &
./src/R/genomic-to-transcriptomic-coord.R --ifile $DATA_DIR/$assembly/fantom5/HeLa.rep2.hg38.ctss.bed \
    --ofile $RES_DIR/common/HeLa.rep2.hg38.ctss.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &
./src/R/genomic-to-transcriptomic-coord.R --ifile $DATA_DIR/$assembly/fantom5/HeLa.rep3.hg38.ctss.bed \
    --ofile $RES_DIR/common/HeLa.rep3.hg38.ctss.genom-to-trans.bed --annot $DATA_DIR/$assembly/genes.gtf &
wait

cat \
    $RES_DIR/common/HeLa.rep1.hg38.ctss.genom-to-trans-ByExon.bed \
    $RES_DIR/common/HeLa.rep2.hg38.ctss.genom-to-trans-ByExon.bed \
    $RES_DIR/common/HeLa.rep3.hg38.ctss.genom-to-trans-ByExon.bed \
    | sort --parallel=$threads -T $RES_DIR/common -k 1,1 -k2,2n \
    | uniq \
    > $RES_DIR/common/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.bed

echo ">> SUBSET FOR DETECTED TRANSCRIPTS <<"

subset-table \
    -i "$RES_DIR/common/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.bed" \
    --nonames \
    -c 1 \
    -l data/HeLa_transcripts.txt \
    -o $sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.bed

echo ">> SELECT PRIMARY AND UNIQUE (=UNAMBIG) POSITONS <<"

./src/R/primary-unique-select.R \
    --ifile "$sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.bed" \
    --type "primary" \
    --ofile "$sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.primary.bed" &
./src/R/primary-unique-select.R \
    --ifile "$sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.bed" \
    --type "unique" \
    --ofile "$sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.unambig.bed" &
wait

echo ">>> COORDINATES OF CAGE TO META-COORDINATES ON MRNAS (CORRECTED AND DEFAULT) <<<"

echo ">> FORMAT FOR PLOTTING (CORRECTED AND DEFAULT) <<"

source $INSTALL/perl-virtualenv/teraseq/bin/activate

src/R/format-meta-coords-bed.R \
    --bed "$sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.primary.bed" \
    --trans "$RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab" \
    --len 200 \
    --ofile stdout \
    | table-paste-col --table - --col-name group1 --col-val CAGE \
    > "$sdir/coords-on-corrected-meta-mrnas.HeLa_transcripts.primary.tab" &
src/R/format-meta-coords-bed.R \
    --bed "$sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.unambig.bed" \
    --trans "$RES_DIR/coords-corrected/genic_elements.mrna.unambig.tab" \
    --len 200 \
    --ofile stdout \
    | table-paste-col --table - --col-name group1 --col-val CAGE \
    > "$sdir/coords-on-corrected-meta-mrnas.HeLa_transcripts.unambig.tab" &

src/R/format-meta-coords-bed.R \
    --bed "$sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.primary.bed" \
    --trans "$DATA_DIR/$assembly/genic_elements.mrna.bed" \
    --len 200 \
    --ofile stdout \
    | table-paste-col --table - --col-name group1 --col-val CAGE \
    > "$sdir/coords-on-meta-mrnas.HeLa_transcripts.primary.tab" &
src/R/format-meta-coords-bed.R \
    --bed "$sdir/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.HeLa_transcripts.unambig.bed" \
    --trans "$DATA_DIR/$assembly/genic_elements.mrna.bed" \
    --len 200 \
    --ofile stdout \
    | table-paste-col --table - --col-name group1 --col-val CAGE \
    > "$sdir/coords-on-meta-mrnas.HeLa_transcripts.unambig.tab" &
wait

echo ">>> PLOT META-COORDINATES OF CAGE ON MRNAS (CORRECTED AND DEFAULT) <<<"

mkdir $sdir/corrected
mkdir $sdir/uncorrected

cat \
    "$sdir/coords-on-corrected-meta-mrnas.HeLa_transcripts.unambig.tab" \
    | src/R/plot-meta-coords.R \
    --ifile stdin \
    --figfile $sdir/corrected/coords-on-corrected-meta-mrnas.HeLa_transcripts.unambig.pdf &

src/R/plot-meta-coords-heatmap.R \
    --ifile "$sdir/coords-on-corrected-meta-mrnas.HeLa_transcripts.unambig.tab" \
    --sub 50000 \
    --figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.HeLa_transcripts.unambig.pdf &

cat \
    "$sdir/coords-on-meta-mrnas.HeLa_transcripts.unambig.tab" \
    | src/R/plot-meta-coords.R \
    --ifile stdin \
    --figfile $sdir/uncorrected/coords-on-meta-mrnas.HeLa_transcripts.unambig.pdf &

src/R/plot-meta-coords-heatmap.R \
    --ifile "$sdir/coords-on-meta-mrnas.HeLa_transcripts.unambig.tab" \
    --sub 50000 \
    --figfile $sdir/uncorrected/coords-on-corrected-meta-mrnas-heatmap.HeLa_transcripts.unambig.pdf &

cat \
    "$sdir/coords-on-corrected-meta-mrnas.HeLa_transcripts.primary.tab" \
    | src/R/plot-meta-coords.R \
    --ifile stdin \
    --figfile $sdir/corrected/coords-on-corrected-meta-mrnas.HeLa_transcripts.primary.pdf &

src/R/plot-meta-coords-heatmap.R \
    --ifile "$sdir/coords-on-corrected-meta-mrnas.HeLa_transcripts.primary.tab" \
    --sub 50000 \
    --figfile $sdir/corrected/coords-on-corrected-meta-mrnas-heatmap.HeLa_transcripts.primary.pdf &

cat \
    "$sdir/coords-on-meta-mrnas.HeLa_transcripts.primary.tab" \
    | src/R/plot-meta-coords.R \
    --ifile stdin \
    --figfile $sdir/uncorrected/coords-on-meta-mrnas.HeLa_transcripts.primary.pdf &

src/R/plot-meta-coords-heatmap.R \
    --ifile "$sdir/coords-on-meta-mrnas.HeLa_transcripts.primary.tab" \
    --sub 50000 \
    --figfile $sdir/uncorrected/coords-on-corrected-meta-mrnas-heatmap.HeLa_transcripts.primary.pdf &
wait

echo "> GET DETECTED GENES, ADD LENGTHS AND DIVIDE TO QUARTILES <"

for a in $sdir/coords-on-corrected-meta-mrnas*tab; do
    outname=$(basename $a)

    src/R/add-gene-length-quart.R \
        --ifile $sdir/$(basename $a) \
        --ofile $sdir/${outname%.*}-inclLen.tab \
        --gtf $gtf
done

for a in $sdir/coords-on-meta-mrnas*tab; do
    outname=$(basename $a)

    src/R/add-gene-length-quart.R \
        --ifile $sdir/$(basename $a) \
        --ofile $sdir/${outname%.*}-inclLen.tab \
        --gtf $gtf
done

echo "> PLOT META MRNA DISTRIB OVER QUARTILES - CONTROL <"

for a in $sdir/coords-on-corrected-meta-mrnas.*.*-inclLen.tab; do
    name=$(basename $a)

    cat $a \
        | src/R/plot-meta-coords-inclLen.R \
        --ifile stdin \
        --figfile $sdir/corrected/${name%.*}-sameQuart.pdf \
        --recalquart \
        --addtransnum &
done

for a in $sdir/coords-on-meta-mrnas.*.*-inclLen.tab; do
    name=$(basename $a)

    cat $a \
        | src/R/plot-meta-coords-inclLen.R \
        --ifile stdin \
        --figfile $sdir/uncorrected/${name%.*}-sameQuart.pdf \
        --recalquart \
        --addtransnum &
done

echo ">>> ALL DONE <<<"
