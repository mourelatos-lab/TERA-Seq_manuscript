#!/bin/bash
#
# Plot relative position distribution of TERA-Seq and Akron5 and Ribo-Seq
#

source ../PARAMS.sh

assembly="hg38"

SPAN=41
OFFSET=100

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> RELATIVE POS DISTRO 5p-5p - AKRON5 VS 5TERA W ADAPTER <<<"

conda activate teraseq
source $INSTALL/perl-virtualenv/teraseq/bin/activate

pairs=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1" "hsa.Akron5Seq.HeLa.whole.2"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1" "hsa.Akron5Seq.HeLa.whole.2"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1" "hsa.Akron5Seq.HeLa.whole.2"
)
tLen=${#pairs[@]}

for (( i=0; i<${tLen}; i++ ));
do
	n1=${pairs[$i]}
	((i++))
	n2=${pairs[$i]}

	sdir1=$SAMPLE_DIR"/"$n1"/"
	sdir2=$SAMPLE_DIR"/"$n2"/"
	echo "Working for" $n1 "-" $n2
	mkdir -p $RES_DIR/$n1-$n2/
	bin/relative-pos-distro-on-feats \
		--db1 $sdir1/db/sqlite.db \
		--table1 transcr \
		--pos1 5p \
		--where1 "$MATE1 AND $UNAMBIG AND $CODINGTRANSCRIPT AND rel5 IS NOT NULL" \
		--db2 $sdir2/db/sqlite.db \
		--table2 transcr \
		--pos2 5p \
		--where2 "$MATE1 AND $PRIMARY" \
		--collapse2 \
		--bed $DATA_DIR/$assembly/genic_elements.cds.bed \
		--span $SPAN \
		--offset $OFFSET \
		| paste-table-col --ifile - --col-name sample --col-val $n1-$n2.rel5 \
			> $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.tab

	src/R/rel-pos-distro.R \
		--ifile $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.tab \
		--figfile $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.pdf \
		--name $n1-$n2".5p-5p.cds.rel5"
done
wait

echo ">>> RELATIVE POS DISTRO 5p-5p - RIBOSEQ VS 5TERA W ADAPTER <<<"

pairs=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1" "hsa.RiboSeq.HeLa.async.2"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1" "hsa.RiboSeq.HeLa.async.2"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1" "hsa.RiboSeq.HeLa.async.2"

)
tLen=${#pairs[@]}

len=29

for (( i=0; i<${tLen}; i++ ));
do
	n1=${pairs[$i]}
	((i++))
	n2=${pairs[$i]}

	sdir1=$SAMPLE_DIR"/"$n1"/"
	sdir2=$SAMPLE_DIR"/"$n2"/"
	echo "Working for" $n1 "-" $n2
	mkdir -p $RES_DIR/$n1-$n2/
	bin/relative-pos-distro-on-feats \
		--db1 $sdir1/db/sqlite.db \
		--table1 transcr \
		--pos1 5p \
		--where1 "$MATE1 AND $UNAMBIG AND $CODINGTRANSCRIPT AND rel5 IS NOT NULL" \
		--db2 $sdir2/db/sqlite.db \
		--table2 transcr \
		--pos2 5p \
		--where2 "$MATE1 AND $PRIMARY AND query_length == $len" \
		--collapse2 \
		--bed $DATA_DIR/$assembly/genic_elements.cds.bed \
		--span $SPAN \
		--offset $OFFSET \
	| paste-table-col --ifile - --col-name sample --col-val $n1 \
		> $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.len$len.tab

	src/R/rel-pos-distro.R \
		--ifile $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.len$len.tab \
		--figfile $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.len$len.pdf \
		--name $n1-$n2"-5p-5p-cds"
done
wait

echo ">>> RELATIVE POS DISTRO 5p-5p - RIBOSEQ VS AKRON5 <<<"

pairs=(
	"hsa.Akron5Seq.HeLa.whole.2" "hsa.RiboSeq.HeLa.async.2"
)
tLen=${#pairs[@]}

for (( i=0; i<${tLen}; i++ ));
do
	n1=${pairs[$i]}
	((i++))
	n2=${pairs[$i]}

	sdir1=$SAMPLE_DIR"/"$n1"/"
	sdir2=$SAMPLE_DIR"/"$n2"/"
	echo "Working for" $n1 "-" $n2
	mkdir -p $RES_DIR/$n1-$n2/
	bin/relative-pos-distro-on-feats \
		--db1 $sdir1/db/sqlite.db \
		--table1 transcr \
		--pos1 5p \
		--where1 "$MATE1 AND $PRIMARY" \
		--db2 $sdir2/db/sqlite.db \
		--table2 transcr \
		--pos2 5p \
		--where2 "$MATE1 AND $PRIMARY AND query_length == 29" \
		--collapse2 \
		--bed $DATA_DIR/$assembly/genic_elements.cds.bed \
		--span $SPAN \
		--offset $OFFSET \
	| paste-table-col --ifile - --col-name sample --col-val $n1 \
		> $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.len$len.tab

	src/R/rel-pos-distro.R \
		--ifile $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.len$len.tab \
		--figfile $RES_DIR/$n1-$n2/rel-pos.5p_5p.cds.rel5.len$len.pdf \
		--name $n1-$n2"-5p-5p-cds"
done
wait

echo ">>> ALL DONE <<<"
