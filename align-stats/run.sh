#!/bin/bash
#
# Get alignment/mapping statistics
#

source ../PARAMS.sh

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> MAPPED LENGTH <<<"

conda activate teraseq

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
)

mkdir $RES_DIR/lengths

for i in "${samples[@]}"; do
    echo $i

    # Longest aligned read - it doesn't mean there is no softclipping or extensive splicing!!!
    samtools view -@ 1 -F 4 -F 256 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.toGenome.sorted.bam \
        | awk '{print $1,$3,length($10)}' | sort -T . -k3,3nr -t ' ' > $RES_DIR/lengths/${i}.alignedReadLen.genome.full.txt &
    samtools view -@ 1 -F 4 -F 256 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
        | awk '{print $1,$3,length($10)}' | sort -T . -k3,3nr -t ' ' > $RES_DIR/lengths/${i}.alignedReadLen.transcriptome.full.txt &

    # Longest aligned section of the read (sum of bases); counts M, X, = and D, doesn't count N, I, ...
    samtools view -@ 1 -F 4 -F 256 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.toGenome.sorted.bam \
        | perl -slane '$l = 0; $F[5] =~ s/(\d+)[MX=D]/$l+=$1/eg; print $F[0],"\t",$F[2],"\t",$l' | sort -T . -k3,3nr \
        > $RES_DIR/lengths/${i}.alignedReadPart.genome.full.txt &
    samtools view -@ 1 -F 4 -F 256 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
        | perl -slane '$l = 0; $F[5] =~ s/(\d+)[MX=D]/$l+=$1/eg; print $F[0],"\t",$F[2],"\t",$l' | sort -T . -k3,3nr \
        > $RES_DIR/lengths/${i}.alignedReadPart.transcriptome.full.txt &

    # Longest read - it doesn't mean there is no softclipping or extensive splicing!!!
    samtools view -@ 1 $SAMPLE_DIR/$i/align/reads.1.sanitize.toGenome.sorted.bam | awk '{print length($10)}' \
        | sort | uniq -c | sort -k2,2n > $RES_DIR/lengths/${i}.allReadLen.genome.txt &
    samtools view -@ 1 $SAMPLE_DIR/$i/align/reads.1.sanitize.noribo.toTranscriptome.sorted.bam \
        | awk '{print length($10)}' | sort | uniq -c | sort -k2,2n > $RES_DIR/lengths/${i}.allReadLen.transcriptome.txt &
    wait
done

# Longest aligned read including softclipped parts
for i in $RES_DIR/lengths/*.alignedReadLen.*.full.txt; do
    cat $i | cut -d ' ' -f3 | sort | uniq -c | sort -k2,2n > ${i%.full.txt}.txt

    lib_name=$(basename $i .txt)
    ref_name=`echo $lib_name | sed 's/.*\.alignedReadLen.//' | sed 's/.full//'`
    lib_name=${lib_name%.aligned*}

    echo -ne "$lib_name\t$ref_name\t" >> $RES_DIR/longest-alignedReadLen.tsv
    cat ${i%.full.txt}.txt | sed 's/ \+/\t/g' | sed 's/^\t//g' | cut -f2 | sort -nr | head -1 >> $RES_DIR/longest-alignedReadLen.tsv
done

# Longest aligned read part excluding softclipped parts - more realistic read alignment stats
for i in $RES_DIR/lengths/*.alignedReadPart.*.full.txt; do
    cat $i | cut -f3 | sort | uniq -c | sort -k2,2n > ${i%.full.txt}.txt

    lib_name=$(basename $i .txt)
    ref_name=`echo $lib_name | sed 's/.*\.alignedReadPart.//' | sed 's/.full//'`
    lib_name=${lib_name%.aligned*}

    echo -ne "$lib_name\t$ref_name\t" >> $RES_DIR/longest-alignedReadPart.tsv
    cat ${i%.full.txt}.txt | sed 's/ \+/\t/g' | sed 's/^\t//g' | cut -f2 | sort -nr | head -1 >> $RES_DIR/longest-alignedReadPart.tsv
done

echo ">>> MAPPING STATS - 5TERA & 5TERA3 (REL5) <<<"

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.4"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.5"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.6"
)

for i in "${samples[@]}"; do
	sdir=$SAMPLE_DIR/$i
	echo "Working for" $i

	echo "All reads mapped to genome (primary)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0));'
	echo "Unambiguously mapped reads to genome (unambiguous, primary flag not considered)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (best_hit == 1 AND number_of_best_hits == 1));'

	echo "All reads mapped to protein-coding transcripts (primary)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
	echo "Unambiguously mapped reads to protein-coding transcripts"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL) AND (best_hit == 1 AND number_of_best_hits == 1));'

	echo "Mapped reads to protein-coding transcripts without REL5 adapter (primary)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel5 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel5 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'

	echo "All mapped protein-coding transcripts without REL5 adapter (from primary mappings)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel5 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
	echo "All mapped protein-coding transcripts with REL5 adapter (from primary mappings)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel5 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
done > $RES_DIR/REL5.stats.txt

echo ">>> MAPPING STATS - 3TERA & 5TERA3 (REL3) <<<"

samples=(
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.4"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.5"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.6"
)

for i in "${samples[@]}"; do
    sdir=$SAMPLE_DIR/$i
    echo "Working for" $i

	echo "All reads mapped to genome (primary)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0));'
	echo "Unambiguously mapped reads to genome (unambiguous, primary flag not considered)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM genome WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (best_hit == 1 AND number_of_best_hits == 1));'

	echo "All reads mapped to protein-coding transcripts (primary)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
	echo "Unambiguously mapped reads to protein-coding transcripts"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (((flag & 1) == 0 OR (flag & 64) == 64) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL) AND (best_hit == 1 AND number_of_best_hits == 1));'

	echo "Mapped reads to protein-coding transcripts without REL3 adapter (primary)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel3 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
	echo "Mapped reads to protein-coding transcripts with REL3 adapter (primary)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT qname) FROM transcr WHERE (rel3 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'

	echo "All mapped protein-coding transcripts without REL3 adapter (from primary mappings)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel3 IS NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
	echo "All mapped protein-coding transcripts with REL3 adapter (from primary mappings)"
	sqlite3 $sdir/db/sqlite.db 'SELECT COUNT (DISTINCT rname) FROM transcr WHERE (rel3 IS NOT NULL AND ((flag & 1) == 0 OR (flag & 64) == 64) AND ((flag & 256) == 0) AND (coding_transcript IS NOT NULL AND noncoding_transcript IS NULL));'
done > $RES_DIR/REL3.stats.txt

echo ">>> ALL DONE <<<"
