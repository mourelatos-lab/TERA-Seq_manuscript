#!/bin/bash
#
# Plot heatmap around annotation points
#

source ../PARAMS.sh

threads=16
assembly="hg38"

RES_DIR="results"
mkdir -p $RES_DIR

####################################################################################################

echo ">>> PREPARE REFERENCES <<<"

mkdir $RES_DIR/common

# Split the input bed by type of signals
feats=(
    "PLS"
)

for i in "${feats[@]}"; do
    grep -w $i $DATA_DIR/$assembly/meth/encodeCcreCombined.genome.bed > $RES_DIR/common/encodeCcreCombined.genome.$i.bed
done

echo ">> MAKE BIGWIG FROM BED <<"
# https://www.biostars.org/p/150036/

conda activate teraseq

for i in $RES_DIR/common/*.bed; do
    awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' $i | sort -k1,1 -k2,2n > ${i%.bed}.bedgraph
    bedGraphToBigWig ${i%.bed}.bedgraph $DATA_DIR/$assembly/chrom.sizes ${i%.bed}.bw &
done
wait

echo ">>>> READS AS CENTER POINT - FROM CCRE BED <<<<"

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
)

echo ">> BAM TO BED <<"

for i in "${samples[@]}"; do
    echo "Working for $i"
    sdir=$RES_DIR/$i
    mkdir -p $sdir

    samtools view -bh -F 4 -F 256 -F 2048 $SAMPLE_DIR/$i/align/reads.1.sanitize.toGenome.sorted.bam \
        | bedtools bamtobed -i stdin \
        > $sdir/reads.1.sanitize.toGenome.sorted.bed &
done
wait

echo ">> SPLIT W/WO REL5 <<"

for i in "${samples[@]}"; do
    echo "Working for $i"
    sdir=$RES_DIR/$i
    mkdir -p $sdir

    subset-table -i $sdir/reads.1.sanitize.toGenome.sorted.bed \
        --nonames -c 4 \
        --ilist $SAMPLE_DIR/$i/fastq/reads.1.sanitize.w_rel5.names.txt \
        -o $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.bed &

    subset-table -i $sdir/reads.1.sanitize.toGenome.sorted.bed \
        --nonames -c 4 \
        --ilist $SAMPLE_DIR/$i/fastq/reads.1.sanitize.wo_rel5.names.txt \
        -o $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.bed &
done
wait

echo "> COLLAPSE POSITIONS - 5' <"

for i in "${samples[@]}"; do
    echo "Working for $i"
    sdir=$RES_DIR/$i
    mkdir -p $sdir

    gzip $sdir/reads.1.sanitize.toGenome.sorted.bed &

    # w rel5
    cat $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.bed | grep -P "\t\+$" | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$2+1,$4,$5,$6}' \
        | sort -k 1,1 -k2,2n -T $sdir > $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.plus.bed &
    cat $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.bed | grep -P "\t\-$" | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$3-1,$3,$4,$5,$6}' \
        | sort -k 1,1 -k2,2n -T $sdir > $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.minus.bed &
    wait
    gzip $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.bed &

    cat $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.plus.bed | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$3,".",".",$6}' \
        | uniq > $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.plus.collaps.bed &
    cat $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.minus.bed | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$3,".",".",$6}' \
        | uniq > $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.minus.collaps.bed &
    wait
    gzip $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.plus.bed &
    gzip $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.minus.bed &

    cat $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.plus.collaps.bed $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.minus.collaps.bed \
        | sort -k 1,1 -k2,2n -T $sdir > $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.collaps.bed
    gzip $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.plus.collaps.bed &
    gzip $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.minus.collaps.bed &

    # wo rel5
    cat $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.bed | grep -P "\t\+$" | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$2+1,$4,$5,$6}' \
        | sort -k 1,1 -k2,2n -T $sdir > $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.plus.bed &
    cat $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.bed | grep -P "\t\-$" | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$3-1,$3,$4,$5,$6}' \
        | sort -k 1,1 -k2,2n -T $sdir > $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.minus.bed &
    wait
    gzip $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.bed &

    cat $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.plus.bed | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$3,".",".",$6}' \
        | uniq > $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.plus.collaps.bed &
    cat $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.minus.bed | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$3,".",".",$6}' \
        | uniq > $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.minus.collaps.bed &
    wait
    gzip $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.plus.bed &
    gzip $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.minus.bed &

    cat $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.plus.collaps.bed $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.minus.collaps.bed \
        | sort -k 1,1 -k2,2n -T $sdir > $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.collaps.bed
    gzip $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.plus.collaps.bed &
    gzip $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.minus.collaps.bed &
done
wait

echo ">> MAKE HEATMAPS - REL5 <<"

conda deactivate
source $INSTALL/deepTools-3.5.0/venv/bin/activate

# Set temp dirs for DeepTools, one of them should work; otherwise it might fail because of an empty tmp dir
mkdir $RES_DIR/tmpdir
export TMPDIR=$RES_DIR/tmpdir
export TMP=$RES_DIR/tmpdir
export TEMP=$RES_DIR/tmpdir

for i in "${samples[@]}"; do
    echo "Working for $i"
    sdir=$RES_DIR/$i
    mkdir -p $sdir

    rnd=$RANDOM

    # split reference bed into smaller chunks to speed it up
    ### w rel5
    split -l 5000 $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.collaps.bed $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.collaps.chunks${rnd}

    for chunk in $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.collaps.chunks${rnd}*; do
        # Rename name column (4) in bed to avoid potential problems which deepTools naming which might happen if the reference position name are not unique
        name=$(basename $chunk)
        name=${name##*.}

        cat $chunk | awk -v name=$name 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$3,name,$5,$6}' > $sdir/tmp && mv $sdir/tmp $chunk
    done

    ### wo rel5
    split -l 5000 $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.collaps.bed $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.collaps.chunks${rnd}

    for chunk in $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.collaps.chunks${rnd}*; do
        # Rename name column (4) in bed to avoid potential problems which deepTools naming which might happen if the reference position name are not unique
        name=$(basename $chunk)
        name=${name##*.}

        cat $chunk | awk -v name=$name 'BEGIN {FS = "\t"; OFS = "\t"} {print $1,$2,$3,name,$5,$6}' > $sdir/tmp && mv $sdir/tmp $chunk
    done

    for feat in "${feats[@]}"; do
        echo $feat

        ### w rel5
        for chunk in $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.collaps.chunks${rnd}*; do
            computeMatrix reference-point \
                --referencePoint TSS \
                -R $chunk \
                -S $RES_DIR/common/encodeCcreCombined.genome.${feat}.bw \
                --samplesLabel ${feat} \
                -b 500 -a 500 \
                --skipZeros \
                --missingDataAsZero \
                --binSize 10 \
                --averageTypeBins median \
                --numberOfProcessors $threads \
                --outFileName ${chunk}.gz
        done

        # merge the chunks back to one file
        computeMatrixOperations rbind -m $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.collaps.chunks${rnd}*.gz \
            -o $sdir/matrix.$feat.w_rel5.gz && rm $sdir/reads.1.sanitize.toGenome.sorted.w_rel5.collaps.chunks${rnd}*.gz
        # relabel
        computeMatrixOperations relabel -m $sdir/matrix.$feat.w_rel5.gz -o $sdir/matrix.$feat.w_rel5.tmp.gz \
            --groupLabels "w_rel5" && mv $sdir/matrix.$feat.w_rel5.tmp.gz $sdir/matrix.$feat.w_rel5.gz

        plotHeatmap \
            -m $sdir/matrix.$feat.w_rel5.gz \
            --sortUsing mean \
            --averageTypeSummaryPlot mean \
            --missingDataColor "#440154" \
            --colorMap viridis \
            --zMax 100 \
            --linesAtTickMarks \
            --whatToShow "heatmap and colorbar" \
            --refPointLabel "5-end Nanopore" \
            --heatmapHeight 20 \
            --heatmapWidth 10 \
            --dpi 300 \
            --outFileName $sdir/heatmap.$feat.w_rel5.png &

        ### wo rel5
        for chunk in $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.collaps.chunks${rnd}*; do
            computeMatrix reference-point \
                --referencePoint TSS \
                -R $chunk \
                -S $RES_DIR/common/encodeCcreCombined.genome.${feat}.bw \
                --samplesLabel ${feat} \
                -b 500 -a 500 \
                --skipZeros \
                --missingDataAsZero \
                --binSize 10 \
                --averageTypeBins median \
                --numberOfProcessors $threads \
                --outFileName ${chunk}.gz
        done

        # merge the chunks back to one file
        computeMatrixOperations rbind -m $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.collaps.chunks${rnd}*.gz \
            -o $sdir/matrix.$feat.wo_rel5.gz && rm $sdir/reads.1.sanitize.toGenome.sorted.wo_rel5.collaps.chunks${rnd}*.gz
       # relabel
        computeMatrixOperations relabel -m $sdir/matrix.$feat.wo_rel5.gz -o $sdir/matrix.$feat.wo_rel5.tmp.gz \
            --groupLabels "wo_rel5" && mv $sdir/matrix.$feat.wo_rel5.tmp.gz $sdir/matrix.$feat.wo_rel5.gz

        plotHeatmap \
            -m $sdir/matrix.$feat.wo_rel5.gz \
            --sortUsing mean \
            --averageTypeSummaryPlot mean \
            --missingDataColor "#440154" \
            --colorMap viridis \
            --zMax 100 \
            --linesAtTickMarks \
            --whatToShow "heatmap and colorbar" \
            --refPointLabel "5-end Nanopore" \
            --heatmapHeight 20 \
            --heatmapWidth 10 \
            --dpi 300 \
            --outFileName $sdir/heatmap.$feat.wo_rel5.png &

        # Merge both to one for easier comparison
        computeMatrixOperations rbind -m $sdir/matrix.$feat.w_rel5.gz $sdir/matrix.$feat.wo_rel5.gz  \
            -o $sdir/matrix.$feat.w_wo_rel5.gz

        plotHeatmap \
            -m $sdir/matrix.$feat.w_wo_rel5.gz \
            --sortUsing mean \
            --averageTypeSummaryPlot mean \
            --missingDataColor "#440154" \
            --colorMap viridis \
            --zMax 100 \
            --linesAtTickMarks \
            --refPointLabel "5-end Nanopore" \
            --heatmapHeight 20 \
            --heatmapWidth 10 \
            --dpi 300 \
            --outFileName $sdir/heatmap.$feat.w_wo_rel5.png &

        plotProfile \
            -m $sdir/matrix.$feat.w_wo_rel5.gz \
            --numPlotsPerRow 1 \
            --colors "#5BBA47" "black" \
            --refPointLabel "5-end Nanopore" \
            --outFileName $sdir/profile.$feat.w_wo_rel5.png &
    done

    rm $sdir/reads.1.sanitize.toGenome.sorted.*_rel5.collaps.chunks${rnd}*

    wait
done

echo ">>> ALL DONE <<<"
