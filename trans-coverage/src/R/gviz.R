#!/usr/bin/env Rscript
#
# Experimental plotting of alignment using Gviz (https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html)
# Book: Statistical Genomics; Ewy Math? ? Sean Davis Editors; Methods and Protocols
# Tutorials: https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/viz_course/Presentations/singlepage/Viz_part_1.html
#   https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/viz_course/Presentations/singlepage/Viz_part_2.html; https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/
#
# Note: This script is organism specific due to biomart inputs
#
# TODO: Change getting the data from biomart to input gtf and genome fasta (I already have it somewhere...)
#

suppressMessages(library("optparse"))
#library(biomaRt)
suppressMessages(library(rtracklayer))
suppressMessages(library(ensembldb))
suppressMessages(library(GenomicFeatures))
#library(IRanges)
suppressMessages(library(Rsamtools))
suppressMessages(library(Gviz))

 option_list <- list(
   make_option(
     c("-n", "--nanopore"), type = "character", metavar = "BAM file",
     help = "Input Nanopore bam file. If you used minimap2 you have to prefilter for primary mappings only."),
   make_option(
     c("-i", "--illumina"), type = "character", metavar = "BAM file",
     help = "Input Illumina bam file. Right now works only for paired-end reads. If single-end reads, change the sam filtering flags."),
   make_option(
     c("-a", "--assembly"), type = "character", default="hg38", metavar = "String",
     help = "Genome reference assembly [hg38 | mm10]"),
   make_option(
     c("-g", "--gtf"), type = "character", metavar = "GTF file",
     help = "Ensembl GTF. Has to have the GTF name."),
   make_option(
     c("-b", "--bed"), type = "character", metavar = "BED file",
     help = "Additional BED annotation. For examples alternative polya database BED."),
   make_option(
     c("-u", "--subset"), type = "character", default="", metavar = "File",
     help = "Read subset list. One read per line with no header. Default: \"\"."),
   make_option(
     c("-s", "--strand"), type = "character", default="reverse", metavar = "String",
     help = "Stradnedness of the Illumina [forward | reverse]. Defaul: reverse."),
   make_option(
     c("-l", "--list"), type = "character", metavar = "File",
     help = "List of transcript ids, gene ids or gene names. One name per line with no header."),
   make_option(
     c("-t", "--transcriptid"), type = "logical", default=FALSE, action="store_true", metavar = "String",
     help = "The --list is Ensembl transcript id list."),
   make_option(
     c("-d", "--geneid"), type = "logical", default=FALSE, action="store_true", metavar = "String",
     help = "The --list is Ensembl gene id list."),
   make_option(
     c("-e", "--genename"), type = "logical", default=FALSE, action="store_true", metavar = "String",
     help = "The --list is Ensembl gene name/symbol."),
   make_option(
     c("-x", "--extend"), type = "integer", default=50, metavar = "Integer",
     help = "Extend the gene region --extend bases to the left and right. Default: 50."),
   make_option(
     c("-o", "--outdir"), type = "character", metavar="String",
     help = "Output directory to save the plots.")
 )
 opt = parse_args(OptionParser(option_list = option_list))

opt$max_depth<-10000 # Max depth for Illumina to avoid RAM draining

### TESTING VARIABLES!!!
#print("TESTING VARIABLES!!!")
#opt<-NULL
## Minimap2 doesn't output SEQ (and QUAL) for secondary alignments and GRanges doesn't like it
#opt$nanopore<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/alignment/reads.1.sanitize.toGenome.sorted.primary.bam" # samtools view -@ 26 -F4 -F256 -bh ~/projects/ribothrypsis/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/align/reads.1.sanitize.toGenome.sorted.bam > CIP.decap.reads.1.sanitize.toGenome.sorted.bam; samtools index CIP.decap.reads.1.sanitize.toGenome.sorted.bam
#opt$illumina<-"/home/joppelt/projects/ribothrypsis/samples/hsa.RNASeq.HeLa.xxx.polyA.ENCSR000CPR.1/align/reads.1.Aligned.sortedByCoord.out.bam"
#opt$assembly<-"hg38"
#opt$gtf<-"/home/joppelt/projects/ribothrypsis/data/hg38/Homo_sapiens.GRCh38.91.gtf.gz"
#opt$bed<-"/home/joppelt/projects/ribothrypsis/analysis/dist-from-annot/results/common/polyasite-2.0.bed"
#opt$subset<-"/home/joppelt/projects/ribothrypsis/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/fastq/reads.1.sanitize.w_rel5.names.txt"
#opt$strand<-"reverse"
#opt$list<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/transcripts_select.CIP.decap.100.txt"
#opt$transcriptid<-TRUE
#opt$geneid<-FALSE
#opt$genename<-FALSE
#opt$extend<-50
#print("TESTING VARIABLES!!!")
### TESTING VARIABLES!!!

inbam<-c(opt$nanopore, opt$illumina) # nanopore, illumina (test)
assembly<-opt$assembly
gtf<-opt$gtf
inbed<-opt$bed
subset_read<-opt$subset # Subset of Nanopore reads (like adapter reads)
strand.ill<-opt$strand # Strandness of the Illumina seq; reverse: R1 on reverse, R2 on forward; forward: R1 on forward, R2 on reverse
inlist<-opt$list # One transcript id/gene id/gene name per line
transcriptid<-opt$transcriptid
geneid<-opt$geneid
genename<-opt$genename
extend<-opt$extend
outdir<-opt$outdir

# Colors
orange<-rgb(250, 164, 26, maxColorValue = 255)
short<-rgb(238, 72, 154, maxColorValue = 255)
long<-rgb(53, 180, 74, maxColorValue = 255)
apa<-rgb(0, 0, 255, maxColorValue = 255)

##########################################
if(nchar(subset_read)>0){
    subset_read<-read.table(subset_read, header=F)
}
inlist<-read.table(inlist, header=F)

# Make TxDb from GTF
txdb<-makeTxDbFromGFF(gtf, format="auto", organism="Homo sapiens", dataSource="Ensembl")

# Load extra BED annotation
bedannot.main <- rtracklayer::import(inbed)

# IMPORTANT: We need to set this to false otherwise Gviz won't show any alignments and will mess up the gene visualization
# This sest ALL chr reference for Gviz to be "1, 2, 3" and not "chr1, chr2, chr3"
# https://support.bioconductor.org/p/78337/
options(ucscChromosomeNames = FALSE)

### Get gene annotation from Ensemlb (or any other GTF)
grtrack.main <- GeneRegionTrack(
  txdb,
  genome=assembly,
  showId = FALSE,
  name = "Gene Annotation"
) # Make it "gene" annotation instead of transcript annotation

# Unfortunately, txdb doesn't hold gene_symbol (name) so we have to get it from somewhere else
ref_dir<-dirname(gtf)
gtf_file<-basename(gtf)
gtf.file<-file.path(ref_dir, gtf_file)
sqlite_file<-gsub(".gtf", ".sqlite", gtf_file, fixed=T) # Name of the SQL file
sqlite_path <- file.path(ref_dir, sqlite_file)
if(!file.exists(sqlite_path)) {
  ensembldb::ensDbFromGtf(gtf=gtf.file, outfile=sqlite_path)
}
EnsDb <- ensembldb::EnsDb(sqlite_path)
#supportedFilters(EnsDb)

dir.create(outdir, recursive=T)

for(selectid in inlist$V1){
tryCatch({

  print(selectid)
  if(transcriptid == TRUE){
    parsedEnsembl <- transcripts(EnsDb,
                          filter = TxIdFilter(selectid),
                          columns = c("gene_id", "seq_name", "seq_length",
                                      "seq_strand", "tx_id", "tx_biotype",
                                      "tx_seq_start", "tx_seq_end", "gene_name"))

    # To get all transcripts for the gene, not just one
    parsedEnsembl <- transcripts(EnsDb,
                           filter = list(GeneIdFilter(parsedEnsembl$gene_id),
                                        TxBiotypeFilter(c("protein_coding", "lincRNA",
                                                            "antisense_RNA", "bidirectional_promoter_lncRNA",
                                                            "processed_transcript", "retained_intron",
                                                            "sense_intronic", "sense_overlapping"))),
                           columns = c("gene_id", "seq_name",
                                       "seq_strand", "tx_id", "tx_biotype",
                                       "tx_seq_start", "tx_seq_end", "gene_name"))
  }else if(geneid == TRUE){
    parsedEnsembl <- transcripts(EnsDb,
                           filter = list(GeneIdFilter(selectid),
                                        TxBiotypeFilter(c("protein_coding", "lincRNA",
                                                            "antisense_RNA", "bidirectional_promoter_lncRNA",
                                                            "processed_transcript", "retained_intron",
                                                            "sense_intronic", "sense_overlapping"))),
                           columns = c("gene_id", "seq_name",
                                       "seq_strand", "tx_id", "tx_biotype",
                                       "tx_seq_start", "tx_seq_end", "gene_name"))

  }else if(genename == TRUE){
    parsedEnsembl <- transcripts(EnsDb,
                           filter = list(SymbolFilter(selectid),
                                         TxBiotypeFilter(c("protein_coding", "lincRNA",
                                                            "antisense_RNA", "bidirectional_promoter_lncRNA",
                                                            "processed_transcript", "retained_intron",
                                                            "sense_intronic", "sense_overlapping"))),
                           columns = c("gene_id", "seq_name",
                                       "seq_strand", "tx_id", "tx_biotype",
                                       "tx_seq_start", "tx_seq_end", "gene_name"))
  }else{
    stop("Please set one of the transcriptid, geneid or genename to TRUE")
  }
  parsedEnsembl<-as.data.frame(parsedEnsembl)

  cov_len<-lengthOf(EnsDb, of="gene", filter = GeneIdFilter(unique(parsedEnsembl$gene_id))) # Length of transcript to get number of reads

  rand.l<-80
  rand.s<-cov_len/10 # how many random alignments we want to get from short read; 200 if we have Nanopore and Illumina in one; we use ALL the reads for coverage
  mychr<-as.character(unique(parsedEnsembl$seqnames))
  if(mychr=="MT"){
    mychr_ucsc<-"chrM"
  }else{
    mychr_ucsc<-paste0("chr", mychr)
    mychr<-as.numeric(mychr)
  }
  afrom<-min(parsedEnsembl$start)-extend
  ato<-max(parsedEnsembl$end)+extend
  astrand<-as.character(unique(parsedEnsembl$strand))
  # Set illumina flag for subsequent filtering
  # For reverse and + genes we take R2 reads, for reverse and - reads we take R1 reads,
  #   for forward and + we take R1 reads, for forward and - reads we take R2 reads
  if(strand.ill=="reverse" & astrand=="+"){
    aflag_R1<-83 # R1, reverse, first in pair, mapped in proper pair
    aflag_R2<-163 # R2, not reverse, mate reverse strand, second in pair, mapped in proper pair
  }else if(strand.ill=="reverse" & astrand=="-"){
    aflag_R1<-99 # R1, not reverse, mate reverse strand, first in pair, mapped in proper pair
    aflag_R2<-147 # R2, reverse, second in pair, mapped in proper pair
  }else if(strand.ill=="forward" & astrand=="+"){
    aflag_R1<-99 # R2, reverse, mate reverse strand, second in pair, mapped in proper pair
    aflag_R2<-147 # R2, reverse, second in pair, mapped in proper pair
  }else if(strand.ill=="forward" & astrand=="-"){
    aflag_R1<-83 # R1, reverse, first in pair, mapped in proper pair
    aflag_R2<-163 # R2, not reverse, mate reverse strand, second in pair, mapped in proper pair
  }else{
    stop("Please use [forward | reverse] for the strandedness of Illumina reads.")
  }

  # Select only gene we won't to visualize to avoid overlapping
  # Better this way than directly in GeneRegionTrack - more robust this way
  grtrack<-grtrack.main
  select_range <- grtrack@range[mcols(grtrack@range)$transcript %in% parsedEnsembl$tx_id]
  grtrack@range <- select_range
  grtrack@name<-unique(parsedEnsembl$gene_name)

  # Aesthetic adjustments
  displayPars(grtrack) <- list(col="transparent", fill=orange, background.title="transparent",
                                showTitle=FALSE, showId=FALSE, min.height=0, stackHeight=0.5, max.height=2, min.height=1) # Remove border around exons

  grtrack.d<-grtrack
  displayPars(grtrack.d)<-list(collapseTranscripts="meta") # Collapse "by gene"
  grtrack.d@name<-unique(parsedEnsembl$gene_name)
  grtrack.d@range$symbol<-grtrack.d@range$gene # Fake transcript id to fit to meta visualization; we could use the above as well

  bedannot<-bedannot.main[seqnames(bedannot.main) == mychr & start(bedannot.main) > afrom & end(bedannot.main) < ato & strand(bedannot.main) == astrand]
  atrack <- AnnotationTrack(bedannot, name = "APA", shape = "box",
                            stacking = "dense", fill = apa, col = apa,
                            background.title="transparent", showTitle=FALSE)

  options(ucscChromosomeNames = TRUE)
  itrack <- IdeogramTrack(
    genome = assembly, chromosome = mychr_ucsc,
    add53=TRUE, add35=TRUE
  )
  itrack@chromosome<-as.character(mychr) # A set of adjustments to fit to Ensembl naming instead of UCSC
  itrack@name<-as.character(mychr)
  levels(itrack@bandTable$chrom)[levels(itrack@bandTable$chrom)==mychr_ucsc] <- mychr
  options(ucscChromosomeNames = FALSE)

  displayPars(itrack)<-list(showAxis=FALSE)

  if(astrand=="+"){
    gtrack <- GenomeAxisTrack(add53=TRUE, add35=FALSE, littleTicks = FALSE, showId = TRUE)
  }else if(astrand=="-"){
    gtrack <- GenomeAxisTrack(add53=FALSE, add35=TRUE, littleTicks = FALSE, showId = TRUE)
  }

  # Load alignments and color by strand https://support.bioconductor.org/p/113503/
  # Note: If you get error: "Error in .normarg_at2(at, x) :
  #                           some ranges in 'at' are off-limits with respect to their corresponding
  #                           sequence in 'x'"
  #        and you have used Minimap2 it's most likely because Minimap2 doesn't output SEQ and QUAL for secondary alignments and GRanges doesn't like it.
  #       You have to get only primary alignments or use different aligner.
  readsTrackL<-Gviz:::.import.bam.alignments(inbam[1], GRanges(seqnames = mychr, ranges = IRanges(afrom, ato)))
  if(nchar(subset_read)>0){
    print("Found subsampling list, filtering long-reads")
    readsTrackL<-readsTrackL[readsTrackL$id %in% subset_read$V1]
  }
  readsTrackL<-readsTrackL[strand(readsTrackL) == astrand] # Further filter BAM GRanges to keep only correct strand mapping
  readsTrackL.all<-readsTrackL # For other but read visualizations
  if(length(unique(readsTrackL$id))>rand.l){
    print(paste0("To many aligned reads, subsampling long-reads to: ", rand.l))
    readnames<-sample(unique(readsTrackL$id), rand.l) # Get random X alignments to cover all possible alignments - full and short otherwise we have only long on the top
    readsTrackL<-readsTrackL[readsTrackL$id %in% readnames]
  }
  readsTrackL <- AlignmentsTrack(readsTrackL, genome=assembly, isPaired=FALSE, name="Nanopore", shape="box") # For reads LESS than the pdf size
  displayPars(readsTrackL)$fill.reads <- c(plus="#ECAFAE", minus="#AFAEED")[as.numeric(sort(ranges(readsTrackL)$readStrand))]
  readsTrackL.all <- AlignmentsTrack(readsTrackL.all, genome=assembly, isPaired=FALSE, name="Nanopore") # For reads LESS than the pdf size

  readsTrackS<-Gviz:::.import.bam.alignments(inbam[2], GRanges(seqnames = mychr, ranges = IRanges(afrom, ato)))
  readsTrackS<-readsTrackS[(readsTrackS$flag == aflag_R1 | readsTrackS$flag == aflag_R2) & readsTrackS$status == "mated"] # Further filter BAM GRanges to keep only reads with correct strand (filtered above) and SAM flag
  readsTrackS.all<-readsTrackS # For other but read visualizations
  if(length(unique(readsTrackS$id))>rand.s){
    print(paste0("To many aligned reads, subsampling short-reads to: ", rand.s))
    readnames<-sample(unique(readsTrackS$id), rand.s) # Get random X alignments to cover all possible alignments - full and short otherwise we have only long on the top
    readsTrackS<-readsTrackS[readsTrackS$id %in% readnames]
  }
  if(length(unique(readsTrackS.all$id))>opt$max_depth){
    print(paste0("Putting a cap to Illumina reads to: ", opt$max_depth))
    readnames<-sample(unique(readsTrackS.all$id), opt$max_depth) # Get random X alignments to cover all possible alignments - full and short otherwise we have only long on the top
    readsTrackS.all<-readsTrackS.all[readsTrackS.all$id %in% readnames]
  }
  readsTrackS <- AlignmentsTrack(readsTrackS, genome=assembly, isPaired=TRUE, name="Illumina", shape="box") # For reads LESS than the pdf size
  displayPars(readsTrackS)$fill.reads <- c(plus="#ECAFAE", minus="#AFAEED")[as.numeric(sort(ranges(readsTrackS)$readStrand))]
  readsTrackS.all <- AlignmentsTrack(readsTrackS.all, genome=assembly, isPaired=TRUE, name="Illumina") # For reads LESS than the pdf size

  # Minor aesthetic adjustments
  displayPars(readsTrackL)<-list(col="transparent", fill.coverage=long, fill.reads=long,
                                 col.axis="black", background.title="transparent", max.height=2, min.height=0)
  displayPars(readsTrackS)<-list(col="transparent", fill.coverage=short, fill.reads=short,
                                 col.axis="black", background.title="transparent", max.height=2, min.height=0)

  displayPars(readsTrackL.all)<-list(col=long)
  displayPars(readsTrackS.all)<-list(col=short)

  # Height of the bottom track
  displayPars(grtrack)<-list(size=length(unique(parsedEnsembl$tx_id))*0.2)

  pdf(paste0(outdir, "/", unique(parsedEnsembl$gene_id), "-alignments.pdf"))
          plotTracks(list(itrack, gtrack, grtrack.d, atrack, readsTrackL, readsTrackS, grtrack), # itrack
                   chromosome=mychr, from=afrom, to=ato,
                   main=paste0(unique(parsedEnsembl$gene_name), " (", unique(parsedEnsembl$gene_id), ")"),
                   sizes=c(0.5, 1, 0.3, 0.15, 7, 7, 2),
                   coverageHeight = 0.08, minCoverageHeight = 0,
                   cex.main=1, cex.title=1, cex.group=0.4,
                   min.width=0,
                   lwd=0.5, lwd.reads=0.5,
                   fontcolor.title="black", fontcolor="black", fontcolor.group="black") # We could remove specific regions

        # Check sashimi plots
        plotTracks(list(itrack, gtrack, grtrack.d, atrack, readsTrackL.all, readsTrackS.all, grtrack),
                   chromosome=mychr, from=afrom, to=ato, type = c("coverage", "sashimi"),
                   sashimiFilterTolerance = 5L, extend.left=500-extend,
                   cex.group=0.4,
                   min.width=0,
                   fontcolor.title="black", fontcolor="black", fontcolor.group="black") # We could remove specific regions
      dev.off()
    }, error=function(e){})
  graphics.off()
}