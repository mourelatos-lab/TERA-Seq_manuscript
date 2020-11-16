#!/usr/bin/env Rscript
#
# Experimental plotting of alignment using Gviz (https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html)
# Book: Statistical Genomics; Ewy Math? ? Sean Davis Editors; Methods and Protocols
# Tutorials: https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/viz_course/Presentations/singlepage/Viz_part_1.html
#   https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/viz_course/Presentations/singlepage/Viz_part_2.html; https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/
#
# Note: This script is organism specific due to biomart inputs
# Note: This script CANNOT do read subset, you have to subset the input BAM
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
    c("-s", "--name"), type = "character", metavar = "String",
    help = "Sample name to use in the plot."),
  make_option(
    c("-a", "--assembly"), type = "character", default="hg38", metavar = "String",
    help = "Genome reference assembly [hg38 | mm10]"),
  make_option(
    c("-g", "--gtf"), type = "character", metavar = "GTF file",
    help = "Ensembl GTF. Has to have the GTF name."),
  make_option(
      c("-f", "--fasta"), type = "character", metavar = "FASTA file",
    help = "Reference FASTA file used for alignment."),
  make_option(
    c("-b", "--bedcage"), type = "character", metavar = "BED file",
    help = "Additional BED annotation. For examples alternative CAGE BED."),
  make_option(
    c("-b2", "--bedcage2"), type = "character", metavar = "BED file",
    help = "Additional BED annotation. For examples alternative NET-CAGE BED."),
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

### TESTING VARIABLES!!!
# # print("TESTING VARIABLES!!!")
# opt<-NULL
# # Minimap2 doesn't output SEQ (and QUAL) for secondary alignments and GRanges doesn't like it
# opt$nanopore<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/alignment/reads.1.sanitize.toGenome.sorted.primary.wo-softclip.w-polya.bam" # samtools view -@ 26 -F4 -F256 -bh ~/projects/ribothrypsis/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/align/reads.1.sanitize.toGenome.sorted.bam > CIP.decap.reads.1.sanitize.toGenome.sorted.bam; samtools index CIP.decap.reads.1.sanitize.toGenome.sorted.bam
# opt$name<-"CIP.decap"
# opt$assembly<-"hg38"
# opt$gtf<-"/home/joppelt/projects/ribothrypsis/data/hg38/Homo_sapiens.GRCh38.91.gtf.gz"
# opt$fasta<-"/home/joppelt/projects/ribothrypsis/data/hg38/genome/genome.fa"
# # opt$fasta<-"/home/jan/playground/ref/hg38/genome.fa"
# opt$bedcage<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/common/HeLa.repAll.hg38.ctss.bed"
# opt$bedcage2<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/common/RepAll-HeLaS3-NETCAGE-1M.hg38.bed"
# opt$list<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/transcripts_select.manual2.txt"
# opt$transcriptid<-TRUE
# opt$geneid<-FALSE
# opt$genename<-FALSE
# opt$extend<-50
# opt$outdir<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/transcript-coverage/results/test"
#print("TESTING VARIABLES!!!")
### TESTING VARIABLES!!!

inbam<-opt$nanopore # nanopore, nanopore2
assembly<-opt$assembly
gtf<-opt$gtf
fasta<-opt$fasta
inbed<-opt$bedcage
inbed2<-opt$bedcage2
inlist<-opt$list # One transcript id/gene id/gene name per line
transcriptid<-opt$transcriptid
geneid<-opt$geneid
genename<-opt$genename
extend<-opt$extend

# Colors
orange<-rgb(250, 164, 26, maxColorValue = 255)
long<-"#F8766D"
long2<-"#00BFC4"
cage<-"black"
cage2<-"violet"

##########################################
inlist<-read.table(inlist, header=F)

# Make TxDb from GTF
txdb<-makeTxDbFromGFF(gtf, format="auto", organism="Homo sapiens", dataSource="Ensembl")

# Load extra BED annotations
bedannot.main <- rtracklayer::import(inbed)
bedannot.main2 <- rtracklayer::import(inbed2)

# IMPORTANT: We need to set this to false otherwise Gviz won't show any alignments and will mess up the gene visualization
# This sest ALL chr reference for Gviz to be "1, 2, 3" and not "chr1, chr2, chr3"
# https://support.bioconductor.org/p/78337/
options(ucscChromosomeNames = FALSE)

### Get gene annotation from Ensemlb (or any other GTF)
grtrack.main <- GeneRegionTrack(
  txdb,
  #  chromosome = mychr, start = afrom-50, end = ato+50,
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

# Import genomic fasta
genome <- Biostrings::readDNAStringSet(fasta, "fasta")

dir.create(opt$outdir, recursive = T)

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
      # parsedEnsembl <- transcripts(EnsDb,
      #                        filter = list(GeneIdFilter(parsedEnsembl$gene_id),
      #                                     TxBiotypeFilter(c("protein_coding", "lincRNA",
      #                                                         "antisense_RNA", "bidirectional_promoter_lncRNA",
      #                                                         "processed_transcript", "retained_intron",
      #                                                         "sense_intronic", "sense_overlapping"))),
      #                        columns = c("gene_id", "seq_name",
      #                                    "seq_strand", "tx_id", "tx_biotype",
      #                                    "tx_seq_start", "tx_seq_end", "gene_name"))
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

    rand.l<-100
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

    # Get fasta sequence track
    if(astrand=="+"){
      strack <- SequenceTrack(sequence = genome, chromosome = mychr, name=assembly, add53 = TRUE)
    }else if(astrand=="-"){
      strack <- SequenceTrack(sequence = genome, chromosome = mychr, name=assembly, add35 = TRUE)
    }

    fcol <- c(A="#2C7BB6", C="#FDAE61", T="#ABD9E9", G="#D7191C", N="#FFFFBF") # Change colors of all letters but A to highlight polyA tail
    displayPars(strack)<-list(fontcolor=fcol)

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
    atrack <- AnnotationTrack(bedannot, name = "CAGE", shape = "box",
                              stacking = "dense", fill = cage, col = cage,
                              background.title="transparent", showTitle=FALSE)

    bedannot2<-bedannot.main2[seqnames(bedannot.main2) == mychr & start(bedannot.main2) > afrom & end(bedannot.main2) < ato & strand(bedannot.main2) == astrand]
    atrack2 <- AnnotationTrack(bedannot2, name = "NET-CAGE", shape = "box",
                               stacking = "dense", fill = cage2, col = cage2,
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

    readsTrackL_simp<-AlignmentsTrack(inbam, genome=assembly, isPaired=FALSE, name=opt$name, shape="box", referenceSequence = strack)
    displayPars(readsTrackL_simp)<-list(showMismatches = TRUE,
										col="transparent", background.title="transparent", col.axis="black",
										col.border.title="transparent", fontcolor.title="black",
                                        max.height=2, min.height=0) 

    # Height of the bottom track
    displayPars(grtrack)<-list(size=length(unique(parsedEnsembl$tx_id))*0.2)

    # Flipping of the visualization in Gviz
    if(astrand=="+"){
      reverseStrand<-FALSE
      complement<-FALSE
    }else if(astrand=="-"){
      reverseStrand<-TRUE
      complement<-TRUE
    }

    pdf(paste0(opt$outdir, "/", unique(parsedEnsembl$gene_id), "-alignments.pdf"))
    plotTracks(list(itrack, gtrack, grtrack.d, atrack, atrack2, readsTrackL_simp, grtrack), # itrack
               chromosome=mychr, from=afrom, to=ato,
               main=paste0(unique(parsedEnsembl$gene_name), " (", unique(parsedEnsembl$gene_id), ")"),
               sizes=c(0.5, 1, 0.75, 0.15, 0.15, 7, 0.75),
               coverageHeight = 0.08, minCoverageHeight = 0,
               reverseStrand = reverseStrand,
               cex.main=1, cex.title=1, cex.group=0.4,
               min.width=0,
               lwd=0.5, lwd.reads=0.5)
    dev.off()
  }, error=function(e){})
  graphics.off()
}
