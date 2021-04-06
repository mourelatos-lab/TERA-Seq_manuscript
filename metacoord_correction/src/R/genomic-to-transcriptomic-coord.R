#!/usr/bin/env Rscript
#
# Conversion of genomic to transcriptomic coordinates
#
# main source: https://www.biostars.org/p/251076/
# for data.frame: https://www.biostars.org/p/230758/
#

library(optparse)
library(tibble)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicFeatures))
library(rtracklayer)

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Input genomic coord BED to be converted to transcriptomic coord. Can be stdin (use 'stdin')", metavar = "File"),
  make_option(
    c("-a", "--annot"), type = "character",
    help = "Gene annotation GTF to use for the conversion.", metavar = "File"),
  make_option(
    c("-l", "--left"), type = "integer", metavar = "Number", default=0,
    help = "Add N nucleotides to the first exon (extend to the \"left\", upstream. Default: 0, no extension."),
  make_option(
    c("-r", "--right"), type = "integer", metavar = "Number", default=0,
    help = "Add N nucleotides to the last exon (extend to the \"right\", downstream. Default: 0, no extension."),       
  make_option(
    c("-o", "--ofile"), type = "character",
    help = "Output genomic to transcriptomic converted BED.", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

### TESTING VARIABLES
# print("USING TESTING VARIABLES!!!")
#opt<-NULL
#opt$ifile<-"test.bed"
#opt$ofile<-"testout.bed"
#opt$annot<-"../../references/Homo_sapiens.GRCh38.91.sorted.gtf"
#opt$left<-10
#opt$right<-0
# print("USING TESTING VARIABLES!!!")
# ### TESTING VARIABLES

##################################################################################
if(opt$ifile == "stdin"){
	genomicLocus <- import.bed(con=file("stdin"))
}else{
	genomicLocus <- import.bed(con=opt$ifile)
}

gencodeTxDb <- makeTxDbFromGFF(file=opt$annot)

# Get coordinates by transcript including introns
gencodeTx <- transcripts(gencodeTxDb)

# Extend annotation left/right
if(opt$left>0 | opt$right>0){
  gencodeTxPlus<-gencodeTx[strand(gencodeTx) == "+"]
  gencodeTxMinus<-gencodeTx[strand(gencodeTx) == "-"]
  
  if(opt$left>0){
    start(gencodeTxPlus)<-start(gencodeTxPlus)-opt$left
    start(gencodeTxPlus)[start(gencodeTxPlus)<1]<-1 # Check for chromosome overflow
  }
  if(opt$right>0){
    end(gencodeTxMinus)<-end(gencodeTxMinus)+opt$right
  }
  gencodeTx<-append(gencodeTxPlus, gencodeTxMinus)
}

names(genomicLocus)<-genomicLocus$name
names(gencodeTx) <- id2name(gencodeTxDb, "tx")
mappedLocusTrans <- mapToTranscripts(x=genomicLocus, transcripts=gencodeTx, ignore.strand = F)

rtracklayer::export(object = mappedLocusTrans, gsub(".bed", "-ByTranscript.bed", opt$ofile), format = "bed")

# Keep only common chromosomes (based on gencodeTxDb)
# This can be an issue if we don't have extra chromsomes in both groups and will result in error
genomicLocus<-genomicLocus[seqnames(genomicLocus) %in% seqlevels(gencodeTxDb)]
seqlevels(genomicLocus) <- as.vector(unique(seqnames(genomicLocus)@values)) # remove unused levels
seqlevels(gencodeTxDb) <- levels(seqnames(genomicLocus)) # remove unused chromosomes from input annotation

# Get coordinates by exons excluding introns
gencodeExons <- exonsBy(gencodeTxDb, "tx", use.names=TRUE)

# Function to extend start/end with accounting for strand information
extendExons<-function(grang, index, extStart = 0, extEnd = 0){
  # Get easy accessible strand column
  strands<-as_tibble(grang) %>% 
    tibble::remove_rownames() %>%
    as.data.frame() %>% 
    dplyr::select(strand)
  strands<-strands$strand
  
  if(length(grang)!=length(index)){
    stop("GRanges object and index do not have the same length, please check your index.")
  }
  
  if(opt$left>0){ # Extend start of the exon
    start(grang[index & (strands == "+")])<-start(grang[index & (strands == "+")])-extStart
    end(grang[index & (strands == "-")])<-end(grang[index & (strands == "-")])+extStart
  }
  if(opt$right>0){ # Extend end of the exon
    end(grang[index & (strands == "+")])<-end(grang[index & (strands == "+")])+extEnd
    start(grang[index & (strands == "-")])<-start(grang[index & (strands == "-")])-extEnd
  }
  return(grang)
}

# Extend annotation left/right
if(opt$left>0 | opt$right>0){
  gencodeExonsUnl <- unlist(gencodeExons, use.names = TRUE)
  gencodeExonsUnl$transcript_id<-names(gencodeExonsUnl)
  
  # Some exnos (same exon_id) are shared among multiple transcripts and they will be duplicated since they could be first exon in one trans, and Nth in another
  singleExIdx<-names(gencodeExonsUnl) %in% names(table(names(gencodeExonsUnl))[table(names(gencodeExonsUnl))==1]) # List of single exon transcripts/ids

  firstExIdx<-gencodeExonsUnl$exon_rank==1 # Get list of first exons

  lastExIdx<-as_tibble(gencodeExonsUnl) %>% 
    tibble::remove_rownames() %>% group_by(transcript_id) %>% mutate(last_ex = max(exon_rank)) %>% 
    as.data.frame()
  lastExIdx<-lastExIdx$exon_rank==lastExIdx$last_ex

  # Extend exons  
  gencodeExonsUnl<-extendExons(grang = gencodeExonsUnl, index = singleExIdx, extStart = opt$left, extEnd = opt$right)
  gencodeExonsUnl<-extendExons(grang = gencodeExonsUnl, index = firstExIdx, extStart = opt$left, extEnd = opt$right)
  gencodeExonsUnl<-extendExons(grang = gencodeExonsUnl, index = lastExIdx, extStart = opt$left, extEnd = opt$right)

  start(gencodeExonsUnl)[start(gencodeExonsUnl)<1]<-1 # Check for chromosome overflow

  gencodeExonsExt<-relist(gencodeExonsUnl, gencodeExons)
  gencodeExons<-gencodeExonsExt
}

mappedLocusExon <- mapToTranscripts(x=genomicLocus, transcripts=gencodeExons, ignore.strand = F)

rtracklayer::export(object = mappedLocusExon, gsub(".bed", "-ByExon.bed", opt$ofile), format = "bed")

sink(file = gsub(".bed", "-check.txt", opt$ofile))
  print("Exon coordinates of ENST00000616125.")
  gencodeExons$ENST00000616125 # These are the exon coordinates
  print("Genomic input bed coordinates of ENST00000616125 first exons.")
  genomicLocus[seqnames(genomicLocus) == "1" & start(genomicLocus) > 925740 & end(genomicLocus) < 944581 & strand(genomicLocus) == "+"] # Check coordinates for one transcript, this is the genomic coordinates to map
  print("Mapped input bed coordinates in ENST00000616125 exons (without introns).")
  mappedLocusExon[seqnames(mappedLocusExon)=="ENST00000616125"] # These are the mapped coordinates in exons
  print("Mapped input bed coordinates in ENST00000616125 transcript (with introns).")
  mappedLocusTrans[seqnames(mappedLocusTrans)=='ENST00000616125']  # These are the mapped coordinates in transcript (with introns most likely!!!)
sink()
