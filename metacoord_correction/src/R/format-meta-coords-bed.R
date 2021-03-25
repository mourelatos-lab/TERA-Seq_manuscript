#!/usr/bin/env Rscript
#
# Format 6-column bed (CAGE/polyA sites, for example) so it fits to the plot-meta-coords.R input data structure
#
# qname	feat	bin5p	bin3p
# dd809ebb-ccb7-41fd-8aa9-31925449756c	ENST00000342066	6	19
#
# Note: This is SOOO inefficient but we'll run it only once (or twice) so it's ok
#

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("rio"))
suppressPackageStartupMessages(library("ggplot2"))
library("optparse")

option_list <- list(
  make_option(
    c("-t", "--trans"), type = "character",
    help = "Transcript annotation (BED format).", metavar = "File"),
  make_option(
    c("-b", "--bed"), type = "character",
    help = "CAGE coordinates in bed format (with header).", metavar = "File"),
  make_option(
    c("-m", "--len"), type = "integer", default = 200,
    help = "Minimum length of transcript to output. Default: 200.", metavar = "Integer"),
  make_option(
    c("-o", "--ofile"), type = "character", default = "stdout",
    help = "Output file (tsv), Default: stdout.", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

### TESTING VARIABLES !!!
#print("USING TESTING VARIABLES!!!")
# opt<-NULL
# opt$trans<-"/home/joppelt/projects/ribothrypsis/analysis/distro-on-meta-gene/results/coords-corrected/genic_elements.mrna.unambig.tab"
# #opt$trans<-"/home/joppelt/projects/ribothrypsis/data/hg38/genic_elements.mrna.bed"
# opt$bed<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/common/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.bed"
# opt$len<-200
# #opt$ofile<-"/home/joppelt/projects/ribothrypsis/analysis/distro-on-meta-gene/results/CAGE/coords-on-corrected-meta-mrnas.tab"
# opt$ofile<-"stdout"
#print("USING TESTING VARIABLES!!!")
### TESTING VARIABLES !!!

# isUnique function from Biobase https://rdrr.io/bioc/Biobase/src/R/tools.R
isUnique = function(x){

  rv = rep(TRUE, length(x))

  if(length(x)>=2) {
    ord = order(x)
    ox = x[ord]
    ## compare consecutive values
    neq = (ox[-length(ox)]!=ox[-1])
    ## a value is unique if neither its predecessor nor successor
    ## in the ordered vector are the same
    rv[ord] = c(neq, TRUE) & c(TRUE, neq)
  }

  return(rv)
}

################################################################################

trans<-rio::import(opt$trans, format="tsv")
colnames(trans)<-c("chr", "start", "end", "name", "score", "strand")
trans<-trans[,!is.na(colnames(trans))] # remove extra columns if the bed file didn't have 6 column

intab_read<-rio::import(opt$bed, format="tsv")
colnames(intab_read)<-c("chr", "start", "end", "name", "score", "strand")

# Add length of the transcript
trans$len <- trans$end - trans$start

intab_read <- intab_read %>%
  left_join(trans, by="chr")

# remove unmapped transcripts and get only positions within the annotation
intab_read<-intab_read %>%
  filter(start.x >= start.y) %>%
  filter(start.x <= end.y) %>%
  filter(end.x <= end.y) %>%
  filter(end.x >= start.y) %>%
  filter(len >= !!opt$len) %>%
  select(c("chr", ends_with(".x")))
colnames(intab_read)<-gsub(pattern = ".x", replacement = "", x = colnames(intab_read), fixed = T)

trans <- trans %>%
  filter(chr %in% unique(intab_read$chr))

tab<-NULL
counter<-1
for(chr_name in unique(trans$chr)){
  trans.t <- trans %>%
    filter(trans$chr %in% chr_name)

  cuts <- cut(seq(trans.t[1, "start"], trans.t[1, "end"], by = 1), breaks = 20, labels = seq(0, 19, by=1), include.lowest = T)
  cuts<-as.data.frame(cuts)
  cuts$cuts <- as.numeric(as.character(cuts$cuts))
  cuts$start<-seq(trans.t[1, "start"], trans.t[1, "end"], by = 1)
  cuts$chr <- chr_name

  intab_read.t <- intab_read %>%
    filter(intab_read$chr %in% chr_name)

  intab_read.t <- intab_read.t %>%
    left_join(cuts, by = c("chr", "start")) %>%
    rename(bin5p = cuts)

  if(is.null(tab)){
    tab<-intab_read.t
  }else{
    tab<-rbind(tab, intab_read.t)
  }

  if(counter %% 1000 == 0){
    message(paste("Processed", counter))
  }
  counter=counter+1
}

tab$bin3p<-tab$bin5p

# Rename if name vector doesn't have unique names
if(sum(isUnique(tab$name))==0){
  tab$name <- paste(tab$chr, seq(1, nrow(tab)), sep=".")
  tab$name <- paste(tab$chr, tab$start, sep=".")
}

tab <- tab[,c("name", "chr", "bin5p", "bin3p")]
colnames(tab) <- c("qname", "feat", "bin5p", "bin3p")

if(opt$ofile == "stdout"){
  rio::export(x = tab, file = "", format = "tsv")
}else{
  dir.create(dirname(opt$ofile), recursive = T)
  rio::export(x = tab, file = opt$ofile, format = "tsv")
}
