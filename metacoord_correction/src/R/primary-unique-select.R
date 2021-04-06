#!/usr/bin/env Rscript
#
# Make the same "primary" (count each position only once) or unique
#   (~unambiguous, only non-overlapping transcripts positions) positions from bed file
# The selection is done on 4th column - the name!
#

library("optparse")
suppressPackageStartupMessages(library("data.table"))

option_list <- list(
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Input file (BED format). The selection is done on 4th column (name).", metavar = "File"),
  make_option(
    c("-t", "--type"), type = "character", default = "primary",
    help = "Type of selection [primary|unique]. Default: primary", metavar = "String"),
  make_option(
    c("-o", "--ofile"), type = "character", default = "stdout",
    help = "Output file (bed), Default: stdout.", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

### TESTING VARIABLES ###
#echo("!!! USING TESTING VARIABLES !!!")
#opt <- NULL
#opt$ifile <- "/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/common/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.bed"
#opt$type <- "primary"
#opt$ofile <- "/home/joppelt/projects/ribothrypsis/analysis/distro-on-meta-gene/results/common/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.primary.bed"
#echo("!!! USING TESTING VARIABLES !!!")
### TESTING VARIABLES ###

# isUnique function from biobase
isUnique <- function(x) {
  rv <- rep(TRUE, length(x))

  if (length(x) >= 2) {
    ord <- order(x)
    ox <- x[ord]
    ## compare consecutive values
    neq <- (ox[-length(ox)] != ox[-1])
    ## a value is unique if neither its predecessor nor successor
    ## in the ordered vector are the same
    rv[ord] <- c(neq, TRUE) & c(TRUE, neq)
  }

  return(rv)
}

# read bed
bed <- fread(opt$ifile, sep = "\t")

# if primary - choose unique and random position if name is duplicated
# if unambing - choose unique only
uniq <- isUnique(bed$V4)
#!uniq are all duplicated

otab <- bed[uniq]

if(opt$type != "unique"){
  # Assign random number to non-unique positions
  dupl <- bed[!uniq]
  setDT(dupl)[, random := sample(.N, .N), V4]
  dupl<-dupl[random==1] # Get only primary (=randomly chosen position)
  dupl[, random:=NULL]
  
  otab <- rbindlist(list(otab, dupl), use.names=TRUE)
}

setorderv(otab, c("V1", "V2", "V3")) # sort

if(opt$ofile == "stdout"){
  fwrite(x = otab, file = stdout(), sep = "\t", col.names = F)
}else{
  fwrite(x = otab, file = opt$ofile, sep = "\t", col.names = F)
}