#!/usr/bin/env Rscript
#
# Originaly designed for Nanopolish, here modified to tailfindr
#
# TODO: Complete merge of nanopolish and tailfindr to one script, rigth now it's somewhere in between
#

library(optparse)
library(ggplot2)
library(data.table)
#library(dplyr) # Loading later

# Read command line options and arguments
option_list <- list(
    make_option(
        c("-i", "--ifile"), type = "character",
        help = "Input table", metavar = "File"),
    make_option(
      c("-s", "--software"), type = "character", default = "tailfindr",
      help = "Which tools was used to make the results [nanopolish; tailfindr]"),
    make_option(
	    c("-l", "--maxL"), type = "integer", default = 300,
	    help = "Maximum poly(A) length to plot", metavar = "Num"),
    make_option(
        c("-f", "--figfile"), type = "character",
        help = "Output pdf file with graphs"),
    make_option(
      c("-r", "--removezero"), action = "store_true", default = FALSE,
      help = "Remove NA or zero values")
)
opt = parse_args(OptionParser(option_list = option_list))

### Testing variables
# print("USING TEST VARIABLES!!!")
# opt<-NULL
# opt$ifile<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/polya-len/results/hsa.dRNASeq.HeLa.total.REL3.2/nanopolish/rna_tails.rel3.tsv"
# opt$software<-"nanopolish"
#opt$removezero<-TRUE
# opt$maxL<-300
# opt$figfile<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/polya-len/results/hsa.dRNASeq.HeLa.total.REL3.2/nanopolish/test.pdf"
# print("USING TEST VARIABLES!!!")
### Testing variables

dt = data.table(read.delim(opt$ifile)) # At the end we converted tailfindr to tsv as well but the default output is csv

# Create data table from file.
if(opt$software=="nanopolish"){
#    dt = data.table(read.delim(opt$ifile))
    dt = dt[qc_tag == "PASS" | qc_tag == "NOREGION"] # Filter rows.
  }else if(opt$software=="tailfindr"){
#    dt = data.table(read.delim(opt$ifile))
#    dt = dt[!(is.na(tail_start)) & !(is.na(tail_end)) & !(is.na(tail_length))] # Filter rows (remove NA) - most likely removes 0-length tails!
    setnames(dt, "tail_length", "polya_length") # Rename tailfindr just to fit to Nanopolish naming
  }else{
    print("Unknown input format, please specify either nanopolish or tailfindr")
}

dt[is.na(dt)]<-0 # Replace all undetected by 0

dt$polya_length<-round(dt$polya_length, 0) # Round values to even numbers

# Get percentage of having/not having polyA tail before we do anything with it
polya_stats.tmp <- dt[, .(
  num_polya = sum(polya_length>=5),
  num_nonpolya = sum(polya_length<5)
  ), by = .(sample)]

polya_stats.tmp$perc_polya<-round((polya_stats.tmp$num_polya/rowSums(polya_stats.tmp[, c("num_polya","num_nonpolya")])*100), 2)

if(opt$removezero == TRUE){
  print("Removing zero (or NA) lengths.")
  dt <- dt[polya_length!=0]
}else{
  print("Keeping zero (or NA) lengths.")
}

# Before we do the length filtering, make boxplots
polya_box<-ggplot(dt, aes(y = polya_length, x = sample, colour = sample)) +
  geom_boxplot() +
  xlab("Libraries") +
  ylab("Poly(A) length (nt)") +
  theme_classic()

# Get statistics
polya_stats <- dt[, .(
  len_max = max(polya_length),
  len_min = min(polya_length),
  len_median = median(polya_length),
  len_mean = round(mean(polya_length), 2),
  len_sd = round(sd(polya_length), 2)
  ), by = .(sample)]

polya_stats <- merge(polya_stats, polya_stats.tmp, all.x=TRUE) # Add perc. of polyA

write.table(polya_stats, file=paste0(sub("^(.*)[.].*", "\\1", opt$figfile), ".tsv"), quote=FALSE, sep='\t', row.names=FALSE)

# Filter rows within poly-A length limit.
if (opt$maxL > 0) {
#	dt = dt[polya_length <= opt$maxL]
    dt$polya_length[dt$polya_length > opt$maxL] <- opt$maxL
}

change_col<-F # Default is to use "common" colors
# Get colors as Fadia likes them
# REL5=Green:
# R:91
# G:186
# B:71
rel5=rgb(91, 186, 71,maxColorValue = 255)
# REL3=Blue
# R:52
# G:105
# B:179
rel3=rgb(52, 105, 179, maxColorValue = 255)

# Change colors only if we have two samples and one of the is REL5 and the other one is REL3 library...uff
if(length(unique(dt$sample))==2 & length(grep("REL5", levels(dt$sample), ignore.case = F))==1 & length(grep("REL3", levels(dt$sample), ignore.case = F))==1){
  change_col<-TRUE # Switch to custom colors
  cols<-c("REL5" = rel5, "REL3" = rel3)
  names(cols)[names(cols)=="REL5"]<-grep("REL5", levels(dt$sample), ignore.case = F, value = T)
  names(cols)[names(cols)=="REL3"]<-grep("REL3", levels(dt$sample), ignore.case = F, value = T)
}

# For ultimate.x sample vs (cip).decap
if(length(unique(dt$sample))==2 & length(grep("decap.REL5", levels(dt$sample), ignore.case = F))==1 & length(grep("total.REL5.long.REL3", levels(dt$sample), ignore.case = F))==1){
  change_col<-TRUE # Switch to custom colors
  cols<-c("REL5" = rel5, "REL3" = rel3)
  names(cols)[names(cols)=="REL5"]<-grep("decap.REL5", levels(dt$sample), ignore.case = F, value = T)
  names(cols)[names(cols)=="REL3"]<-grep("total.REL5.long.REL3", levels(dt$sample), ignore.case = F, value = T)
}

# For REL3 only
if(length(grep("REL3", levels(dt$sample), ignore.case = F))==1 & length(grep("REL3", levels(dt$sample), ignore.case = T, invert = T))==0){
  change_col<-TRUE # Switch to custom colors
  cols<-c("REL3" = rel3)
  names(cols)[names(cols)=="REL3"]<-grep("REL3", levels(dt$sample), ignore.case = F, value = T)
}

# Create plots
pdf(opt$figfile, width = 7)
  if(change_col==TRUE){
    a<-ggplot(dt, aes(x = polya_length, stat(density), color = sample)) +
    	geom_freqpoly(binwidth = 1, alpha = 0.50) +
      scale_colour_manual(values = cols) +
      scale_fill_manual(values = cols) +
    	xlab("Poly(A) length (nt)") +
    	theme_classic()
  }else{
    a<-ggplot(dt, aes(x = polya_length, stat(density), color = sample)) +
      geom_freqpoly(binwidth = 1, alpha = 0.50) +
      xlab("Poly(A) length (nt)") +
      theme_classic()
  }
  print(a)

  if(change_col==TRUE){
   a<-ggplot(dt, aes(x=polya_length, color = sample, fill=sample)) +
     geom_histogram(aes(y=..density..), position="identity", binwidth=10, alpha=.25) + # For ggplot >=3.0.0
     geom_vline(data=polya_stats, aes(xintercept=len_median, color=sample),
                linetype="dashed") +
    scale_colour_manual(values = cols) +
     scale_fill_manual(values = cols) +
     theme(legend.position="top") +
     xlab("Poly(A) length (nt)")  +
     ylab("Density")  +
     ggtitle("Poly(A) length (bin 10 nt) densities\n with median") +
     theme_classic()
  }else{
    a<-ggplot(dt, aes(x=polya_length, color = sample, fill=sample)) +
      geom_histogram(aes(y=..density..), position="identity", binwidth=10, alpha=.25) + # For ggplot >=3.0.0
      geom_vline(data=polya_stats, aes(xintercept=len_median, color=sample),
                 linetype="dashed") +
      theme(legend.position="top") +
      xlab("Poly(A) length (nt)")  +
      ylab("Density")  +
      ggtitle("Poly(A) length (bin 10 nt) densities\n with median") +
      theme_classic()
  }
  print(a)

  if(change_col==TRUE){
    ggplot(dt, aes(x=polya_length, color = sample, fill=sample)) +
      geom_histogram(position="identity", binwidth=10, alpha=.25)  +
      geom_vline(data=polya_stats, aes(xintercept=len_median, color=sample),
                 linetype="dashed") +
      scale_colour_manual(values = cols) +
      scale_fill_manual(values = cols) +
      theme(legend.position="top") +
      xlab("Poly(A) length (nt)")  +
      ggtitle("Poly(A) length (bin 10 nt) total counts\n with median") +
      theme_classic()
  }else{
    ggplot(dt, aes(x=polya_length, color = sample, fill=sample)) +
      geom_histogram(position="identity", binwidth=10, alpha=.25)  +
      geom_vline(data=polya_stats, aes(xintercept=len_median, color=sample),
                 linetype="dashed") +
      theme(legend.position="top") +
      xlab("Poly(A) length (nt)")  +
      ggtitle("Poly(A) length (bin 10 nt) total counts\n with median") +
      theme_classic()
  }

  if(change_col==TRUE){
    print(polya_box + scale_colour_manual(values = cols))
  }else{
    print(polya_box)
  }
graphics.off()
