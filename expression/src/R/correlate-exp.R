#!/usr/bin/env Rscript
#
# Correlate two or more expression inputs
#
# Expects input to be tab separated with following columns: lib_name\ttranscript_id\tTPM\tfeature1\tfeature2\tfeature2\tfeatureN
#   lib_name, transcript, expression column names should be as stated or at least in the same order
#   feature1/2/3/... can have any column name and should be numeric
#

library("optparse")
library("data.table")
library("ggplot2")
library("plyr")
library("corrplot")

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Input table", metavar = "File"),
  make_option(
    c("-l", "--exp_limit"), type = "numeric", default = 999999,
    help = "Maximum expression for the correlation plots (rows with >exp_limit genes are removed)", metavar = "integer"),
  make_option(
    c("-f", "--figfile"), type = "character",
    help = "Output pdf file", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

################################################################################
exp_table<-fread(opt$ifile)

exp_table<-exp_table[TPM>0] # Also remove 0 expressed transcripts - most likely from the Salmon quantification as it reports all transcript in the annotation and not just the detected ones

head(exp_table)
unique(exp_table$lib_name)

### Expressions
exp_stats<-exp_table[, .(
  expressed_tx = sum(TPM>0),
  expressed_tx_tpm_1 = sum(TPM>1),
  average_tx_exp = mean(TPM),
  median_tx_exp = median(TPM),
  max_tx_exp = max(TPM)
), by = .(lib_name)]

fwrite(exp_stats, paste0(gsub(".pdf", "", opt$figfile), ".exp_stats.tsv"), sep = "\t")

exp_table<-exp_table[!grep("cytosolic|nuclear", exp_table$lib_name, perl = T),] # Temporary remove cytosolic and nuclear Illumina libraries

exp_table$lib_type<-sub('.*\\.', '', exp_table$lib_name) # Get point types for the plot - expects lib_name with ".xxx" at the end where ".xxx" will be the point type
if(length(grep(".RNASeq.", exp_table$lib_name, fixed=T))>0){
    exp_table[grep(".RNASeq.", exp_table$lib_name, fixed = T)]$lib_type<-paste(exp_table[grep(".RNASeq.", exp_table$lib_name, fixed = T)]$lib_type, "ILL", sep=".")
}
exp_table$lib_name<-factor(exp_table$lib_name, levels = unique(exp_table$lib_name)) # Fix factors
exp_table$lib_type<-factor(exp_table$lib_type, levels = unique(exp_table$lib_type)) # Fix factors

min_tpm<-min(exp_table$TPM)
exp_table$log2TPM<-log2(exp_table$TPM+min_tpm)

movmed <- ddply(exp_table, .(lib_name), summarise, med = median(log2TPM)) # Median points for the libraries
movmed$lib_type<-sub('.*\\.', '', movmed$lib_name) # Get point types for the plot - expects lib_name with ".xxx" at the end where ".xxx" will be the point type
if(length(grep(".RNASeq.", movmed$lib_name, fixed=T))>0){
    movmed[grep(".RNASeq.", movmed$lib_name, fixed = T),]$lib_type<-paste(movmed[grep(".RNASeq.", movmed$lib_name, fixed = T),]$lib_type, "ILL", sep=".")
}
movmed$lib_type<-factor(movmed$lib_type, levels = unique(movmed$lib_type)) # Fix factors

### Correlations
# Put a limit on expression is required
exp_table.cor<-exp_table
exp_table.cor <- exp_table.cor[, TPM:=as.double(TPM)]

if(opt$exp_limit!=999999){
  print(paste("Putting limit on expression for the correlation of:", opt$exp_limit))
  exp_table.cor[exp_table.cor$TPM>opt$exp_limit, "TPM"]<-NA
}else{
  print("Using all expression levels for the correlation plots.")
}

corr_table.in<-dcast(data = exp_table.cor, formula = transcript_id+cdna_len~lib_name, value.var = "log2TPM")
rownames(corr_table.in)<-corr_table.in$transcript_id
corr_table.in$transcript_id<-NULL

corr_table.in$cdna_len<-NULL
corr_table.in<-corr_table.in[rowSums(corr_table.in)!=0] # Remove empty rows
corr_table.in<-corr_table.in[rowSums(is.na(corr_table.in))<ncol(corr_table.in)] # Remove NA rows

### Correlation plots
## Pearson
corr_table.in.bckp<-corr_table.in
corr_table.in<-corr_table.in[rowSums(corr_table.in)!=0] # Remove empty rows
corr_table.in<-corr_table.in[rowSums(is.na(corr_table.in))<ncol(corr_table.in)] # Remove NA rows

Mp.pairwise<-cor(corr_table.in, method="pearson", use = "na.or.complete") # Spearman = rank, pearson = value

pdf(paste0(gsub(".pdf", "", opt$figfile), ".samples_correlation-pearson.pdf"), width = 12, height = 10)
  corrplot(Mp.pairwise, addCoef.col = "black", number.digits = 2, number.cex = 2,
           type = "upper", cl.pos = "r", tl.col = "black",  cl.cex = 1.2, tl.cex = 1.2, tl.srt = 45,
           main=paste0("\nPearson - all vs. all (exp. limit: " , opt$exp_limit, ")"))
  corrplot(Mp.pairwise, method="number", main=paste0("\nPearson - all vs. all (exp. limit: " , opt$exp_limit, ")"))
dev.off()

### Correlation plots
corr_table.in<-corr_table.in.bckp

Mp.pairwise<-cor(corr_table.in, method="spearman", use = "na.or.complete") # Spearman = rank, pearson = value

pdf(paste0(gsub(".pdf", "", opt$figfile), ".samples_correlation-spearman.pdf"), width = 12, height = 10)
  # Publication read plot
  corrplot(Mp.pairwise, method="circle", addCoef.col = "black", number.digits = 2, number.cex = 2,
           type = "upper", cl.pos = "r", tl.col = "black",  cl.cex = 1.2, tl.cex = 1.2, tl.srt = 45,
           main=paste0("\nSpearman - all vs. all (exp. limit: " , opt$exp_limit, ")"))
  corrplot(Mp.pairwise, method="number", main=paste0("\nSpearman - all vs. all (exp. limit: " , opt$exp_limit, ")"))
dev.off()

### Do Pearson with Illumina
corr_table.in<-corr_table.in.bckp

Mp.pairwise<-cor(corr_table.in, method="pearson", use = "na.or.complete") # Spearman = rank, pearson = value

### Scatter plots
# Correlation panel; http://www.sthda.com/english/wiki/scatter-plot-matrices-r-base-graphs
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method="pearson", use = "na.or.complete"), digits=2)
  txt <- paste0("R = ", r)
  text(0.5, 0.5, txt, cex = 5)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = "black", cex=0.2) # col = my_cols[iris$Species
}

# Create the plots

corr_table.in<-corr_table.in.bckp

if(length(grep("REL5", colnames(corr_table.in)))>1){
  corr_table.in.rel5<-corr_table.in[, grep("REL5", colnames(corr_table.in)), with=FALSE]

  png(paste0(gsub(".pdf", "", opt$figfile), ".samples_pairs-REL5.png"), width = 1200, height = 1200)
  pairs(corr_table.in.rel5,
        lower.panel = panel.cor,
        upper.panel = upper.panel,
        cex.labels = 4,
        main="Pearson correlation and pairwise plots")
  dev.off()
}
if(length(grep("REL3", colnames(corr_table.in)))>1){
  corr_table.in.rel3<-corr_table.in[, grep("REL3", colnames(corr_table.in)), with=FALSE]

  png(paste0(gsub(".pdf", "", opt$figfile), ".samples_pairs-REL3.png"), width = 1200, height = 1200)
  pairs(corr_table.in.rel3,
        lower.panel = panel.cor,
        upper.panel = upper.panel,
        cex.labels = 5,
        main="Pearson correlation and pairwise plots")
  dev.off()
}
