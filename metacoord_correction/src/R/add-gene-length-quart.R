#!/usr/bin/env Rscript
#
# Add gene lengths and divide them by quartiles
#

suppressMessages(library("ensembldb"))
suppressMessages(library("dplyr"))
library("readr")
library("optparse")

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Input table", metavar = "File"),
  make_option(
    c("-g", "--gtf"), type = "character",
    help = "Path to a Ensembl gtf file (please you the original name which includes organism and annotation version).", metavar = "File"),
  make_option(
    c("-o", "--ofile"), type = "character",
    help = "Output tsv file", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

### Testing variables
#opt<-NULL
#opt$gtf<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/data/hg38/Homo_sapiens.GRCh38.91.gtf.gz"
#opt$ifile<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/distro-on-meta-gene/results/hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1/coords-on-corrected-meta-mrnas.rel5.tab"
#opt$ofile<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/distro-on-meta-gene/results/hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1/coords-on-corrected-meta-mrnas-inclLen.rel5.tab"
### Testing variables

ref_dir <- dirname(opt$gtf)
gtf_file <- basename(opt$gtf)

### Parse GTF 
gtf.file<-file.path(ref_dir, gtf_file)
sqlite_file<-gsub(".gtf", ".sqlite", gtf_file, fixed=T) # Name of the SQL file
sqlite_path <- file.path(ref_dir, sqlite_file)

# This part sometimes does not put sqlite database to the 'path', needs to be transfered manually
if(!file.exists(sqlite_path)) {
  ## generate the SQLite database file
  ensembldb::ensDbFromGtf(gtf=gtf.file, outfile=sqlite_path)
}
EnsDb <- ensembldb::EnsDb(sqlite_path)

### Read the input table(s)
df<-read_tsv(opt$ifile)

# Get lengths only for the detected transcripts (much faster)
len <- transcriptLengths(EnsDb, 
                         with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE, 
                         filter=TxIdFilter(df$feat))

len<-len %>%
  filter(tx_len!=0) # Remove zero lengths

df.temp <- df %>%
  select(feat) %>%
  unique()

df.temp <- len %>% 
  rename(feat = tx_id) %>%
  right_join(df.temp) %>% 
  mutate(quartile = ntile(tx_len, 4)) # Add quartiles based on length

df <- df.temp %>%
  select(feat, quartile, tx_len) %>%
  right_join(df)

# Keep only the quartiles and tx length as the additional information
df <- df %>%
  dplyr::select(group1, qname, feat, bin5p, bin3p, tx_len, quartile)

write.table(df, file=opt$ofile, quote=FALSE, sep='\t', row.names=FALSE)

# Get stats
len_stats <- df %>%
  group_by(quartile) %>%
  summarize(Read_count=n(), Tx_count=length(unique(feat)), Tx_mean_len=mean(tx_len), Tx_med_len=median(tx_len))

# Remove everything after last "." and replace https://stackoverflow.com/questions/17847189/remove-text-after-final-period-in-string
write.table(len_stats, file=paste0(sub("^(.*)[.].*", "\\1", opt$ofile), "_quartStats.tsv"), quote=FALSE, sep='\t', row.names=FALSE)
