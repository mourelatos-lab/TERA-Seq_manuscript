#!/usr/bin/env Rscript
#
# Integrate Manolis corrections to the original annotation
# Manolis "corrected" the coordinates on some of the transcipts based on the real observation but didn't use this information to improve the orignal annotation
#
# We can call this as a "bed update" script which takes an original bed and updates whatever is in the new bed and saves out the "old with correction" bed
#

library(optparse)
suppressMessages(library(dplyr))
library(ggplot2)

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ifile5p"), type="character",
    help="Coordinate for the replacement of the 5' end"),
  make_option(
    c("-i2", "--ifile3p"), type="character",
    help="Coordinate for the replacement of the 3' end"),
  make_option(
    c("-b", "--bed"), type="character",
    help="Coordinates where the replacement should happen (expecting bed with at least 4 columns)"),
  make_option(
    c("-c1", "--column_ifile"), type="character", default = 1,
    help="Number of the matching column in ifile (expecting bed with at least 3 columns)"),
  make_option(
    c("-c2", "--column_bed"), type="character", default = 1,
    help="Number of the matching column in bed"),
  make_option(
    c("-o", "--ofile"), type="character",
    help="Output the updated file (bed with original bed+1 columns)")
)
opt = parse_args(OptionParser(option_list = option_list))

### Testing variables
 # opt<-NULL
 # opt$ifile5p<-"/home/joppelt/projects/ribothrypsis/analysis/distro-on-meta-gene/results/coords-corrected/new-mrna-coords.5p.primary.tab" # Corrected coordinates
 # opt$ifile3p<-"/home/joppelt/projects/ribothrypsis/analysis/distro-on-meta-gene/results/coords-corrected/new-mrna-coords.3p.primary.tab" # Corrected coordinates
 # opt$bed<-"/home/joppelt/projects/ribothrypsis/data/hg38/genic_elements.mrna.bed" # Where to correct the coordinates
 # opt$column_ifile<-1 # What is the common feature (transcript,gene,...)
 # opt$column_bed<-1 # What is the common feature (transcript,gene,...)
 # opt$ofile<-"/home/joppelt/projects/ribothrypsis/analysis/distro-on-meta-gene/results/coords-corrected/genic_elements.mrna.rel5.corrected.bed"
### Testing variables

#new_coor<-read.delim(opt$ifile5p, header=F, stringsAsFactors = F)
new_coor_5p<-read.delim(opt$ifile5p, header=F, stringsAsFactors = F)
new_coor_3p<-read.delim(opt$ifile3p, header=F, stringsAsFactors = F)
bed<-read.delim(opt$bed, header=F, stringsAsFactors = F)

no_col<-ncol(bed)

# Merge and select 5p and 3p coordinates
new_coor_5p$V3<-NA # annulate end position
new_coor_3p$V2<-NA # annulate start position

new_coor<-new_coor_5p %>% select(c("V1", "V2")) %>%
  full_join(select(new_coor_3p, c("V1", "V3")), by = "V1")

# Merge - dplyr cannot do by.x/by.y
bed<-bed %>%
  left_join(new_coor, by = c(colnames(new_coor)[as.numeric(opt$column_ifile)], colnames(bed)[as.numeric(opt$column_bed)]))

# Add note about changing the coordinates (also testing if we use the correct condition for the next part)
bed$corr<-""
bed[(!(is.na(bed$V2.y)) & (bed$V2.x != bed$V2.y)), "corr"]<-"5p" # Add a note
bed[(!(is.na(bed$V3.y)) & (bed$V3.x != bed$V3.y)), "corr"]<-paste0(bed[(!(is.na(bed$V3.y)) & (bed$V3.x != bed$V3.y)), "corr"], "3p") # Add a note
bed$corr[bed$corr==""]<-"orig"

# Replace if not NA with dplyr
bed <- bed %>%
  mutate(V2.x = ifelse((!is.na(V2.y) & (V2.x != V2.y)), V2.y, V2.x)) %>%
  mutate(V3.x = ifelse((!is.na(V3.y) & (V3.x != V3.y)), V3.y, V3.x))

# Make bed great again
bed<-bed[, c(1:no_col, ncol(bed))]

# Adjust length to reflect the corrections
bed[,5]<-bed[,3]-bed[,2]

# Summary of number of changes
table(bed$corr)

# Export
write.table(x = bed, file = opt$ofile, sep="\t", row.names=F, col.names=F, quote = F)

## Get some statistics about the number of correction and the length of transcript
# Remove very short
bed <- bed %>%
  filter(V5 >= 200)

# Get basic stats
bed$len_bin<-cut(bed$V5, breaks=seq(0, max(bed$V5)+200, by=200), labels=seq(200, max(bed$V5)+200, by=200))

corr_stats<-as.data.frame(table(bed[, c("corr", "len_bin")]))
corr_stats$len_bin<-as.numeric(as.character(corr_stats$len_bin))
corr_stats$Freq<-as.numeric(as.character(corr_stats$Freq))
colnames(corr_stats)<-c("Ends", "Length", "Count")

# Add the total sum
corr_stats<-corr_stats %>%
  group_by(Length) %>%
  summarize(Count = sum(Count)) %>%
  mutate(Ends="total") %>%
  bind_rows(corr_stats)

corr_stats$Ends<-factor(corr_stats$Ends, levels=unique(corr_stats$Ends)) # Put back the factors

print(corr_stats)

#pdf(paste0(sub("^(.*)[.].*", "\\1", opt$ofile), "-stats.pdf"), height=6, width = 8)
#  ggplot(data=corr_stats, aes(x=Length, y=Count, colour=Ends)) +
#    geom_line() +
#    xlim(0, 15000)
#
#  # Cut xaxis and show just a zoom of one section https://stackoverflow.com/questions/7194688/using-ggplot2-can-i-insert-a-break-in-the-axis
#  ggplot() +
#    aes(x = Length, y = Count, colour=Ends) +
#    geom_line(data = corr_stats %>% mutate(subset = "all")) +
#    geom_line(data = corr_stats %>% filter(Count <= 200) %>% mutate(subset = "short")) +
#    facet_wrap(~ subset, scales = "free_y") +
#    xlim(0, 5000) +
#	xlab("Transcript length (bins by 200 nt)")
#dev.off()
