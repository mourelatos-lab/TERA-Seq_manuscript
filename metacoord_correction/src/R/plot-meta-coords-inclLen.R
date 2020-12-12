#!/usr/bin/env Rscript
#
# Author: Manolis
# Here changed to confirm that the Minion doesn't prefer shorter reads over longer ones (added tx_len and quartile to the input file)
#

library(optparse)
library(ggplot2)
suppressMessages(library(data.table))
library(grid)
library(gtable)
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
#library(scales)

# dens calculates the density for numeric vector x.
dens <- function(x) {
	z = x/sum(x)
	return(z)
}

# quant calculates the quantile value at pct for the non-zero values of numeric
# vector x.
quant <- function(x, pct) {
	q = quantile(x[x > 0], probs = pct)[1]
	return(q)
}

# Get colors as Fadia likes them - wo adapter - black, rel5/w adapter - green, rel3 - blue
color.rel5=rgb(91, 186, 71,maxColorValue = 255) # REL5/w adapter=Green; # R:91; # G:186; # B:71
color.wAd=color.rel5 # REL5/w adapter=Green; # R:91; # G:186; # B:71
color.rel3=rgb(52, 105, 179, maxColorValue = 255) # REL3=Blue; R:52; # G:105; # B:179
color.woAd=rgb(0, 0, 0, maxColorValue = 255) # wo adapter=Black; # R:0; # G:0; # B:0

# Read command line options and arguments
option_list <- list(
	make_option(
	    c("-i", "--ifile"), type = "character",
	    help = "Input table", metavar = "File"),
	make_option(
	    c("-s", "--subset"), type = "double", default=1.0,
	    help = "Percent of top expressed transcripts to keep [0-1]. Default: 1 (keep all)", metavar = "Num"),
	make_option(
            c("-f", "--figfile"), type = "character",
	    help = "Output pdf file with graphs"),
	make_option(
	    c("-r", "--recalquart"), action = "store_true", default = FALSE,
	    help = "Recalculate quartiles (usually for merged inputs)"),
	make_option(
	    c("-a", "--addtransnum"), action = "store_true", default = FALSE,
	    help = "Add number of transcripts to the main plot. Will fail for multiple groups!")
)
opt = parse_args(OptionParser(option_list = option_list))

### Testing variables
# print("USING TESTING VARIABLES!!!")
# opt<-NULL
# opt$ifile<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/distro-on-meta-gene/results/hsa.dRNASeq.HeLa.polyA.REL5OH.long.1/coords-on-corrected-meta-mrnas.unambig.norel5-inclLen.tab"
# opt$subset<-0.3
# opt$recalquart<-TRUE
# opt$addtransnum<-TRUE
# opt$figfile<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/distro-on-meta-gene/results/hsa.dRNASeq.HeLa.polyA.REL5OH.long.1/test.pdf"
# print("USING TESTING VARIABLES!!!")
### Testing variables

qprob = 1-as.numeric(opt$subset)

# Create data frame from file.
df = read.delim(opt$ifile)

if(opt$recalquart == TRUE){
  print("Recalculating quartiles.")

# Split to equal n-tiles - 4=quartiles
#  df <- df %>%
#    select(feat, tx_len) %>%
#    unique() %>%
#    mutate(quartile = ntile(tx_len, 4)) # Add quartiles based on length

# Split by 500 bp
  df$quartile <- cut(df$tx_len, seq(0, max(df$tx_len)+500, by=500), labels = seq(500, max(df$tx_len)+500, by=500),
                          include.lowest = TRUE) # Add bins by 500 bp

  len_stats <- df %>%
    group_by(quartile) %>%
    summarize(Read_count=n(), Tx_count=length(unique(feat)), Tx_mean_len=mean(tx_len), Tx_med_len=median(tx_len))

  # Export recalculated quartile stats
  write.table(len_stats, file=paste0(sub("^(.*)[.].*", "\\1", opt$figfile), "_quartStats.tsv"), quote=FALSE, sep='\t', row.names=FALSE)
}else{
  print("Keeping quartiles as they are.")
}

# Create data table.
dt = data.table(df)

#dt<-dt[sample(nrow(dt), 50000), ] # Subset for 50000 reads

# Calculate meta-length. Meta-length is the difference of the 3' and 5' bin.
dt[, metalen := bin3p - bin5p]

# Filter data to keep the most highly expressed transcripts.
# 1. Measure the number of reads for each transcript.
# 2. Find a limit for the counts to filter lowly expressed transcripts.
# 3. Remove lowly-expressed from dt by merging the expressed ones.
# 4. Remove extra columns.
dt.tr.count = dt[, .(count = .N), by = .(group1, feat)]
dt.tr.count[, q := quant(count, qprob), by = .(group1)]
dt.tr.count.expressed = dt.tr.count[count >= q]
dt = merge(dt, dt.tr.count.expressed, by = c("group1", "feat"))
dt[, c("count", "q") := NULL]

# Calculate meta-length distribution. Meta-length is the difference of the 3'
# and 5' bins.
dt.metalen = dt[, .(count = .N), by = .(group1, quartile, metalen)]
dt.metalen = dt.metalen[order(group1, quartile, metalen),]
dt.metalen[, percent := (count/sum(count))*100, by = .(group1, quartile)]
dt.metalen[, cumpercent := cumsum(percent), by = .(group1, quartile)]

# Melt dt to have long format for 5' and 3' bins.
dt.melt = melt(dt, id.vars = c("group1", "feat", "qname", "metalen", "quartile", "tx_len"), variable.name = "type", value.name = "bin")
dt.melt$type <- factor(dt.melt$type, levels = c("bin5p", "bin3p"))

# Count reads in each bin of each transcript.
dt.bin.counts = dt.melt[, .(count = .N), by = .(group1, feat, type, bin, quartile)]

# Bins without reads are missing from the table. Add them with 0s.
# https://stackoverflow.com/questions/22462405/add-missing-rows-to-a-data-table
dt.bin.counts = dt.bin.counts[CJ(group1 = group1, feat = feat, type = type, bin = bin, quartile=quartile, unique = TRUE), on = .(group1, feat, type, bin, quartile)][is.na(count), count := 0L][]

# Count percent count over total library in each bin in each quartile.
dt.bin.counts[, percent := 100*count/sum(count), by = .(group1, quartile, type)]

# Measure read density in each bin per transcript in each quartile.
dt.bin.counts[, dens := percent/sum(percent), by = .(group1, quartile, feat, type)]

# Filling-in missing bins has a sideffect and adds extra transcripts that
# exist in one group1 but not another. Remove these transcripts to avoid NAs.
dt.bin.counts = na.omit(dt.bin.counts)

# Aggregate transcripts in each quartile.
dt.bin.counts.transcr_agg = dt.bin.counts[, .(
				percent_sum = sum(percent),
				count_sum = sum(count),
        count_mean = mean(count),
        count_sd = sd(count),
        dens_mean = mean(dens),
        dens_sd = sd(dens),
				n_transcripts = paste0("Tn=", uniqueN(feat)))
            , by = .(group1, quartile, type, bin)]
dt.bin.counts.transcr_agg[, n_reads := paste0("Rn=", sum(count_sum)), by = .(group1, quartile, type)] # Add read number; we have to to it by type (3p and 5p) otherwise we count each read twice

dt.bin.counts.transcr_agg$group1<-gsub("hsa.dRNASeq.HeLa.", "", dt.bin.counts.transcr_agg$group1) # Just some renaming
dt.bin.counts.transcr_agg$group1<-factor(dt.bin.counts.transcr_agg$group1, levels=unique(dt.bin.counts.transcr_agg$group1))

# Add transcript and read number
if(opt$addtransnum==TRUE){
  dt.bin.counts.transcr_agg$quartile_name<-paste0(dt.bin.counts.transcr_agg$quartile,
                                                  "\n(", dt.bin.counts.transcr_agg$n_transcripts, "\n", dt.bin.counts.transcr_agg$n_reads, ")")
}else{
  dt.bin.counts.transcr_agg$quartile_name<-dt.bin.counts.transcr_agg$quartile
}

dt.bin.counts.transcr_agg$quartile_name<-factor(dt.bin.counts.transcr_agg$quartile_name,
                                                levels=unique(dt.bin.counts.transcr_agg$quartile_name[order(dt.bin.counts.transcr_agg$quartile)])) # Make sure the levels are sorted by the bin length

### Assign colors for Fadia
# For dt.bin.counts.transcr_agg
adapt.g<-unique(levels(dt.bin.counts.transcr_agg$group1)[grep(pattern = "-rel5", x = levels(dt.bin.counts.transcr_agg$group1))])
noadapt.g<-unique(levels(dt.bin.counts.transcr_agg$group1)[grep(pattern = "-norel5", x = levels(dt.bin.counts.transcr_agg$group1))])

if(length(adapt.g)==1 & length(noadapt.g)==1){ # colors if we have both w/wo rel5
  cols<-c("bin5p" = color.wAd, "bin3p" = color.woAd)
  names(cols)[1]<-adapt.g
  names(cols)[2]<-noadapt.g
}else if(length(adapt.g)==1 & length(noadapt.g)==0){ # colors if we have w rel5 only
  cols<-c("bin5p" = color.wAd)
  names(cols)[1]<-adapt.g
}else if(length(adapt.g)==0 & length(noadapt.g)==1){ # colors if we have wo rel5 only
  cols<-c("bin3p" = color.woAd)
  names(cols)[1]<-noadapt.g
}else{ # Default ggplot colors https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  cols<-scales::hue_pal()(length(levels(dt.bin.counts.transcr_agg$group1))) # See the colors: scales::show_col(scales::hue_pal()(length(levels(dt.bin.counts.transcr_agg$group1))))
  names(cols)<-levels(dt.bin.counts.transcr_agg$group1)
}
# For dt.metalen
adapt.g<-unique(levels(dt.metalen$group1)[grep(pattern = "-rel5", x = levels(dt.metalen$group1))])
noadapt.g<-unique(levels(dt.metalen$group1)[grep(pattern = "-norel5", x = levels(dt.metalen$group1))])

if(length(adapt.g)==1 & length(noadapt.g)==1){ # colors if we have both w/wo rel5
  cols.2<-c("bin5p" = color.wAd, "bin3p" = color.woAd)
  names(cols.2)[1]<-adapt.g
  names(cols.2)[2]<-noadapt.g
}else if(length(adapt.g)==1 & length(noadapt.g)==0){ # colors if we have w rel5 only
  cols.2<-c("bin5p" = color.wAd)
  names(cols.2)[1]<-adapt.g
}else if(length(adapt.g)==0 & length(noadapt.g)==1){ # colors if we have wo rel5 only
  cols.2<-c("bin3p" = color.woAd)
  names(cols.2)[1]<-noadapt.g
}else{ # Default ggplot colors https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  cols.2<-scales::hue_pal()(length(levels(dt.metalen$group1))) # See the colors: scales::show_col(scales::hue_pal()(length(levels(dt.bin.counts.transcr_agg$group1))))
  names(cols.2)<-levels(dt.metalen$group1)
}

# Create plots
pdf(opt$figfile, width = 9)
  ggplot(data = subset(dt.bin.counts.transcr_agg, quartile %in% c(seq(500, 3000, by=500))),
         aes(x = bin, y = percent_sum)) +
  	geom_line(aes(colour = group1)) +
  	coord_cartesian(ylim = c(0, 100)) +
  	scale_y_continuous(breaks = seq(0, 100, by=10), labels = c(0, "", 20, "", 40, "", 60, "", 80, "", 100)) +
  	facet_grid(quartile_name ~ type) +
  	ylab("read count (%)") +
    scale_colour_manual(values = cols) +
  	theme_classic() +
    guides(color = guide_legend(ncol=2)) +
    theme(legend.position="bottom")

  ggplot(data = subset(dt.metalen, quartile %in% c(seq(500, 3000, by=500))),
         aes(x = metalen, y = percent, color = group1)) +
      geom_line() +
      geom_point() +
      facet_grid(. ~ quartile) +
      xlab("meta-length (bins)") +
      ylab("read count (%)") +
      scale_colour_manual(values = cols.2) +
      theme_classic() +
    guides(color = guide_legend(ncol=2)) +
    theme(legend.position="bottom")

  ggplot(data = subset(dt.metalen, quartile %in% c(seq(500, 3000, by=500))),
         aes(x = metalen, y = cumpercent, color = group1)) +
      geom_line() +
      geom_point() +
      facet_grid(. ~ quartile) +
      xlab("meta-length (bins)") +
      ylab("cumulative read count (%)") +
      scale_colour_manual(values = cols.2) +
      theme_classic() +
    guides(color = guide_legend(ncol=2)) +
    theme(legend.position="bottom")

  ggplot(data = subset(dt.bin.counts.transcr_agg, quartile %in% c(seq(500, 3000, by=500))),
         aes(x = bin, y = dens_mean, ymin = dens_mean - dens_sd, ymax = dens_mean + dens_sd)) +
  	geom_line(aes(colour = group1)) +
    coord_cartesian(ylim = c(0, 1.0)) +
    scale_y_continuous(breaks = seq(0.0, 1.0, by=0.1), labels = c(0, "", 0.2, "", 0.4, "", 0.6, "", 0.8, "", 1.0)) +
    geom_ribbon(aes(fill = group1), alpha = 0.1) +
  	facet_grid(quartile ~ type) +
  	ylab("average read density (avg. across transcripts)") +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
  	theme_classic() +
    guides(color = guide_legend(ncol=2), fill = guide_legend(ncol=2)) +
    theme(legend.position="bottom")
graphics.off()
