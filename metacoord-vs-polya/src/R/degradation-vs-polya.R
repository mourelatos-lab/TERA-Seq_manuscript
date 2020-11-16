#!/usr/bin/env Rscript
#
# Compare polyA length with different aspects, mainly with degradation
#

library(optparse)
suppressMessages(library(data.table))
library(ggplot2)
suppressMessages(library(reshape2))
suppressMessages(library(ggpubr))
suppressMessages(library(plyr))

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Input table", metavar = "File"),
  make_option(
    c("-g", "--tx_list"), type = "character",
    help = "List of transcripts for which to do separate plots - one transcript per line without any header, or tsv separate file where transcript is the first column and again no header", metavar = "File"),
  make_option(
    c("-o", "--outdir"), type = "character",
    help = "Output directory for the plots. Multiple plots will be generated.")
)

opt = parse_args(OptionParser(option_list = option_list))

### Test variables
# print("USING TEST VARIABLES!!!")
# opt<-NULL
# opt$ifile<-"/home/joppelt/projects/ribothrypsis/analysis/degradation-levels/results/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/coords-on-meta-mrnas.unambig.rel5-inclPolya.tab"
# opt$tx_list<-"/home/joppelt/projects/ribothrypsis/data/hg38/HeLa_transcripts.all-perGene-one2one.transcript-gene.txt"
# opt$outdir<-paste(dirname(opt$ifile), "test", sep="/")
# print("USING TEST VARIABLES!!!")
# ### Test variables

# Get colors
# REL5=Green:
# R:91; G:186; B:71
rel5=rgb(91, 186, 71,maxColorValue = 255)
# REL3=Blue
# R:52; G:105; B:179
rel3=rgb(52, 105, 179, maxColorValue = 255)

tx_list<-read.table(opt$tx_list, sep="\t", header=F, stringsAsFactors = F)$V1

dt<-fread(opt$ifile)

dt = dt[qc_tag == "PASS" | qc_tag == "NOREGION"] # Filter rows.

## Plot polyA vs degradation
# Make polyA bins
b <- seq(0, max(dt$polya_length)+25, by=25)
dt$polya_bin<-cut(dt$polya_length, 
                           breaks = b, 
                           labels = seq(25, max(dt$polya_length)+25, by=25), 
                           include.lowest = T)

# Get polya stats
# If the corr. test fails it gets 999.9 assigned 
polya_stats<-dt[, .(
  reads = .N,
  num_polya = sum(polya_length>=5),
  num_nonpolya = sum(polya_length<5),
  max_polya = max(polya_length, na.rm = T),
  min_polya = min(polya_length, na.rm = T),
  polya_length.median = median(polya_length, na.rm = T),
  polya_length.mean = mean(polya_length, na.rm = T),
  polya_length.sd = sd(polya_length, na.rm = T)
), by = .(feat)]

# Add RPM
polya_stats[, RPM := ((reads+1) / sum(reads+1)) * 10^6]

# How many transcripts have median length of polyA smaller than 1?
print("This many transcripts has polyA length median <1:")
nrow(polya_stats[polya_stats$polya_length.median < 1, "feat"])

# bin5p should be negative correlation
print("Total features to measure the correlation:")
nrow(polya_stats)

## Make plot expression vs polyA length
# Make bins for polyA
b <- seq(0, max(polya_stats$polya_length.median)+25, by=25)
polya_stats$polya_bin<-cut(polya_stats$polya_length.median, 
                           breaks = b, 
                           labels = seq(25, max(polya_stats$polya_length.median)+25, by=25), 
                           include.lowest = T)
# Make bins for expression
polya_stats$RPM.log10<-log10(polya_stats$RPM+1)
c <- seq(0, max(polya_stats$RPM.log10)+0.1, by=0.1)
polya_stats$rpm_bin<-cut(polya_stats$RPM.log10, 
                         breaks = c, 
                         labels = seq(0.1, max(polya_stats$RPM.log10)+0.1, by=0.1), 
                         include.lowest = T)

## Heatmaps
dt.heatmap<-dt

# Merge with the main table for heatmaps (and other plots)
dt.heatmap <- merge(dt.heatmap, polya_stats[, c("feat", "reads", "RPM")], all.x = T)

## Subset data.table
# Get top X most expressed genes
top_genes<-polya_stats[order(-polya_stats$RPM)[1:nrow(polya_stats)], "feat"][1:20] # top 20 genes by expression

# Make the list shorter otherwise it doesn't look nice in the plot - get only random 50 transcripts
tx_list<-tx_list[tx_list %in% subset(dt.heatmap, reads > 10)$feat] # Get only transcripts which are in the list and in the data and have more than 10 reads
if(length(tx_list)>50){
  print("tx_list is too long, randomly subsetting to 50 transcripts.")
  tx_list<-tx_list[sample(1:length(tx_list), 50, replace=F)]
}

dt.heatmap.subset<-dt.heatmap[feat %in% unique(c(top_genes$feat, tx_list))] # most expressed + transcripts in a list
# Sort the data.table by the top20 at the beginning and then the rest of the transcripts
setDT(dt.heatmap.subset)[ , feat := factor(feat, levels = unique(c(top_genes$feat, tx_list)))] # https://stackoverflow.com/questions/11977102/order-data-frame-rows-according-to-vector-with-specific-order
dt.heatmap.subset<-dt.heatmap.subset[with(dt.heatmap.subset, order(feat)), ] # https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-multiple-columns

# Make a new order per group (transcript) - expression, transcript name, completness of 5p, completness of 3p
sort.index <- dt.heatmap.subset[order(-dt.heatmap.subset$RPM, -dt.heatmap.subset$bin5p, dt.heatmap.subset$bin3p), "qname"] # dt.heatmap.subset$feat, 

dt.heatmap.subset.m <- reshape2::melt(dt.heatmap.subset, measure.vars=c("bin5p", "bin3p")) # Parially melt the table

dt.heatmap.subset.m$value<-factor(dt.heatmap.subset.m$value, levels = seq(0, 19, by = 1)) # Make the bin factors to keep them in the plot

### Plots
## Heatmap - this works fine only for a few transcripts
dt.heatmap.subset.m$feat<-factor(dt.heatmap.subset.m$feat, levels=unique(c(top_genes$feat, tx_list))) # Make factors

dt.heatmap.subset.m.bckp<-dt.heatmap.subset.m

pdf(file = paste(opt$outdir, paste0(sub(".[^.]+$", "", basename(opt$ifile)), "-transcript-heatmap-polya.pdf"), sep="/"), onefile=T)
  for(tx_name in unique(dt.heatmap.subset.m.bckp$feat)){

    dt.heatmap.subset.m<-dt.heatmap.subset.m.bckp[feat == tx_name]

    # Change colors only if we have two samples and one of the is REL5 and the other one is REL3 library...uff
    if(length(grep("bin5p", levels(dt.heatmap.subset.m$variable), ignore.case = F))==1 & length(grep("bin3p", levels(dt.heatmap.subset.m$variable), ignore.case = F))==1){
      cols<-c("bin5p" = rel5, "bin3p" = rel3)
      names(cols)[names(cols)=="bin5p"]<-grep("bin5p", levels(dt.heatmap.subset.m$variable), ignore.case = F, value = T)
      names(cols)[names(cols)=="bin3p"]<-grep("bin3p", levels(dt.heatmap.subset.m$variable), ignore.case = F, value = T)
    }
      
# If we want separately 5p and 3p plot next to each other. Also usefull if we want more transcripts in a single plot and not each of them separately
    # main.plot.5p<-ggplot(subset(dt.heatmap.subset.m, variable %in% "bin5p"), 
    #   aes(x = value, y = factor(qname, levels = sort.index$qname), color = feat, fill = feat)) + # 
    #   geom_tile() + 
    #   geom_tile(subset(dt.heatmap.subset.m, variable %in% "bin3p"), mapping = aes(x = value, y = factor(qname, levels = sort.index$qname), 
    #                                                                               fill = feat, color=feat), alpha=0.1) + 
    #   scale_x_discrete(breaks = seq(0, 19, by=1), drop=FALSE) + # Keep empty factors in the plot https://stackoverflow.com/questions/10834382/ggplot2-keep-unused-levels-barplot
    #   theme(
    #     panel.grid.minor=element_blank(),
    #     panel.grid.major=element_blank(),
    #     axis.text.y = element_blank()
    #   )

    # main.plot.3p<-ggplot(subset(dt.heatmap.subset.m, variable %in% "bin3p"), 
    #   aes(x = value, y = factor(qname, levels = sort.index$qname), color = feat, fill = feat)) + # 
    #   geom_tile() + 
    #   geom_tile(subset(dt.heatmap.subset.m, variable %in% "bin5p"), mapping = aes(x = value, y = factor(qname, levels = sort.index$qname), 
    #                                                                               fill = feat), alpha=0.1) +
    #   scale_x_discrete(breaks = seq(0, 19, by=1), drop=FALSE) + # Keep empty factors in the plot https://stackoverflow.com/questions/10834382/ggplot2-keep-unused-levels-barplot
    #   theme(
    #     panel.grid.minor=element_blank(),
    #     panel.grid.major=element_blank(),
    #     axis.line.y = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.title.y = element_blank(), 
    #     axis.ticks.y = element_blank()
    #   )

# If we want both 5p and 3p plot in the same plot and we have one transcript per plot
     main.plot<-ggplot(dt.heatmap.subset.m, 
       aes(x = value, y = factor(qname, levels = sort.index$qname), color = variable, fill = variable)) + # 
       geom_tile() + 
       scale_x_discrete(breaks = seq(0, 19, by=1), drop=FALSE) + # Keep empty factors in the plot https://stackoverflow.com/questions/10834382/ggplot2-keep-unused-levels-barplot
       xlab("Degradation bin") +
       ylab("Reads") +
       labs(fill = "Read end", color = "Read end") +
       scale_colour_manual(values = cols) +
       scale_fill_manual(values = cols) +
       theme_classic() + 
       theme(
         panel.grid.minor=element_blank(),
         panel.grid.major=element_blank(),
         axis.text.y = element_blank()
       )

    dt.heatmap.subset.m$polya_length[dt.heatmap.subset.m$polya_length>300]<-300 # Cheat and make a cap at the max. polyA length - using ylim causes NOT to show values over the limit making them 0 which we don't want; we can do the same with closing the xlim/ylim to coord_cartesian() but it throws an error
  
    polya.bar <- 
      ggplot(subset(dt.heatmap.subset.m, variable %in% "bin5p"), # We have to take just one variable (bin5p or bin3p) otherwise we get a sum for the polyA columns messing up the plot
      aes(x = factor(qname, levels = sort.index$qname), y = polya_length)) +
      geom_col(width=1) +
      coord_flip() +
      ylab("PolyA length") +
      theme_classic() + 
      theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank() 
        )
 
# If we want separately 5p and 3p plots or have more than one transcript per page
    # p <- ggarrange(main.plot.5p + rremove("y.text"), 
    #             main.plot.3p, 
    #             polya.bar, 
    #             ncol=3, nrow=1,
    #             align = "h",
    #             widths = c(3, 3, 1),
    #             common.legend = T) +
    #             ggtitle(tx_name)

# If we want both 5p and 3p plots together and we have one transcript per plot
    p <- ggarrange(main.plot + rremove("y.text"), 
                   polya.bar, 
                   ncol=2, nrow=1,
                   align = "h",
                   widths = c(3, 1),
                   common.legend = T)
    p<-annotate_figure(p, fig.lab = tx_name) # https://github.com/kassambara/ggpubr/issues/78
    print(p)
  }
dev.off()

graphics.off()