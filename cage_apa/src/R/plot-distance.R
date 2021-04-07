#!/usr/bin/env Rscript
#
# Plot overlap bed from bedtools closest with added first column as sample name(s)
# Assumes columns: c("dataset", "chr", "start", "end", "id", "score", "strand", "database", "chr_d", "start_d", "end_d", "id_d", "score_d", "strand_d", "distance")
#

library(optparse)
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Input table to subset from. For streaming type \"stdin\"", metavar = "File"),
  make_option(
    c("-d", "--distance"), type = "integer", default = 10,
    help = "Maximum distance of the \"closest\" to consider"),
  make_option(
    c("-s", "--direction"), type = "character",
    help = "Consider upstream or downstream positions? [up|down|both]. Note: based on genomic distances (upstream on + strand is lower number, upstream on - strand is higher number, etc.)"),
  make_option(
    c("-c", "--collapse"), action = "store_true", default = FALSE,
    help = "Collapse the positions (count one overlap position only once)"),
  make_option(
    c("-z", "--zoom"), action = "store_true", default = FALSE,
    help = "Zoom to --distance from the total distribution for the plot. Default: Consider only the --distance for the plot."),
  make_option(
    c("-k", "--keepChr"), action = "store_true", default = FALSE,
    help = "Keep all the chromosomes, not just the \"nice ones\""),
  make_option(
    c("-o", "--ofile"), type = "character",
    help = "Output subset table [pdf|ps]", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

### Testing variables
  # print("USING TEST VARIABLES!!!")
  # opt<-NULL
  # opt$ifile<-"/home/joppelt/projects/ribothrypsis/analysis/dist-from-annot/results/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/reads.1.sanitize.toGenome.sorted.3p-apa-upstream.wAdapter.bed" # Table to subset from
  # opt$ifile<-"reads.1.sanitize.toGenome.sorted.5p-cage-upstream.wAdapter.bed" # Table to subset from
  # opt$distance<-10
  # opt$keepChr<-FALSE
  # opt$direction<-"up"
  # opt$collapse<-FALSE
  # opt$zoom<-FALSE
  # opt$ofile<-"/home/joppelt/projects/ribothrypsis/analysis/dist-from-annot/results/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/reads.1.sanitize.toGenome.sorted.3p-apa.TEST.pdf" # Table to subset by
  # opt$ofile<-"test.pdf"
  # print("USING TEST VARIABLES!!!")
### Testing variables

if(opt$ifile=="stdin"){
  input_bed<-read.table(file("stdin"), sep="\t", header=F, stringsAsFactors = T)
}else{
  input_bed<-read.table(opt$ifile, sep="\t", header=F, stringsAsFactors = T)
}
colnames(input_bed)<-c("dataset", "chr", "start", "end", "id", "score", "strand", "database", "chr_d", "start_d", "end_d", "id_d", "score_d", "strand_d", "distance")

if(opt$keepChr==FALSE){
  print("Keeping only the nice chromosomes.")
  input_bed<-input_bed %>%
    filter(chr %in% c(seq(1,22,by=1),"X","Y","MT"))
  input_bed$chr<-droplevels(input_bed$chr) # remove unused chromosome levels
}

if (opt$collapse == TRUE) {
  print("Collapsing positions.")
  input_bed$dupl <- paste(input_bed$dataset, input_bed$chr, input_bed$start, input_bed$end, input_bed$strand, sep = ".")
  input_bed <- input_bed[!duplicated(input_bed$dupl), ]
  input_bed$dupl <- NULL
} else {
  print("Keeping all positions.")
}

plot_bed<-NULL # Object for plotting
for(dataset.u in unique(input_bed$dataset)){
  print(dataset.u)
  tmp<-input_bed %>%
    filter(dataset %in% dataset.u)

  nreads<-tmp %>%
    select(id) %>%
    unique() %>%
    nrow()

  # Replace "-1" for reads with no overlap with 999999; bedtools stores "-1" in distance for no overlap column which could be a problem
  tmp<-tmp %>%
    mutate(distance=replace(distance, database==".", 999999))

  if(opt$direction == "up"){
    print("Counting only upstream or overlapping features.")
    xlims<-seq(-opt$distance, 0, by=1)
    tmp.f<-tmp %>%
      filter(database != ".") %>%
      filter(distance %in% xlims)
  }else if(opt$direction == "down"){
    print("Counting only downstream or overlapping features.")
    xlims<-seq(0, opt$distance, by=1)
    tmp.f<-tmp %>%
      filter(database != ".") %>%
      filter(distance %in% xlims)
  }else if(opt$direction == "both"){
    print("Counting both upstream and downstream or overlapping features.")
    xlims<-seq(-opt$distance, opt$distance, by=1)
    tmp.f<-tmp %>%
      filter(database != ".") %>%
      filter(distance %in% xlims)
  }else{
    stop("Please specify correct direction - up|down|both!")
  }

  nreads.f<-tmp.f %>%
    select(id) %>%
    unique() %>%
    nrow()

  print(paste("Stats for dataset:", dataset.u))
  print(paste("Total number of positions is:", nreads))
  print(paste("Distance for overlap is:",  paste(xlims, collapse = ", ")))
  print(paste("Number of overlapping positions is:", nreads.f))
  print(paste("Perc. of overlap is:", round(nreads.f/(nreads/100),2)))

  tmp<-as.data.frame(table(tmp$distance), stringsAsFactors = F)
  tmp$perc<-tmp$Freq/(sum(tmp$Freq)/100)
  tmp$dataset<-dataset.u

  if(is.null(plot_bed)){
    plot_bed<-tmp
  }else{
    plot_bed <- rbind(plot_bed, tmp)
  }
}

plot_bed$Var1<-as.numeric(as.character(plot_bed$Var1))
plot_bed$dataset<-factor(plot_bed$dataset, levels = unique(plot_bed$dataset))

# Add cummulative sum only for the selected span
if (opt$direction == "up") {
  plot_bed_cum <- plot_bed %>%
    group_by(dataset) %>%
    filter(if (!!opt$zoom == FALSE) Var1 <= 0 & Var1 >= (-opt$distance) else Var1 <= 0) %>%
    mutate(cumsum = cumsum(perc)) %>%
    ungroup()
} else if (opt$direction == "down") {
  plot_bed_cum <- plot_bed %>%
    group_by(dataset) %>%
    filter(if (!!opt$zoom == FALSE) Var1 >= 0 & Var1 <= opt$distance else Var1 >= 0) %>%
    mutate(cumsum = cumsum(perc)) %>%
    ungroup()
} else if (opt$direction == "both") {
  plot_bed_cum <- plot_bed %>%
    group_by(dataset) %>%
    filter(if (!!opt$zoom == FALSE) Var1 <= opt$distance & Var1 >= (-opt$distance)) %>%
    mutate(cumsum = cumsum(perc)) %>%
    ungroup()
} else {
  stop("Don't know the direction. Please select one of the [\"up\", \"down\", \"both\"]")
}

# Uncomment this if you want to extend the plotting limits (not the stats limits, just the plotting limits)
# if(opt$direction == "up"){
#   xlims_plot<-seq(-(ceiling(opt$distance/2))-opt$distance, 0, by=1)
# }else if(opt$direction == "down"){
#   xlims_plot<-seq(0, opt$distance+(ceiling(opt$distance/2)), by=1)
# }else if(opt$direction == "both"){
#   xlims_plot<-seq(-(ceiling(opt$distance/2))-opt$distance, opt$distance+(ceiling(opt$distance/2)), by=1)
# }else{
#   xlims_plot<-xlims
# }
xlims_plot<-xlims
if(length(xlims_plot)>=45 & length(xlims_plot)<=100){
  text_size<-375/length(xlims_plot)
}else{
  text_size<-10
}

# Auto check for output file (pdf or ps/eps)
# Note: eps cannot do alpha
output_type<-sub('.*\\.', '', opt$ofile)
if(output_type=="pdf"){
  pdf(opt$ofile, height = 6, width = 8)
  alpha_plot<-0.6
}else if(output_type=="ps" | output_type=="eps"){
  postscript(opt$ofile, height = 6, width = 8, paper="special")
  alpha_plot<-1.0
}else{
  pdf(opt$ofile, height = 6, width = 8)
  alpha_plot<-0.6
}
  # "raw"
  p <- ggplot(plot_bed, aes(x = Var1, y = perc, color = dataset, shape = dataset)) +
    geom_line(size = 0.8, alpha = alpha_plot) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    ggtitle("Distance to the annotated feature(s)") +
    xlab("Distance (bp)") +
    ylab("Percentage of reads (%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 8)) +
    theme(
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
      axis.text.x = element_text(size = text_size)
    )

  # cummulative
  p1 <- ggplot(plot_bed_cum, aes(x = Var1, y = cumsum, color = dataset, shape = dataset)) +
    geom_line(size = 0.8, alpha = alpha_plot) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    ggtitle("Distance to the annotated feature(s)") +
    xlab("Distance (bp)") +
    ylab("Cummulative percentage of reads (%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 8)) +
    theme(
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
      axis.text.x = element_text(size = text_size)
    )

  if(opt$zoom == TRUE){
    p <- p +
      coord_cartesian(ylim = c(0, 100),
                      xlim = c(min(xlims_plot), max(xlims_plot)))
    p1 <- p1 +
      coord_cartesian(ylim = c(0, 100),
                      xlim = c(min(xlims_plot), max(xlims_plot)))
  }else{
    p <- p +
      coord_cartesian(xlim = c(min(xlims_plot), max(xlims_plot)))
    p1 <- p1 +
      coord_cartesian(xlim = c(min(xlims_plot), max(xlims_plot)))
  }

	if(length(xlims_plot)<=100){
		p <- p +
				geom_point(size=1.5, alpha=alpha_plot) +
				scale_x_continuous(breaks = xlims_plot, labels = xlims_plot, limits = c(min(xlims_plot), max(xlims_plot)))
		p1 <- p1 +
				geom_point(size=1.5, alpha=alpha_plot) +
				scale_x_continuous(breaks = xlims_plot, labels = xlims_plot, limits = c(min(xlims_plot), max(xlims_plot)))
 	}else{
 	  scale_x <- c(min(xlims_plot), round(min(xlims_plot)/4, 0)*3, round(min(xlims_plot)/2, 0), round(min(xlims_plot)/4, 0)*1,
 	               0,
 	               round(max(xlims_plot)/4, 0)*1, round(max(xlims_plot)/2, 0), round(max(xlims_plot)/4, 0)*3, max(xlims_plot))

		p <- p +
				scale_x_continuous(breaks = scale_x,
									         labels = scale_x, limits = c(min(xlims_plot), max(xlims_plot)))
		p1 <- p1 +
				scale_x_continuous(breaks = scale_x,
								           labels = scale_x, limits = c(min(xlims_plot), max(xlims_plot)))
	}

  if (opt$direction == "both") {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.5)
    p1 <- p1 + geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.5)
  }

	print(p)
	print(p1)
dev.off()
