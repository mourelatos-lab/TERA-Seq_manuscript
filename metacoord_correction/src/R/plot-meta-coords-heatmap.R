#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
suppressMessages(library(data.table))
suppressMessages(library(grid))
library(gtable)
suppressMessages(library(ggpubr))

# Read command line options and arguments
option_list <- list(
    make_option(
        c("-i", "--ifile"), type = "character",
        help = "Input table", metavar = "File"),
	make_option(
	    c("-s", "--sub"), type = "integer", default = 50000,
	    help = "keep a random subset of reads for meta-coords scatterplot (prevents overplotting)", metavar = "Num"),
    make_option(
        c("-f", "--figfile"), type = "character",
        help = "Output pdf file with graphs")
)
opt = parse_args(OptionParser(option_list = option_list))

# Create data table from file.
dt = data.table(read.delim(opt$ifile))
if (nrow(dt) < opt$sub) {
    opt$sub = nrow(dt)
}

# Heatmap: Measure the number of reads for each [bin5p, bin3p] pair.
dt.heatmap = dt[, .(count = .N), by = .(bin5p, bin3p)]

# Create plots

scatter_geom = geom_point(
    data = dt[sample(nrow(dt), opt$sub), ],
    aes(x = jitter(bin5p, 1.6), y = jitter(bin3p, 1.6)), size = .05, alpha = .05)

pl.heatmap <- ggplot(dt.heatmap, aes(x = bin5p, y = bin3p )) +
    geom_tile(aes(fill = log10(count)), color = "white") +
    scale_fill_gradient(low = "white", high = "firebrick") +
    xlab("5' end position (bin)") +
    ylab("3' end position (bin)") +
    theme_classic()

pl.bin5p.bar <- ggplot(dt, aes(x = bin5p)) +
    geom_bar(aes(y = (100*..count..)/sum(..count..)),alpha = .5) +
    ylab("read percent") +
    theme_classic() +
    theme(
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

pl.bin3p.bar <- ggplot(dt, aes(x = bin3p)) +
    geom_bar(aes(y = (100*..count..)/sum(..count..)),alpha = .5) +
    ylab("read percent") +
    coord_flip() +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

pdf(file = NULL) # Hack to prevent empty first page; see https://github.com/wilkelab/cowplot/issues/24
p <- ggarrange(pl.bin5p.bar, NULL, pl.heatmap + scatter_geom, pl.bin3p.bar,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(4, 1), heights = c(1, 4),
          common.legend = TRUE)
graphics.off()

pdf(opt$figfile, width = 7)
p
graphics.off()
