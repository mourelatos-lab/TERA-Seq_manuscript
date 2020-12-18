#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gtable))
suppressMessages(library(ggpubr))

# zscore calculates z-score for numeric vector x. Each column in x will have 0
# mean and standard deviation 1 after scaling.
zscore <- function(x) {
	z = scale(x, center = TRUE, scale = TRUE)
	return(z)
}

# Read command line options and arguments
option_list <- list(
	make_option(
		c("-i", "--ifile"), type = "character",
		help = "Input table", metavar = "File"),
	make_option(
		c("-f", "--figfile"), type = "character",
		help = "Output pdf file with graphs"),
	make_option(
		c("-u", "--posMin"), type = "integer",
		help = "The minimum position to be kept in the data", metavar = "Num"),
	make_option(
		c("-d", "--posMax"), type = "integer",
		help = "The maximum position to be kept in the data", metavar = "Num"),
	make_option(
		c("-s", "--subset"), type = "double",
		help = "Percent of top expressed transcripts to keep", metavar = "Num"),
	make_option(
	    c("-n", "--name"), type = "character",
	    help = "Plot title", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

name = opt$name
posMin = opt$posMin # for data filtering
posMax = opt$posMax # for data filtering
subsetPct = opt$subset

# Create data.table
dt = data.table(read.delim(opt$ifile))

# Filter by pos
dt = dt[pos >= posMin & pos <= posMax]

# Aggregate count at the transcript and position level (sum samples).
dt = dt[, .(count = sum(count)), by = .(ref, pos)]

# Aggregate count at the transcript level.
dt.tr = dt[, .(count = sum(count)), by = .(ref)]

# Normalize counts for each transcript by calculating z-score values.
dt[, zscore := zscore(count), by = .(ref)]

# Find a limit for the counts in order to filter lowly expressed transcripts.
countLim = quantile(dt.tr[count > 0]$count, probs = 1 - subsetPct)[1]

# How many transcripts above limit?
trCount = dim(dt.tr[count > countLim])[1]

# Keep only highly expressed transcripts.
expressed_transcripts = dt.tr[count > countLim]$ref
dt.tr = dt.tr[ref %in% expressed_transcripts]
dt = dt[ref %in% expressed_transcripts]

# Sort transcripts by count.
sorted_transcripts = dt.tr$ref[order(dt.tr$count)]

# Reorder transcript levels by expression.
dt$ref <- factor(dt$ref, levels = sorted_transcripts)
dt.tr$ref <- factor(dt.tr$ref, levels = sorted_transcripts)

# Reorder transcript levels by clustering.
dt.dcast = dcast(dt, ref ~ pos, value.var = "zscore")
hc = hclust( dist(dt.dcast[, !"ref"], method = "euclidean"), method = "ward.D" )
dt$ref <- factor(dt$ref, levels = dt.dcast[hc$order]$ref)
dt.tr$ref <- factor(dt.tr$ref, levels = dt.dcast[hc$order]$ref)

# Reorder transcript levels by position with maximum zscore.
dt.whichMaxZscore = dt[dt[, .I[which.max(zscore)], by = ref]$V1]
dt$ref <- factor(dt$ref, levels = dt.whichMaxZscore[order(pos)]$ref)
dt.tr$ref <- factor(dt.tr$ref, levels = dt.whichMaxZscore[order(pos)]$ref)

# Aggregate count and zscore at the position level and calculate mean, sd and se.
dt.pos = dt[, .(count_sum = sum(count),
                count_mean = mean(count),
                count_sd = sd(count),
                count_se = sd(count)/length(count),
                zscore_sum = sum(zscore),
                zscore_mean = mean(zscore),
                zscore_sd = sd(zscore),
                zscore_se = sd(zscore)/length(zscore))
            , by = .(pos)]

# Calculate density.
dt.pos[, count_density := count_sum / sum(count_sum)]

# Calculate the non-zero position with the highest count and zscore.
maxZScorePos = dt.pos$pos[which.max(dt.pos[pos != 0]$zscore_mean)]
maxCountPos = dt.pos$pos[which.max(dt.pos[pos != 0]$count_density)]
if (maxZScorePos > 0) {
    maxZScorePos = maxZScorePos + 1 # accounts for excluded index as pos 0
}
if (maxCountPos > 0) {
    maxCountPos = maxCountPos + 1 # accounts for excluded index as pos 0
}

# Create plots
pl.zscore_line <- ggplot(data = dt.pos, aes(x = pos, y = zscore_mean, ymin = zscore_mean - zscore_sd, ymax = zscore_mean + zscore_sd)) +
	geom_vline(xintercept = 0, colour = "darkgrey", linetype = "dotted", alpha = 1) +
	geom_vline(xintercept = maxZScorePos, colour = "orange", linetype = "dotted", alpha = 0.8) +
    geom_text(aes(x = maxZScorePos, y = max(dt.pos$zscore_mean)), label = maxZScorePos, hjust = -.1,  colour = "orange", size = 2) +
	geom_line() +
	ylab("mean z-score") +
	theme_classic() +
    theme(
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

pl.zscore_heatmap <- ggplot(dt, aes(x = pos, y = ref )) +
    geom_tile(aes(fill = zscore)) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    theme_classic() +
	theme(
	    legend.key.height = unit(5, "pt"),
	    axis.text.y = element_text(size = 1),
	    axis.ticks.y = element_blank())

pl.count_col <- ggplot(dt.tr, aes(x = ref, y = count)) +
    geom_col() +
    ylab("count") +
    coord_flip() +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

pdf(NULL) # Prevent ggarrange creating an empty pdf for unknown reason
p1 <- ggarrange(pl.zscore_line, NULL, pl.zscore_heatmap, pl.count_col,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(4, 1), heights = c(1, 4), common.legend = TRUE, legend = "bottom")

pl.count_line <- ggplot(data = dt.pos, aes(x = pos, y = count_density)) +
	geom_vline(xintercept = 0, colour = "darkgrey", linetype = "dotted", alpha = 1) +
	geom_vline(xintercept = maxCountPos, colour = "orange", linetype = "dotted", alpha = 0.8) +
    geom_text(aes(x = maxCountPos, y = max(dt.pos$count_density)), label = maxCountPos, hjust = -.1,  colour = "orange", size = 2) +
	geom_line() +
	ylab("count density") +
	theme_classic() +
    theme(
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

pl.count_heatmap <- ggplot(dt, aes(x = pos, y = ref )) +
    geom_tile(aes(fill = count)) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    theme_classic() +
	theme(
	    legend.key.height = unit(5, "pt"),
	    axis.text.y = element_text(size = 1),
	    axis.ticks.y = element_blank())

p2 <- ggarrange(pl.count_line, NULL, pl.count_heatmap, pl.count_col,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(4, 1), heights = c(1, 5), common.legend = TRUE, legend = "bottom")

pdf(opt$figfile, height = 14) # pdf(file = NULL) # Hack to prevent empty first page; see https://github.com/wilkelab/cowplot/issues/24
    annotate_figure(p1, top = sprintf("%s: transcrs: %d (%.1f%%)", name, trCount, subsetPct*100))
    annotate_figure(p2, top = sprintf("%s: transcrs: %d (%.1f%%)", name, trCount, subsetPct*100))
dev.off()
graphics.off()
