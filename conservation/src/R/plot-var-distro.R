#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(ggplot2)

# Read command line options and arguments
option_list <- list(
	make_option(
		c("-i", "--ifile"), type="character",
		help="Input table", metavar="File"),
	make_option(
		c("-f", "--figfile"), type="character",
		help="Output pdf file with graphs"),
	make_option(
		c("-n", "--name"), type="character",
		help="Plot title"),
	make_option(
		c("-m", "--negate"), action="store_true", default=FALSE,
		help="Make variable values negative")
)
opt = parse_args(OptionParser(option_list = option_list))

# Read data
df = read.delim(opt$ifile)

# Crete data.table.
dt = data.table(df)

# Aggregate all columns by position (in case a sample column exists).
dt = dt[, .(val = sum(val), records = sum(records)), by = .(pos, cond)]

# Normalize val by the records.
dt[, norm := val / records]

# Negate values if option is set.
ylabel = ""
if (opt$negate){
	dt$norm = -dt$norm
	ylabel = "(negative)"
}

# Calculate z-scores
dt[, centered := scale(norm, center = TRUE, scale = FALSE), by = .(cond)]

# Create plots
pdf(opt$figfile)
ggplot(dt, aes(x = pos, y = norm, colour = cond)) +
	geom_vline(xintercept=0, colour="darkgrey", linetype="dotted", alpha=1) +
	geom_line() +
	xlab("position") +
	ylab(paste(sep=" ", "normalized value", ylabel)) +
	ggtitle(sprintf("%s", opt$name)) +
	theme_classic() +
    theme(legend.position="bottom")
ggplot(dt, aes(x = pos, y = centered, colour = cond)) +
	geom_vline(xintercept=0, colour="darkgrey", linetype="dotted", alpha=1) +
	geom_line() +
	xlab("position") +
	ylab(paste(sep=" ", "centered normalized value", ylabel)) +
	ggtitle(sprintf("%s", opt$name)) +
	theme_classic() +
    theme(legend.position="bottom")
ggplot(dt, aes(x = pos, y = centered, colour = cond)) +
	geom_vline(xintercept=0, colour="darkgrey", linetype="dotted", alpha=1) +
	geom_smooth(method = "loess", se = FALSE, span = .1) +
	xlab("position") +
	ylab(paste(sep=" ", "centered normalized value", ylabel)) +
	ggtitle(sprintf("%s", opt$name)) +
	theme_classic() +
    theme(legend.position="bottom")
graphics.off()
