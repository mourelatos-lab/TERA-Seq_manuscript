#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

corr_square_eqn <- function(x, y, digits = 2) {
	corr_coef <- round(cor(x, y, use="pairwise.complete.obs"), digits = digits)
	paste("italic(r^2) == ", corr_coef^2)
}

corr_eqn <- function(x, y, digits = 2) {
	corr_coef <- round(cor(x, y, use="pairwise.complete.obs"), digits = digits)
	paste("italic(r) == ", corr_coef)
}

positive_non_zero_min <- function(x) {
    return(min(x[x > 0]))
}

# Read command line options and arguments
option_list <- list(
	make_option(
		c("-i", "--ifile"), type="character",
		help="Input table file"),
	make_option(
		c("-f", "--figfile"), type="character",
		help="Output pdf file with graphs")
)
opt = parse_args(OptionParser(option_list = option_list))

### TESTING VARIABLES ###
#print("USING TESTING VARIABLES!!!")
#opt<-NULL
#opt$ifile<-"counts.corrected-mrna.unambig.rel5.tab"
#opt$figfile<-"compare-counts.pdf"
#print("USING TESTING VARIABLES!!!")
### TESTING VARIABLES ###

# Read data
dt = data.table(read.delim(opt$ifile))

# Calculate RPM (appropriate for long reads).
dt[, rpm_total := (total / sum(total)) * 10^6, by = .(sample)]
dt[, rpm_full := (full / sum(total)) * 10^6, by = .(sample)]
dt[, rpm_fullcds := (fullcds / sum(total)) * 10^6, by = .(sample)]

# Calculate RPKM (appropriate for short reads).
dt[, rpkm_total := ((total / len) / sum(total)) * 10^9, by = .(sample)]
dt[, rpkm_full := ((full / len) / sum(full)) * 10^9, by = .(sample)]
dt[, rpkm_fullcds := ((fullcds / cds_len) / sum(fullcds)) * 10^9, by = .(sample)]

# Calculate percent of full length and full length CDS.
dt[, percent_full := (full / total)*100]
dt[, percent_fullcds := (fullcds / total)*100]

# Aggregate all columns
dt = dt[, .(
  total = mean(total),
  full = mean(full),
  fullcds = mean(fullcds),
  len = mean(len),
  cds_len = mean(cds_len),
  rpm_total = mean(rpm_total),
  rpm_full = mean(rpm_full),
  rpm_fullcds = mean(rpm_fullcds),
  rpkm_total = mean(rpkm_total),
  rpkm_full = mean(rpkm_full),
  rpkm_fullcds = mean(rpkm_fullcds),
  percent_full = mean(percent_full),
  percent_fullcds = mean(percent_fullcds)
), by = .(ref)]

# Sort table by total read count.
dt <- dt[order(-total),]

# Write output file
write.table(file = gsub(".pdf", ".tsv", opt$figfile), dt, quote = FALSE, sep = "\t", row.names = FALSE)

# Discard all zeros.
dt = dt[rpm_total > 0]

# Add a small positive non-zero value to avoid problems with log transformations.
dt$rpm_total = dt$rpm_total + positive_non_zero_min(dt$rpm_total)
dt$rpm_full = dt$rpm_full + positive_non_zero_min(dt$rpm_full)
dt$rpm_fullcds = dt$rpm_fullcds + positive_non_zero_min(dt$rpm_fullcds)
dt$rpkm_total = dt$rpkm_total + positive_non_zero_min(dt$rpkm_total)
dt$rpkm_full = dt$rpkm_full + positive_non_zero_min(dt$rpkm_full)
dt$rpkm_fullcds = dt$rpkm_fullcds + positive_non_zero_min(dt$rpkm_fullcds)

# Find a limit for the expression in order to filter lowly expressed transcripts.
rpkmLim = quantile(dt[rpkm_total > 0]$rpkm_total, probs = .25)[1]
rpmLim = quantile(dt[rpm_total > 0]$rpm_total, probs = .25)[1]

# Keep only highly expressed transcripts.
expressed_transcripts = dt[rpm_total > rpmLim]$ref
dt = dt[ref %in% expressed_transcripts]

# Create plots.
pdf(opt$figfile)
  ggplot(dt, aes(x = rpm_total, y = rpm_fullcds, color = rpm_full)) +
  	geom_point(size = 1.5, alpha = .6, shape = 16) +
  	stat_cor() +
  	scale_x_log10() +
  	scale_y_log10() +
  	theme_classic() +
    scale_color_continuous(high = "#132B43", low = "#56B1F7", name="full\nmRNAPM") +
    xlab("log10(RPM)") + ylab("log10(full CDSPM)")
  ggplot(dt, aes(x = rpm_total, y = rpm_full, color = rpm_fullcds)) +
  	geom_point(size = 1.5, alpha = .6, shape = 16) +
  	stat_cor() +
  	scale_x_log10() +
  	scale_y_log10() +
  	theme_classic() +
    scale_color_continuous(high = "#132B43", low = "#56B1F7", name="full\nCDSPM") +
    xlab("log10(RPM)") + ylab("log10(full mRNAPM)")
  ggplot(dt, aes(x = rpm_total, y = percent_full, color = percent_fullcds)) +
  	geom_point(size = 1.5, alpha = .6, shape = 16) +
  	stat_cor() +
  	scale_x_log10() +
  	theme_classic() +
    scale_color_continuous(high = "#132B43", low = "#56B1F7", name="Percent\nfull CDSPM") +
    xlab("log10(RPM)") + ylab("Percentage of full mRNA")
  ggplot(dt, aes(x = percent_fullcds, y = percent_full, color = log2(rpm_total))) +
    geom_point(size = 1.5, alpha = .7, shape = 16) +
    stat_cor() +
    theme_classic() +
    xlab("Percentage of full CDS") + ylab("Percentage of full mRNA") + 
    scale_colour_viridis_b(option = "inferno", name="log2(RPM)", direction = -1) +
    theme(legend.direction = "horizontal", legend.position = "bottom")
  ggplot(dt, aes(x = percent_full)) +
    geom_density() +
    theme_classic() +
    xlab("Density") + ylab("Percentage of full mRNA")
  ggplot(dt, aes(x = log2(rpm_total), y = log2(percent_full))) +
    geom_point() +
    geom_smooth(method = "lm")
    theme_classic() +
    xlab("log2(RPM)") + ylab("log2(Percentage of full mRNA")
graphics.off()
