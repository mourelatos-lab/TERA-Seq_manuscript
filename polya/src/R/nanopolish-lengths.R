#
# Comparing poly(A) tails with
#
# TODO: Remove hard-coded section and replace them with command-line inputs
#

library(rio)
library(dplyr)
library(ggplot2)
library(patchwork)
#library(reshape2) # used only once to melt

odir <- "results"
dir.create(odir, recursive = T)

min_count <- 10 # Min. count of reads to include the transcript

paths <- c(
   "../samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab",
   "../samples/hsa.dRNASeq.HeLa.polyA.REL5.long.1/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab",
   "../samples/hsa.dRNASeq.HeLa.total.REL3.3/align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.tab"
)

subsets <- c(
   "results/other/nanopolish-checks/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/coords-on-corrected-meta-mrnas.primary.19bin.names.txt",
   "results/other/nanopolish-checks/hsa.dRNASeq.HeLa.polyA.REL5.long.1/coords-on-corrected-meta-mrnas.primary.19bin.names.txt",
   "results/other/nanopolish-checks/hsa.dRNASeq.HeLa.total.REL3.3/coords-on-corrected-meta-mrnas.primary.19bin.names.txt"
)

samples <- c(
  "CAP-Poly(A)",
  "5P-Poly(A)",
  "TERA3"
)

groups <- c(
  rep("5TERA", 2),
  rep("TERA3", 1)
)

samples_table <- data.frame(
  polya_path = paths,
  sample_name = samples,
  group = groups
)

polya_all <- NULL
for (i in 1:length(paths)) {
  polya <- rio::import(paths[i], format = "tsv")
  reads <- rio::import(subsets[i], format="tsv", header=F) # subsample read to whatever we need, here only reads reaching the end of the molecule (3' end at bin 19)
  polya <- polya %>%
    filter(readname %in% reads$V1) %>%
    mutate(library = !!samples[i], group = groups[i]) %>%
    select(contig, polya_length, qc_tag, library, group)
  if (is.null(polya_all)) {
    polya_all <- polya
  } else {
    polya_all <- rbind(polya_all, polya)
  }
}

polya_all <- polya_all %>%
  filter(qc_tag == "PASS" | qc_tag == "NOREGION")

polya_all <- polya_all %>%
  group_by(contig, library) %>%
  mutate(reads = n()) %>%
  filter(reads >= !!min_count) %>%
  mutate(med = median(polya_length)) %>%
  ungroup()

# Make scatterplots
# TODO: Consider changing to ggpairs https://www.r-graph-gallery.com/199-correlation-matrix-with-ggally.html
pairs <- t(combn(unique(polya_all$library), 2)) # Make all pairs
pairs <- apply(pairs, 2, as.character) # remove factors

for (i in 1:nrow(pairs)) {
  pair_names <- pairs[i, ]
  # Check we don't compare the same samples
  if (sum(table(pair_names) == 2) > 0) {
    next
  }
  xval <- subset(polya_all, polya_all$library %in% pair_names[1], c("contig", "med"))
  xval <- unique(xval)
  colnames(xval)[ncol(xval)] <- "xval"
  yval <- subset(polya_all, polya_all$library %in% pair_names[2], c("contig", "med"))
  yval <- unique(yval)
  colnames(yval)[ncol(yval)] <- "yval"
  pair <- full_join(x = xval, y = yval)

  #  pair[is.na(pair)] <- 0
  #  pair <- pair[rowSums(pair[-1])>0, ] # Remove not expressed genes in neither of samples
  pair <- pair[rowSums(is.na(pair)) == 0, ] # Remove transcripts which are not in both libraries

  sink(paste(odir, "nanopolish-lengths_stats.txt", sep = "/"), append = T)
    print("Number of transcripts in common:")
    print(pair_names)
    print(nrow(pair))
  sink()

  #corr_plot <- cor(pair$xval, pair$yval, method = "spearman")
  #corr_plot <- cor(pair$xval, pair$yval, method = "kendall")
  corr_plot <- cor(pair$xval, pair$yval, method = "pearson") # data must have normal distribution hist(pair$xval); hist(pair$yval)

  p <- ggplot(pair, aes(x = xval, y = yval)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(intercept = 0, slope = 1, col = "darkgrey", linetype = "dashed", size = 1.5) +
#    geom_smooth(method='lm', se=TRUE, fullrange=TRUE, formula= y~x, col = "blue", size = 1.5, alpha = 0.75) +
    coord_cartesian(xlim = c(0, round(max(polya_all$med)/50)*50), ylim = c(0, round(max(polya_all$med)/50)*50)) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 22)
    ) +
#    ggtitle("Median poly(A) tail length scatterplot") +
    xlab(paste(pair_names[1], "Median poly(A)")) +
    ylab(paste(pair_names[2], "Median poly(A)"))

  p <- p +
    annotate("text", x = 40, y = 300, label = paste0("Pearson: ", round(corr_plot, 3)), size = 10)

  pair.melt<-reshape2::melt(pair)
  p2 <- ggplot(data = pair.melt, aes(x=variable, y=value, fill=variable)) +
    geom_boxplot() +
    theme_bw() +
    theme(
      axis.title.x=element_blank(),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 17),
      legend.position="bottom",
      legend.text = element_text(size = 16)
      )  +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    scale_fill_discrete(name = "", labels = c(pair_names[1], pair_names[2])) +
#    xlab("") +
    ylab("Median poly(A)")

  p <- p +
    inset_element(p2, left = 0.7, bottom = 0.01, right = 0.99, top = 0.4)

  # Heatmap version https://stackoverflow.com/questions/7714677/scatterplot-with-too-many-points; alt. https://slowkow.com/notes/ggplot2-color-by-density/
  # ggplot(pair, aes(xval, yval)) +
  #   stat_density_2d(aes(fill = stat(density)), geom = 'raster', contour = FALSE) +
  #   scale_fill_viridis_c() +
  #   coord_cartesian(expand = FALSE) +
  #   geom_point(col = 'white', size = 0.5)  + # shape = '.',
  #   geom_abline(intercept = 0, slope = 1, col = "white", linetype = "dashed", size = 1.5) +
  #   theme_bw() +
  #   theme(
  #     axis.text = element_text(size = 18),
  #     axis.title = element_text(size = 22)
  #   ) +
  #   xlab(paste(pair_names[1], "Median poly(A)")) +
  #   ylab(paste(pair_names[2], "Median poly(A)"))

  png(paste(odir, paste0(pair_names[1], "-vs-", pair_names[2], ".png"), sep = "/"), width = 1000, height = 1000)
    print(p)
  dev.off()
}

system(paste0("convert ", odir, "/*.png ", odir, "/nanopolish-lengths.pdf"), intern = TRUE)
