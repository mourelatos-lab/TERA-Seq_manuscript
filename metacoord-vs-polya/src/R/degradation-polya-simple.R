#!/usr/bin/env Rscript
#
# Get dependence of poly(A) tail length on degradation status
#

args = commandArgs(trailingOnly=TRUE)

opt<-NULL
opt$ifile<-args[1]

tab<-read.table(opt$ifile, sep="\t", stringsAsFactors = F, header=T)

tab$degrad<-tab$bin5p # +(19-tab$bin3p) # get overall degradation

tab<-tab[tab$qc_tag %in% c("PASS", "NOREGION"),]

# Testing
tab.bckp<-tab

#http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
#suppressPackageStartupMessages(library("dplyr"))

# Note: The position of the formula stat_cor & stat_regline is hardcoded to fit to 5P Poly(A) sample.
#       If you want to place the formula for other libraries you have to adjust it (or make it relative). 
print("Read-based results")

p <- ggplot(tab, aes(degrad, polya_length)) +
        geom_violin(aes(fill = factor(degrad))) +
        geom_jitter(size=0.075, alpha=0.5, color="black") +
        geom_smooth(method = "lm", color="red", se=TRUE, level = 0.95) +
        scale_x_continuous(breaks = 0:19) +
        scale_y_continuous(breaks = seq(50, max(tab$polya_length), by=50)) +
        coord_cartesian(ylim = c(0, 600)) +
        theme_classic() +
        theme(legend.position = "none") +
        stat_cor(method = "kendall", label.y.npc = 0.7, label.x.npc = 0.1) +
        stat_regline_equation( label.y.npc = 0.68, label.x.npc = 0.1)

p <- p  + xlab("Bin") + ylab("Poly(A) length (nt)") # xlab("5' end degradation level (bin)")

pdf(file = paste(args[2], paste0(sub(".[^.]+$", "", basename(opt$ifile)), "-polya-vs-degrad.reads.pdf"), sep="/"), width = 8)
        print(p)
dev.off()
print("Transcripts-based results")

suppressPackageStartupMessages(library("dplyr"))

tabsum <- tab %>%
  group_by(sample, feat, degrad) %>%
  summarize(polya_length = mean(polya_length)) %>%
  ungroup()

tab <- tabsum

p <- ggplot(tab, aes(degrad, polya_length)) +
  geom_violin(aes(fill = factor(degrad))) +
  geom_jitter(size=0.1, alpha=0.5, color="black") +
  geom_smooth(method = "lm", color="red", se=TRUE, level = 0.95) +
  scale_x_continuous(breaks = 0:19) +
  scale_y_continuous(breaks = seq(50, max(tab$polya_length), by=50)) +
  coord_cartesian(ylim = c(0, 600)) +
  theme_classic() +
  theme(legend.position = "none") +
  stat_cor(method = "kendall", label.y.npc = 0.9, label.x.npc = 0.1) +
  stat_regline_equation(label.y.npc = 0.88, label.x.npc = 0.1)

p <- p  + xlab("Bin") + ylab("Poly(A) length (nt)") # xlab("5' end degradation level (bin)")

pdf(file = paste(args[2], paste0(sub(".[^.]+$", "", basename(opt$ifile)), "-polya-vs-degrad.transcripts.pdf"), sep="/"), width = 8)

  print(p)
dev.off()