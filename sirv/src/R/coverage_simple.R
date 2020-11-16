#!/usr/bin/env Rscript
#
# Very simple coverage plot for example from samtools depth
# One mandatory argument (input tab) and one optional (output pdf)
#

args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
suppressPackageStartupMessages(library(dplyr))

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste0(gsub("\\..*", "", intab), ".pdf")
}

intab<-args[1]

tab<-read.table(intab, sep="\t", header=T, stringsAsFactors = F)

# This helps us to remove "empty" regions if the coverage is from samtools depth
# Using the ggplot2 afterwards would plot a long region with an even coverage which is caused by the "jump" in the coords
if(length(unique(tab$Library))>1){
  tab_final<-tab[0,]
  tab_final$coverage<-NULL
  tab_final$Library<-NULL
  for(lib in unique(tab$Library)){
    print(lib)
    tab1<-tab %>% filter(Library==!!lib)
    colnames(tab1)[colnames(tab1)=="coverage"]<-lib
    tab1$Library<-NULL
    tab_final<-tab_final %>% full_join(tab1, by = c("chr", "coord"))
  }
}else{
  tab_final<-tab
}

tab_final$coord<-seq(1, nrow(tab_final), by=1)

if(length(unique(tab$Library))>1){
  tab_final<-reshape2::melt(tab_final, variable.name="Library", measure.vars=colnames(tab_final)[!(colnames(tab_final) %in% c("chr", "coord"))], value.name="coverage")
}

tab_final$Library<-factor(tab_final$Library, levels = unique(sort(tab_final$Library)))

tab_final<-tab_final %>% group_by(Library) %>%
  mutate(coverage_rel=coverage/sum(coverage))

len <- tab_final %>% group_by(Library) %>% count()
len <- max(len$n)
name<-paste0("SIRV",
             strsplit(
               strsplit(x = basename(intab), split = "SIRV")[[1]][2],
               ".", fixed = T)[[1]][1],
             ": ", len, "nt"
)


p<-ggplot(tab_final, aes(coord, coverage_rel, color=Library)) +
  geom_line(aes(color = Library), stat = "summary_bin", binwidth = 10) +
  theme_bw()+
  theme(strip.text.y = element_text(size = 8, colour = "black", angle = 0),
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  ggtitle(name)

pdf(file = args[2], width = 8, height = 6)
  print(p)
dev.off()
