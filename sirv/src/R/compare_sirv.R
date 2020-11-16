#!/usr/bin/env Rscript
#
# Compare levels of the sirv spikeins
# The only input argument is a pointer to SIRV annotation in tsv format with molarities
#
##################
# PLEASE, make sure you are using the correct SIRV mix (here E2, but could be E0 or E1)
##################
#

args = commandArgs(trailingOnly=TRUE)

library("ggplot2")
library("ggrepel")

ifile<-list.files(pattern="reads.1.sanitize.noribo.toTranscriptome-.*.sorted.idxstats")
idxstats<-read.table(ifile, sep="\t", header=T)[,1:4]
annot<-read.table(args[1], sep="\t", header=T)

### TESTING VARIABLES!!!
#print("USING TEST VARIABLES!!!")
#idxstats<-read.table("/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/spikein/results/hsa.dRNASeq.SIRV.polyA.REL5.long.1/reassign/reads.1.sanitize.noribo.toTranscriptome.sorted.idxstats", sep="\t", header=T)[,1:4]
#print("USING TEST VARIABLES!!!")
### TESTING VARIABLES!!!

counts<-merge(annot, idxstats, by="transcript_id", all = T)
counts<-counts[!(is.na(counts$E2)),] # Remove transcripts not in E2
counts<-counts[order(counts$E2),]
counts$E2<-as.character(as.numeric(counts$E2))

counts<-counts[!is.na(counts$library),]
counts$count<-counts$count+1

# Do a simple normalization, just based on the input sum
for(libr in unique(counts$library)){
  counts[counts$library==libr, "count_norm"]<-counts[counts$library==libr, "count"]/(sum(counts[counts$library==libr, "count"])/1000000)
}

colnames(counts)[colnames(counts)=="transcript_id"]<-"Transcript"
colnames(counts)[colnames(counts)=="library"]<-"Library"

# Add length bins
counts<-counts[!(is.na(counts$length)),] # Remove the extra 0 column from the idxstats
max.val<-ceiling(max(counts$length)/500)*500 # round to closest 500 since we'll do 500 bins
counts$length_bin<-cut(counts$length, breaks = seq(0, max.val, by=500),
                       labels = paste(seq(0, max.val-500, by=500), seq(500, max.val, by=500), sep="-"), include.lowest = T)

#http://molbiol.edu.ru/eng/scripts/01_07.html
if(length(grep("cDNASeq", counts$Library))>1 && length(grep("dRNASeq", counts$Library))>1){
  counts.in<-counts[grep("dRNASeq", counts$Library),]
  pdf("molarity_vs_count_dRNASeq.pdf", width=8)
    for(molar in unique(counts.in$E2)){
      print(molar)
       p<-ggplot(data = subset(counts.in, E2 %in% molar), aes(x=length_bin, y=count_norm)) +
        geom_boxplot() +
        geom_jitter(aes(shape=Library, color=Transcript), position=position_jitter(0.2), size=3.5) +
        coord_trans(y = "log2") +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
        ggtitle(paste("Molarity:", molar)) +
        xlab("Length") + ylab("Normalized counts (log2 RPM)") +
        guides(col = guide_legend(ncol = 4), shape = guide_legend(ncol = 1)) +
        guides(col = F) +
        theme_classic() +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

      p <- p + geom_hline(yintercept = median(subset(counts.in, E2 %in% molar)$count_norm), color="darkgrey")

      print(p)
    }
  dev.off()

  counts.in<-counts[grep("cDNASeq", counts$Library),]
  pdf("molarity_vs_count_cDNASeq.pdf", width=8)
    for(molar in unique(counts$E2)){
      p<-ggplot(data = subset(counts.in, E2 %in% molar), aes(x=length_bin, y=count_norm)) +
        geom_boxplot() +
        geom_jitter(aes(shape=Library, color=Transcript), position=position_jitter(0.2), size=3.5) +
        coord_trans(y = "log2") +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
        ggtitle(paste("Molarity:", molar)) +
        xlab("Length") + ylab("Normalized counts (log2 RPM)") +
        guides(col = guide_legend(ncol = 4), shape = guide_legend(ncol = 1)) +
        guides(col = F) +
        theme_classic() +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

      p <- p + geom_hline(yintercept = median(subset(counts.in, E2 %in% molar)$count_norm), color="darkgrey")

      print(p)
    }
  dev.off()

}

counts$Library_type<-counts$Library
counts$Library_type<-as.character(counts$Library_type)

if(length(grep("cDNASeq", counts$Library))>1){
  counts[grep("cDNASeq", counts$Library), "Library_type"]<-"cDNASeq"
}
if(length(grep("dRNASeq", counts$Library))>1){
  counts[grep("dRNASeq", counts$Library), "Library_type"]<-"dRNASeq"
}

counts$Library_type<-factor(counts$Library_type, levels=unique(counts$Library_type))

# Error "In trans$transform(out$range) : NaNs produced" for out calculation here is not an error of the data but
#   a ggplot lo2 error https://rstudio-pubs-static.s3.amazonaws.com/59196_30f45c3acd6d4036a9506f5d83170d44.html
# If works well on ggplot2_3.2.0 & ggrepel_0.8.1, R 3.6.3 but doesn't on ggplot2_3.1.1 & ggrepel_0.8.2, R 3.5.1
pdf("molarity_vs_count.pdf", width=8)
  for(molar in unique(counts$E2)){
    tryCatch({
       p<-ggplot(data = subset(subset(counts, E2 %in% molar)), aes(x=length_bin, y=count_norm)) +
        coord_trans(y = "log2") +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
         geom_boxplot() +
         geom_jitter(aes(shape=Library, color=Transcript), position=position_jitter(0.2), size=3.5) +
        scale_shape_manual(values=c(15,16,17, 4:nlevels(subset(subset(counts, E2 %in% molar))$Library))) +
        ggtitle(paste("Molarity:", molar)) +
        xlab("Length") + ylab("Normalized counts (log2 RPM)") +
        guides(col = guide_legend(ncol = 3), shape = guide_legend(ncol = 1)) +
        guides(col = F) +
        theme_classic() +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

      p <- p +
        geom_hline(yintercept = median(subset(counts, E2 %in% molar)$count_norm), color="darkgrey")

      print(p)
    }, error=function(e){})
  }
dev.off()
