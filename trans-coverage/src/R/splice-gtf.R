#!/usr/bin/env Rscript
#
# Get coordinates of spliced exons (and start/stop codon if present) from gtf like they would look like in cDNA
#
# Note: Sometimes, the 3UTR and 5UTR are not annotated as exon (they should be) and right now, we ignore this. This might
#	cause a discrepancy at the ends of the transcripts (transcripts ends vs exon ends)
#

suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("dplyr"))

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-g", "--gtf"), type = "character",
    help = "GTF file with transcripts to be spliced (only exons considered).", metavar = "File"),
  make_option(
    c("-s", "--subset"), type = "character", default="",
    help = "List of transcripts to subset from the GTF (one transcript per line). Doesn't work for long lists (better use grep on gtf before). Default: use all transcripts.", metavar = "File"),
  make_option(
    c("-o", "--ofile"), type = "character",
    help = "Output file for the spliced transcript exon (BED).", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

### TESTING VARIABLES!!!
# print("USING TEST VARIABLES!!!")
# opt<-NULL
# opt$gtf<-"test.gtf"
# opt$subset<-""
# # opt$subset<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/transcripts_select.txt"
# opt$ofile<-"/home/joppelt/projects/ribothrypsis/analysis/transcript-coverage/results/transcripts_select-ultimate.spliced.bed"
# print("USING TEST VARIABLES!!!")
### TESTING VARIABLES!!!

#coords_main<-read.table(opt$gtf, sep="\t", header=F, stringsAsFactors = F) # slow
coords_main<-read_tsv(opt$gtf, col_names = F, comment = "#", col_types = cols(X1 = col_character())) # better, but needs readr
colnames(coords_main)<-gsub("^X", "V", colnames(coords_main))

coords_main <- coords_main %>%
  filter(V3 %in% c("start_codon", "stop_codon", "exon"))

# Subset the GTF
if(nchar(opt$subset)>1){
  trans_subset<-read.table(opt$subset, sep="\t", header=F, stringsAsFactors = F)

  coords_main <- coords_main[grep(paste(trans_subset$V1, collapse="|"),
       coords_main$V9), ]
  coords_main$transcript_id<-sapply(strsplit(sapply(strsplit(coords_main$V9, "transcript_id "), "[[", 2), "; "), "[[", 1) # Get only transcript ids
  coords_main$transcript_id<-gsub("\"", "", coords_main$transcript_id)

  coords_main <- coords_main %>%
    filter(transcript_id %in% trans_subset$V1)
}else{
  coords_main$transcript_id<-sapply(strsplit(sapply(strsplit(coords_main$V9, "transcript_id "), "[[", 2), "; "), "[[", 1) # Get only transcript ids
  coords_main$transcript_id<-gsub("\"", "", coords_main$transcript_id)
}

coords_main$len<-coords_main$V5-coords_main$V4 # Get length

coords_all<-NULL

# plus transcripts
if(sum(coords_main$V7=="+")>0){
  coords_sub<-coords_main %>%
    filter(V7=="+")

  for(id in unique(coords_sub$transcript_id)){
#    print(id)

    coords<-coords_sub %>%
      filter(transcript_id == id)

    # BED formated - 0 based
    # Check if start/stop codon is annotated and if it adjust the start position; stop codon adjustment not used
    adj_start<-FALSE
    if("start_codon" %in% coords$V3){
      adj_start<-TRUE
    }
    # adj_stop<-FALSE
    # if("stop_codon" %in% coords$V3){
    #   adj_stop<-TRUE
    # }

    coords <- coords %>%
      filter(V3 == "exon") %>%
      arrange(transcript_id, V4)

    # First exon + adjustment if necessary
    if(adj_start==TRUE){
      coords[1, "start"] <- 3
    }else{
      coords[1, "start"] <- 0
    }
    coords[1, "end"] <- coords[1, "start"] + coords[1, "len"]

    # All the other exons (if more than one)
    if(nrow(coords)>1){ # Multi exon transcripts
      for(line in 2:nrow(coords)){
        coords[line, "start"] <- coords[line-1, "end"]
        coords[line, "end"] <- coords[line-1, "end"] + coords[line, "len"]
        coords <- coords[, c("transcript_id", "start", "end", "transcript_id", "len", "V7")] # Make BED
        colnames(coords)[4]<-"name"
        coords$name <- paste(coords$name, "exon", sep="-")
      }
    }else if(nrow(coords)==1){ # Single exon
      coords <- coords[, c("transcript_id", "start", "end", "transcript_id", "len", "V7")] # Make BED
      colnames(coords)[4]<-"name"
      coords$name <- paste(coords$name, "exon", sep="-")
    }

    if(is.null(coords_all)){
      coords_all<-coords
    }else{
      coords_all<-rbind(coords_all, coords)
    }
  }
}

# Minus strand transcripts
if(sum(coords_main$V7=="-")>0){
  coords_sub<-coords_main %>%
    filter(V7 == "-")

  for(id in unique(coords_sub$transcript_id)){
#    print(id)

    coords<-coords_sub %>%
      filter(transcript_id == id)

    # BED formated - 0 based
    # Check if start/stop codon is annotated and if it adjust the start position
    adj_start<-FALSE
    if("start_codon" %in% coords$V3){
      adj_start<-TRUE
    }
    adj_stop<-FALSE
    if("stop_codon" %in% coords$V3){
      adj_stop<-TRUE
    }

    coords <- coords %>%
      filter(V3 == "exon") %>%
      arrange(transcript_id, -V5)

    # First exon + adjustment if necessary
    if(adj_start==TRUE){
      coords[1, "start"] <- 3
    }else{
      coords[1, "start"] <- 0
    }
    coords[1, "end"] <- coords[1, "start"] + coords[1, "len"]

    # All the other exons (if more than one)
    if(nrow(coords)>1){ # Multi exon transcripts
      for(line in 2:nrow(coords)){
        coords[line, "start"] <- coords[line-1, "end"]
        coords[line, "end"] <- coords[line-1, "end"] + coords[line, "len"]
        coords <- coords[, c("transcript_id", "start", "end", "transcript_id", "len", "V7")] # Make BED
        colnames(coords)[4]<-"name"
        coords$name <- paste(coords$name, "exon", sep="-")
      }
    }else if(nrow(coords)==1){ # Single exon
      coords <- coords[, c("transcript_id", "start", "end", "transcript_id", "len", "V7")] # Make BED
      colnames(coords)[4]<-"name"
      coords$name <- paste(coords$name, "exon", sep="-")
    }
    
    if(is.null(coords_all)){
      coords_all<-coords
    }else{
      coords_all<-rbind(coords_all, coords)
    }
    
  }
}

#write.table(x = coords_all, file = opt$ofile,
#            quote = F, col.names=F, row.names=F, sep="\t")
readr::write_delim(x = coords_all, path = opt$ofile, delim = "\t", na = "NA", append = FALSE,
            col_names = F) # , quote_escape = F
