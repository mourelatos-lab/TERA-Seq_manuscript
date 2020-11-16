#!/usr/bin/env Rscript
#
# Attempt to reassign reads to gene/transcript reference for ONT data
# It's using the main "feature" of direct RNA-Seq ONT data which are more complete at the 3' end so
#   it goes from the end and assigns the reads by exon overlaps
#
# IMPORTANT: Right now I am taking as the first/last exon only first/last exon MAPPED, not in the reference annotation (we have different criteria for first start and last end)
#            We can add first/last annotated if it behaves strange
#

suppressMessages(library("optparse"))
library("stringr")
suppressMessages(library("dplyr"))
#library("reshape2") # Needed later
#library(ggplot2) # Needed later

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ibed"), type = "character",
    help = "Input table end-to-end (BED format from bedtools bamtobed)", metavar = "File"),
  make_option(
    c("-b", "--ibed_split"), type = "character",
    help = "Input table split (BED format from bedtools bamtobed & intersect)", metavar = "File"),
  make_option(
    c("-u", "--subset_rem"), type = "character", default=NULL,
    help = "Read name to transcript/gene id mapping table for UNAMBIGUOUS mappings (TAB; no header)", metavar = "File"),
  make_option(
    c("-a", "--trans_list"), type = "character", default=NULL,
    help = "Read name to transcript/gene id mapping table for ALL mappings (TAB; no header)", metavar = "File"),
  make_option(
    c("-m", "--molarity"), type = "character",
    help = "Gene/transcript id to molarity (or any other group) mapping.", metavar = "File"),
  make_option(
    c("-o", "--outdir"), type = "character",
    help = "Name of the output directory", metavar = "Directory"),
  make_option(
    c("-e", "--tolerance_end"), type = "integer", default=20,
    help = "Tolerance (bp) for mapping around last exon end boundary. Default: 20.", metavar = "Num"),
  make_option(
    c("-p", "--tolerance_splice"), type = "integer", default=10,
    help = "Tolerance (bp) for mapping around internal exon start/end boundary (splice). Default: 10.", metavar = "Num"),
  make_option(
    c("-v", "--overlap"), type = "integer", default=3,
    help = "Minimum overlap (bp) between exon and mapping. Default: 3", metavar = "Num")
)
opt = parse_args(OptionParser(option_list = option_list))

ibed<-opt$ibed
ibed_split<-opt$ibed_split
molarity<-opt$molarity
outdir<-opt$outdir
tolerance_end<-opt$tolerance_end
tolerance_splice<-opt$tolerance_splice
min_overlap<-opt$overlap

### TESTING VARIABLES!!!
# print("USING TEST VARIABLES!!!")
# opt<-NULL
# ibed<-"/home/joppelt/projects/ribothrypsis/analysis/spikein/results/hsa.dRNASeq.SIRV.polyA.REL5.long.2/unambig/reads.1.sanitize.toGenome.sorted.bed"
# ibed_split<-"/home/joppelt/projects/ribothrypsis/analysis/spikein/results/hsa.dRNASeq.SIRV.polyA.REL5.long.2/unambig/reads.1.sanitize.toGenome.sorted.split.overlaps.bed"
# # opt$subset_rem<-"E:/Results/zissimos/spikeins/sirv/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.unambig-tab.txt" # Read names to subset from the original bed (=remove); for example those uniquely assignable to transcripts from transcriptome bam; can be NULL or empty
# opt$subset_rem<-NULL; opt$trans_list<-NULL
# # opt$trans_list<-"E:/Results/zissimos/spikeins/sirv/reads.1.sanitize.noribo.toTranscriptome-polya.sorted.trans_tab.txt" # list of reads assigned to genes from transcriptome bam, we'll use these to "purify" the assigned reads to avoid overmapping issues; use ALL mappings, including secondary ones
# molarity<-"/home/joppelt/projects/ribothrypsis/data/spikein/sirv/molarity-E2.tab"
# outdir<-"/home/joppelt/projects/ribothrypsis/analysis/spikein/results/hsa.dRNASeq.SIRV.polyA.REL5.long.2/reassign/genome"
# # Tolerance: Based on IGV we see the basecalling+alignment is not always 100%
# # The downstream tolerance is important because of mapping of the RTA and problem with deciding how to assign this
# tolerance_end<-20 # maximum tolerance to count the read assignment; right now both up and downstream are the same
# tolerance_splice<-10 # maximum tol. in internal splice sites; igv shows this is very precise so we can go low
# min_overlap<-3 # Minimum overlap
# print("USING TEST VARIABLES!!!")
### TESTING VARIABLES!!!

####################################################################################################

bed<-read.table(ibed, sep="\t", stringsAsFactors = F, header=T)
bed_split<-read.table(ibed_split, sep="\t", stringsAsFactors = F, header=T)

if(!(is.null(opt$trans_list))){
  trans_refine<-read.table(opt$trans_list, sep="\t", header=F, stringsAsFactors = F)
  colnames(trans_refine)<-c("transcript_id", "name")
}
molar<-read.table(molarity, sep="\t", header=T, stringsAsFactors = F)

bed_assign<-NULL # Initiate final table

####################################################################################################
### Remove already assigned reads
if(!(is.null(opt$subset_rem))){
  print(paste("Removing reads from file:", opt$subset_rem))
  sub_rem<-read.table(opt$subset_rem, sep="\t", stringsAsFactors = F, header=F)
  colnames(sub_rem)<-c("transcript_id", "name")
  bed<-bed[!(bed$name %in% sub_rem$name),]
  bed_split<-bed_split[!(bed_split$name %in% sub_rem$name),]
  if(!(is.null(opt$trans_list))){
    trans_refine<-trans_refine[!(trans_refine$name %in% sub_rem$name), ]
  }
}

### Remove non-unique alignments
# From BED we can see that the duplicated alignment is usually very short but we would have to solve this at the BAM level
print(paste("Removing duplicated reads. Duplicated reads are the reads which are >1 in the ibed file:", ibed))
dupl<-bed$name[duplicated(bed$name)]
bed<-bed[!(bed$name %in% dupl), ]
bed_split<-bed_split[!(bed_split$name %in% dupl), ]

bed_split$size_ex<-bed_split$end_ex-bed_split$start_ex # Exon sizes, we might use it later for short exons filtering

bed_split$transcript_id <- bed_split$name_ex %>%
  str_split_fixed(pattern = "_", n = 2) %>%
  as.data.frame() %>%
  select(1) %>%
  unlist() # Add transcript name

####################################################################################################
####################################################################################################
### Get distances from starts and ends
bed_split$starts <- bed_split$start - bed_split$start_ex
bed_split$ends <- bed_split$end - bed_split$end_ex

### Split reads to plus and minus
# Plus
bed_split_readId <- bed_split %>%
  filter(strand=="+") %>%
  select(name, start, end) %>%
  unique() %>%
  group_by(name) %>%
  arrange(start) %>%
  mutate(read_id = row_number()) %>%
  mutate(read_max = max(read_id)) %>%
  ungroup()

bed_split_exonId <- bed_split %>%
  filter(strand=="+") %>%
  select(transcript_id, name_ex, start_ex, end_ex) %>%
  unique() %>%
  group_by(transcript_id) %>%
  arrange(start_ex) %>%
  mutate(exon_id = row_number()) %>%
  mutate(exon_max = max(exon_id)) %>%
  ungroup()

bed_plus <- bed_split %>%
  filter(strand=="+") %>%
  left_join(bed_split_readId) %>%
  left_join(bed_split_exonId) %>%
  filter(transcript_id!=".") # Remove non-overlapping reads (we already counted them)

# Remove reads where at least one part doesn't overlap the transcript, or to be more precise - keep only where all the segments map to the gene(s)
bed_plus <- bed_plus %>%
  select(start, end, name, transcript_id) %>%
  unique() %>%
  group_by(name, transcript_id) %>%
  mutate(read_in_trans=n()) %>%
  ungroup() %>%
  right_join(bed_plus) %>%
  filter(read_max==read_in_trans)

###
# Minus
bed_split_readId <- bed_split %>%
  filter(strand=="-") %>%
  select(name, start, end) %>%
  unique() %>%
  group_by(name) %>%
  arrange(desc(end)) %>%
  mutate(read_id = row_number()) %>%
  mutate(read_max = max(read_id)) %>%
  ungroup()

bed_split_exonId <- bed_split %>%
  filter(strand=="-") %>%
  select(transcript_id, name_ex, start_ex, end_ex) %>%
  unique() %>%
  group_by(transcript_id) %>%
  arrange(desc(end_ex)) %>%
  mutate(exon_id = row_number()) %>%
  mutate(exon_max = max(exon_id)) %>%
  ungroup()

bed_minus <- bed_split %>%
  filter(strand=="-") %>%
  left_join(bed_split_readId) %>%
  left_join(bed_split_exonId) %>%
  filter(transcript_id!=".") # Remove non-overlapping reads (we already counted them)

# Remove reads where at least one part doesn't overlap the transcript, or to be more precise - keep only where all the segments map to the gene(s)
bed_minus<-bed_minus %>%
  select(start, end, name, transcript_id) %>%
  unique() %>%
  group_by(name, transcript_id) %>%
  mutate(read_in_trans=n()) %>%
  ungroup() %>%
  right_join(bed_minus) %>%
  filter(read_max==read_in_trans)

####################################################################################################
### Get overlaps and score

assign_reads<-function(bed_run){
#  bed_run<-bed_minus
  # IMPORTANT: We have already solved this by changing the read_max==read_in_trans filtering...but we still keep this here just in case
  # Get reads which do not end with the first exon and span further outside the first exon and remove them
  # THIS NOTE NOT VALID ANYMORE: This removes also "splicing" within the exon - if aligner decides it won't to insertion/deletetion but would to splicing instead this remove such reads. You have to trust the aligner!
  #   Now, we check if the first overlap start or the last end is inside (ok) - these should be the spliced reads
  fp<-bed_run %>%
    filter(exon_id==1) %>%
    filter(read_id!=1) %>%
    filter(starts<0) %>% # if it's inside and it's not the first fragment, the starts difference is positive (we want to keep it)
    select(name, transcript_id)
  # Do the same as above but for last exon and fragment
  fp<-bed_run %>%
    filter(exon_id==exon_max) %>%
    filter(read_id!=read_max) %>%
    filter(ends>0) %>% # if it's inside and it's not the last fragment, the ends difference is negative (we want to keep it)
    select(name, transcript_id) %>%
    rbind(fp) %>%
    unique()

  bed_run <- bed_run[!(paste(bed_run$name, bed_run$transcript_id, sep=".") %in% unique(paste(fp$name, fp$transcript_id, sep="."))), ]

  # Add row index for easier filtering
  bed_run <- bed_run %>%
    mutate(row_id = row_number())

  # Filter min. overlaps
  bed_run <- bed_run %>%
    filter(overlap>=!!min_overlap)

  # Conditionaly change start of first and end of last exon within tolerance_end to tolerance_splice
  bed_run <- bed_run %>%
    mutate(ends = if_else(condition = read_id==read_max & exon_id==exon_max & abs(ends)<=tolerance_end,
                          true = as.integer(tolerance_splice), false = ends))

  bed_run <- bed_run %>%
    mutate(starts = if_else(condition = read_id==1 & exon_id==1 & abs(starts)<=tolerance_end,
                          true = as.integer(tolerance_splice), false = starts))

  ####
  # Check if a read goes all the way end to end (ete) of the transcript
  bed_run$ete_first<-bed_run$read_id==1 & bed_run$exon_id==1 & abs(bed_run$starts) <= tolerance_splice # Start of first
  bed_run$ete_last<-bed_run$read_id==bed_run$read_max & bed_run$exon_id==bed_run$exon_max & abs(bed_run$ends) <= tolerance_splice # End of last

  bed_run<-bed_run %>%
#    filter((ete_first + ete_last)!=0) %>%
    group_by(name, transcript_id) %>%
    mutate(etes=n()) %>%
    ungroup()

  # Get only end-to-end mappings
  bed_final_tmp <- bed_run %>%
    filter(etes==2) %>%
    select(name, transcript_id)

  bed_final<-bed_final_tmp[!(bed_final_tmp$name %in% bed_final_tmp$name[duplicated(bed_final_tmp$name)]), ] # Get only unique

  # Remove assigned reads
  bed_run <- bed_run %>%
    filter(!(name %in% unique(bed_final$name)))

  # Compare all starts and all ends (we already have changed last end and first start) - we basicaly give points for each match
  bed_run$full_match <- abs(bed_run$ends)<=tolerance_splice & abs(bed_run$starts)<=tolerance_splice # both ends matching = full-match

  # Compare partial matches
  bed_run$partial_match_end <- abs(bed_run$ends)<=tolerance_splice # Partial matches
  bed_run$partial_match_start <- abs(bed_run$starts)<=tolerance_splice

  bed_run <- bed_run %>%
    group_by(name, transcript_id) %>%
    mutate(full_matches=sum(full_match)) %>%
    mutate(partial_match_ends=sum(partial_match_end)) %>%
    mutate(partial_match_starts=sum(partial_match_start)) %>%
    mutate(overlaps=sum(overlap)) %>%
    ungroup()

  bed_run$partial_matches <- bed_run$partial_match_ends + bed_run$partial_match_starts # ?Sum them?

  bed_run<-bed_run %>%
    group_by(name) %>%
    mutate(full_matches_max = max(full_matches)) %>%
    mutate(partial_matches_max = max(partial_matches)) %>%
    mutate(overlaps_max = max(overlaps)) %>%
    ungroup()

  # If most of the full matches is unique, move it to a final table
  bed_final_tmp <- bed_run %>%
    dplyr::filter(full_matches > 0) %>%
    dplyr::filter(full_matches == full_matches_max) %>%
    select(name, transcript_id) %>%
    unique()

  bed_final<-bed_final_tmp[!(bed_final_tmp$name %in% bed_final_tmp$name[duplicated(bed_final_tmp$name)]),] # Get only unique max

  # Remove assigned reads
  bed_run <- bed_run %>%
    filter(!(name %in% unique(bed_final$name)))

  # Compare the partial matches for reads with non-unique max. full-matches
  bed_final_tmp <- bed_run %>%
    dplyr::filter(full_matches == full_matches_max) %>%
    dplyr::filter(partial_matches>0) %>%
    dplyr::filter(partial_matches == partial_matches_max) %>%
    select(name, transcript_id) %>%
    unique()

  # Add to the final table if unique
  bed_final<-rbind(bed_final, bed_final_tmp[!(bed_final_tmp$name %in% bed_final_tmp$name[duplicated(bed_final_tmp$name)]),]) # Get only unique max

  # Remove assigned reads
  bed_run <- bed_run %>%
    filter(!(name %in% unique(bed_final$name)))

  # If full-matches and partial-matches cannot decide look at the overlap
  bed_final_tmp <- bed_run %>%
    dplyr::filter(full_matches == full_matches_max) %>%
    dplyr::filter(partial_matches>0) %>%
    dplyr::filter(partial_matches == partial_matches_max) %>%
    dplyr::filter(overlaps>0) %>%
    dplyr::filter(overlaps == overlaps_max) %>%
    select(name, transcript_id) %>%
    unique()

  bed_final<-rbind(bed_final, bed_final_tmp[!(bed_final_tmp$name %in% bed_final_tmp$name[duplicated(bed_final_tmp$name)]),])

  # Remove assigned reads
  bed_run <- bed_run %>%
    filter(!(name %in% unique(bed_final$name)))

  return(bed_final)
}

bed_final_plus<-assign_reads(bed_run = bed_plus)
bed_final_minus<-assign_reads(bed_run = bed_minus)

bed_final<-rbind(bed_final_plus, bed_final_minus)

###
# Get the remaining reads and decide what to do with them
bed_plus <- bed_plus %>%
  filter(!(name %in% unique(bed_final$name)))

bed_minus <- bed_minus %>%
  filter(!(name %in% unique(bed_final$name)))

##############################################
###
# From the back - if unique overlap end -> final; if unique overlap start -> final; otherwise keep note and move to next exon

# Since we already changed the ORDER (NOT coordinates of start/end) of exons and read fragments for - transcripts to look like + we can make one function both both
#   We have to change columns starts/ends for - transcripts so we can use the function. If we keep it we compare not the end of the last exon but the start in the first round and
#     then start of the last exon instead of end of the last exon which is not what we want

assign_reads2<-function(bed_run){
  for(iter in 0:(max(bed_run$exon_max)-1)){
    print(paste("Working on exon/block (from the end):", iter))

    if(iter==0){
      tolerance<-tolerance_end
      print(paste("Going with reads:", length(unique(bed_run$name))))
    }else{
      tolerance<-tolerance_splice
      print(paste("Going with reads:", length(unique(bed_run$name))))
    }

    # Get the last exon-last segment pairs and get only those which have end-end <= tolerance_splice
    bed_iter <- bed_run %>%
      filter(exon_max>!!iter) %>% # Test
      filter(exon_id==(exon_max-!!iter)) %>%
      filter(read_max>!!iter) %>% # Test
      filter(read_id==(read_max-!!iter)) %>%
      ungroup() %>%
      filter(abs(ends)<=tolerance)

    bed_iter2 <- bed_iter %>%
      select(name, transcript_id) %>%
      unique()

    # Get unique overlaps only
    if(iter==0){
      bed_final <- bed_iter2[!(bed_iter2$name %in% bed_iter2$name[duplicated(bed_iter2$name)]),]
    }else{
      bed_final <- rbind(bed_final, bed_iter2[!(bed_iter2$name %in% bed_iter2$name[duplicated(bed_iter2$name)]),])
    }

    # Remove assigned reads
    bed_iter <- bed_iter %>%
      filter(!(name %in% unique(bed_final$name)))

    bed_iter <- bed_iter %>%
      filter(abs(starts)<=tolerance_splice)

    bed_iter2 <- bed_iter %>%
      select(name, transcript_id) %>%
      unique()

    bed_final <- rbind(bed_final, bed_iter2[!(bed_iter2$name %in% bed_iter2$name[duplicated(bed_iter2$name)]),])

    # Remove assigned reads
    bed_iter <- bed_iter %>%
      filter(!(name %in% unique(bed_final$name)))

    # Get only pairs passing the previous round
    bed_run <- bed_run[paste(bed_run$name, bed_run$transcript_id, sep=".") %in% unique(paste(bed_iter$name, bed_iter$transcript_id, sep=".")), ]
  }
    return(bed_final)
}

bed_final_plus2<-assign_reads2(bed_run = bed_plus)
# We have to swap starts/ends in minus (-) transcripts otherwise we don't compare what we want (please see the note above the function)
bed_minus_swapped<-bed_minus
bed_minus_swapped$starts.tmp<-bed_minus_swapped$starts
bed_minus_swapped$starts<-bed_minus_swapped$ends
bed_minus_swapped$ends<-bed_minus_swapped$starts.tmp
bed_minus_swapped$starts.tmp<-NULL
bed_final_minus2<-assign_reads2(bed_run = bed_minus_swapped)

bed_final<-rbind(bed_final, rbind(bed_final_plus2, bed_final_minus2))

###
# Get the remaining reads and decide what to do with them
bed_plus <- bed_plus %>%
  filter(!(name %in% unique(bed_final$name)))

bed_minus <- bed_minus %>%
  filter(!(name %in% unique(bed_final$name)))

################################################################################
### Refine assigned reads

bed_final<-as.data.frame(bed_final)
bed_final$transcript_id<-as.character(bed_final$transcript_id) # Remove factor

dir.create(outdir, recursive = T)

# Get unassigned reads for debugging
unassigned<-bed_split %>%
  filter(!(name %in% unique(bed_final$name))) %>%
  select(transcript_id, name)
unassigned$transcript_id<-"SIRVunassigned"
unassigned <- unique(unassigned)

print(paste("Total reads:", length(unique(bed_split$name))))
print(paste("Total assigned:", length(unique(bed_final$name))))
print(paste("Total assigned percentage:", round(length(unique(bed_final$name))/(length(unique(bed_split$name))/100), 2)))
print(paste("Total unassigned:", length(unique(unassigned$name))))

write.table(x = bed_final %>% select(transcript_id, name) %>% rbind(unassigned), file=paste(outdir, "assign_read-to-trans.tsv", sep="/"), sep="\t", row.names=F, quote = F)

tab_compare<-bed_final %>% left_join(molar)

tab_compare <- tab_compare %>%
  group_by(transcript_id) %>%
  mutate(total=n()) %>%
  ungroup()

tab_comp_m<-tab_compare %>%
  select(transcript_id, molarity, total)
tab_comp_m$molarity<-factor(tab_comp_m$molarity, levels=sort(unique(tab_comp_m$molarity)))

sink(paste(outdir, "assign_summary.txt", sep="/"))
  tab_compare %>%
    group_by(molarity) %>%
    summarize(median=median(total))
sink()

library(ggplot2)

pdf(paste(outdir, "assign_summary.pdf", sep="/"))
  ggplot(data = tab_comp_m, aes(x=molarity, y=total, colour=molarity)) +
    theme_bw() +
    geom_boxplot(position=position_dodge(0.8), width = 0.8) +
    geom_smooth(aes(x=as.integer(molarity), y=total, color=molarity, fill=molarity), method=lm, se = F)

  ggplot(data = tab_comp_m, aes(x=molarity, y=total, colour=molarity)) +
    theme_bw() +
    geom_violin(position = position_dodge(width = 0.9))
dev.off()
