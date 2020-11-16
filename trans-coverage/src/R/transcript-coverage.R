#!/usr/bin/env Rscript
#
# Plot RCS (enolase) and chromosome/transcript coverage
#

# Rounding for the coords to the plot based on a number
mround <- function(x,base){
  base*ceiling(x/base)
}

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Input table - read ends (BED format from bedtools bamtobed). For stdin use 'stdin'.", metavar = "File"),
  make_option(
    c("-d", "--depth"), type = "character", default="", 
    help = "Input table - coverage (tab format from samtools depth) with colnames.", metavar = "File"),
  make_option(
    c("-u", "--illumina"), type = "character", default = "",
    help = "Add Illumina coverage (same format as samtools depth).", metavar = "File"),  
  make_option(
    c("-o", "--suffix"), type = "character", default = "",
    help = "Output file SUFFIX The individual plots will be named accoring to transcript id + suffix.", metavar = "File"),  
  make_option(
    c("-c", "--coords"), type = "character", default="",
    help = "BED with coordinates to plot if the input BED is genomic. If empty plots by chromosome names in first column.", metavar = "File"),
  make_option(
    c("-b", "--bin"), type = "integer", default=5,
    help = "Bin size for the coverage. Default: 5.", metavar = "Num"),
  make_option(
    c("-s", "--subset"), type = "integer", default=9999,
    help = "Subset the reads for each chromosome/transcript/line in coords bed. Default: 9999.", metavar = "Num"),
  make_option(
    c("-l", "--addlines"), type = "character", default="",
    help = "Add lines (exon boundaries, for example) to the plots (spliced BED format).", metavar = "File"),
  make_option(
    c("-p", "--addpoints"), type = "character", default="",
    help = "Add points (CAGE, for example) to the plots (spliced BED format).", metavar = "File"),
  make_option(
    c("-p2", "--addpoints2"), type = "character", default="",
    help = "Add points (NET-CAGE, for example) to the plots (spliced BED format).", metavar = "File"),  
  make_option(
    c("-r", "--addregions"), type = "character", default="",
    help = "Add regions (TSS), for example) to the plots (spliced BED format).", metavar = "File"),  
  make_option(
    c("-n", "--normdepth"), action="store_true", default=FALSE,
    help = "Normalize to library specific depth only. Default is to normalize do depth and all-libraries maximum covered region.")
)
opt = parse_args(OptionParser(option_list = option_list))

### TESTING VARIABLES
#print("USING TEST VARIABLES!!!!")
#opt<-NULL
#opt$ifile<-"E:/Results/zissimos/coverage/coverage.100.w_wo_rel5.bed"
#opt$coords<-""
#opt$bin=5
#opt$subset<-10000
#opt$addlines<-"E:/Results/zissimos/coverage/transcripts_select.CIP.decap.100.spliced.bed"
#opt$addpoints<-"E:/Results/zissimos/coverage/HeLa.repAll.hg38.ctss.genom-to-trans-ByExon.bed"
#opt$addpoints2<-"E:/Results/zissimos/coverage/RepAll-HeLaS3-NETCAGE-1M.hg38.genom-to-trans-ByExon.bed"
#opt$suffix<-""
#opt$addregions<-"E:/Results/zissimos/coverage/atlas.clusters.hg38.2-0.genom-to-trans-ByExon.bed"
#print("USING TEST VARIABLES!!!!")
### TESTING VARIABLES

input<-opt$ifile # legacy
coords_in<-opt$coords
bin_by<-opt$bin
subset<-opt$subset

if(opt$ifile == "stdin"){
  coverage_in <- read.table(file("stdin"), sep="\t", header=T, stringsAsFactors = T)
}else{
  coverage_in <- read.table(opt$ifile, sep="\t", header=T, stringsAsFactors = F)
}

coverage_depth<-NULL
if(nchar(opt$depth)>0){
    coverage_depth <- read.table(opt$depth, sep="\t", header=T, stringsAsFactors = F)
}

illumina_depth<-NULL
if(nchar(opt$illumina)>0){
    illumina_depth <- read.table(opt$illumina, sep="\t", header=T, stringsAsFactors = F)
}

if(nchar(coords_in)>0){
  coords_main<-read.table(coords_in, sep="\t", header=F, stringsAsFactors = F)
}else{ # Make "fake" bed so we can use both transcriptome and genome to make the coverage plots
  coords_main<-NULL
  coords_main<-as.data.frame(unique(coverage_in$chr))
  colnames(coords_main)<-"V4"
}

if(nchar(opt$addlines)>0){
  add_lines<-read.table(opt$addlines, sep="\t", header=F, stringsAsFactors = F)
}

if(nchar(opt$addpoints)>0){
  add_points<-read.table(opt$addpoints, sep="\t", header=F, stringsAsFactors = F)
}

if(nchar(opt$addpoints2)>0){
  add_points2<-read.table(opt$addpoints2, sep="\t", header=F, stringsAsFactors = F)
}

if(nchar(opt$addregions)>0){
  add_regions<-read.table(opt$addregions, sep="\t", header=F, stringsAsFactors = F)
}

for(chr_name in unique(coords_main$V4)){
#  chr_name<-"ENST00000215754" # testing
  print(chr_name)

  coords<-coords_main %>%
    filter(V4 == chr_name)
    
  # Assuming these are genomic coordinates get the region with a bit of tolerance on the sides
  if(nchar(coords_in)>0){
    coverage_main <- coverage_in %>%
      filter(chr%in%coords$V1) %>%
      filter(start>=coords$V2-50) %>%
      filter(end<=coords$V3+50) 
  }else{
    coverage_main <- coverage_in %>%
      filter(chr == chr_name)
  }
  
  if(nchar(opt$addlines)>0){
    add_line<-add_lines %>% 
      filter(V1 == chr_name)
  }else{
    add_line<-NULL
  }
  
  if(nchar(opt$addpoints)>0){
    add_point<-add_points %>% 
      filter(V1 == chr_name)
  }else{
    add_point<-NULL
  }
  
  if(nchar(opt$addpoints2)>0){
    add_point2<-add_points2 %>% 
      filter(V1 == chr_name)
  }else{
    add_point2<-NULL
  }
  
  if(nchar(opt$addregions)>0){
    add_region<-add_regions %>% 
      filter(V1 == chr_name)
  }else{
    add_region<-NULL
  }
    
  all_forplot<-NULL
  all_coverage<-NULL
  
  len_for_norm<-max(coverage_main$end)-min(coverage_main$start) # Get covered transcript region for normalization - all libraries together
  len_for_norm<-len_for_norm/opt$bin # Make it bins instead of coords
  
  for(lib_select in unique(coverage_main$Library)){
    print(lib_select)
    
    coverage<-coverage_main %>%
      filter(Library==lib_select)
      
    # Subset to get random subset reads
    if(nrow(coverage)>subset){
      coverage<-sample_n(coverage, subset, replace=F)
    }
    
    # Get aligned lengths
    print("Mean and median read length.")
    print(mean(coverage$end-coverage$start))
    print(median(coverage$end-coverage$start))
      
    full_len<-coverage
    full_len$coord<-paste(full_len$chr, full_len$start, full_len$end, sep=".")
    full_len<-table(full_len$coord)
    
    sum(full_len) # How many reads we have
    full_len[full_len==max(full_len)] # what is the actual full length? probably the most frequent one (+/-)
    
    p5<-as.data.frame(table(coverage$start))
    p3<-as.data.frame(table(coverage$end))
            
    p5<-p5 %>% 
      mutate(end="5p") %>%
      mutate(Library=unique(coverage$Library))
    p3<-p3 %>% 
      mutate(end="3p") %>%
      mutate(Library=unique(coverage$Library))
    
    p5$Var1<-as.numeric(as.character(p5$Var1))
    p5$Freq<-as.numeric(as.character(p5$Freq))
    p3$Var1<-as.numeric(as.character(p3$Var1))
    p3$Freq<-as.numeric(as.character(p3$Freq))  
    
    # Make bins by 5 nt
    p5$bin<-cut(p5$Var1, breaks = seq(floor(min(coverage$start)/10)*10, ceiling(max(coverage$end)/10)*10, by=bin_by), 
                labels = seq((floor(min(coverage$start)/10)*10)+bin_by, ceiling(max(coverage$end)/10)*10, by=bin_by), include.lowest = T)
    p3$bin<-cut(p3$Var1, breaks = seq(floor(min(coverage$start)/10)*10, ceiling(max(coverage$end)/10)*10, by=bin_by), 
                labels = seq((floor(min(coverage$start)/10)*10)+bin_by, ceiling(max(coverage$end)/10)*10, by=bin_by), include.lowest = T)
      
    # Sum bins and add percentage
    p5.cover<-p5 %>%
      group_by(Library, bin) %>%
      summarize(Coverage=sum(Freq)) %>%
      group_by(Library) %>%
      mutate(Percentage = (Coverage/sum(Coverage)*100)) %>%
      mutate(End="5P")
    
    p3.cover<-p3 %>%
      group_by(Library, bin) %>%
      summarize(Coverage=sum(Freq)) %>%
      group_by(Library) %>%
      mutate(Percentage = (Coverage/sum(Coverage)*100)) %>%
      mutate(End="3P")
    
    forplot<-rbind(p5.cover, p3.cover)
    
    forplot$bin<-as.numeric(as.character(forplot$bin))
    
    if(is.null(all_forplot)){
      all_forplot<-forplot
    }else{
      all_forplot<-rbind(all_forplot, forplot)
    }
    
    ### Get "classic" coverage
    all_ranges<-NULL
    
    # Get coverage from bamtobed if samtools depth is not available
    if(is.null(coverage_depth)){
        for(rows in 1:nrow(coverage)){
          all_ranges<-c(all_ranges, seq(coverage$start[rows], coverage$end[rows], by=1))
        }
        
        coverage_ranges<-as.data.frame(table(all_ranges))
        coverage_ranges[]<-lapply(coverage_ranges, function(x) as.numeric(as.character(x))) # Convert to numeric
        
        coverage_ranges<-coverage_ranges %>%
            rename(Position="all_ranges", Coverage="Freq")
    }else{ # Get coverage if we have samtools depth
        coverage_ranges <- coverage_depth %>% 
            filter(chr == !!chr_name) %>%
            filter(Library == !!lib_select) %>%
            rename(Position="pos", Coverage="depth")
    }

    # Make fake zero coverage in case we don't have any info about the chromosome
    if(nrow(coverage_ranges)==0){
      coverage_ranges[1, ]<-c(lib_select, chr_name, 1, 0)
      coverage_ranges$Position<-as.numeric(coverage_ranges$Position)
      coverage_ranges$Coverage<-as.numeric(coverage_ranges$Coverage)
    }
    
    coverage_ranges$Library<-lib_select
    
    coverage_ranges$bin<-cut(coverage_ranges$Position, breaks = seq(floor(min(coverage_ranges$Position)/10)*10, ceiling(max(coverage_ranges$Position)/10)*10, by=bin_by), 
        labels = seq((floor(min(coverage_ranges$Position)/10)*10)+bin_by, ceiling(max(coverage_ranges$Position)/10)*10, by=bin_by), include.lowest = T)
    
    coverage.cover<-coverage_ranges %>%
      group_by(Library, bin) %>%
      summarize(Coverage=sum(Coverage)) %>%
      group_by(Library) %>%
      mutate(Percentage = (Coverage/sum(Coverage)*100))
    
    if(opt$normdepth==FALSE){
      coverage.cover$Percentage<-coverage.cover$Percentage/len_for_norm # Normalize to length - experimental
    }
    
    coverage.cover$bin<-as.numeric(as.character(coverage.cover$bin))
    
    if(is.null(all_coverage)){
        all_coverage<-coverage.cover
    }else{
        all_coverage<-rbind(all_coverage, coverage.cover)
    }
  }  
  
  # Add Illumina if we have it
  if(!is.null(illumina_depth)){
      illumina_ranges <- illumina_depth %>% 
      filter(chr == !!chr_name) %>%
      rename(Position="pos", Coverage="depth")  
    
    # Fake coverage is illumina coverage is none
    if(nrow(illumina_ranges)==0){
      illumina_ranges[1, ]<-c("illumina-nocover", chr_name, 1, 0)
      illumina_ranges$Position<-as.numeric(illumina_ranges$Position)
      illumina_ranges$Coverage<-as.numeric(illumina_ranges$Coverage)
    }
    
    illumina_ranges$bin<-cut(illumina_ranges$Position, breaks = seq(floor(min(illumina_ranges$Position)/10)*10, ceiling(max(illumina_ranges$Position)/10)*10, by=bin_by), 
        labels = seq((floor(min(illumina_ranges$Position)/10)*10)+bin_by, ceiling(max(illumina_ranges$Position)/10)*10, by=bin_by), include.lowest = T)
    
    illumina.cover<-illumina_ranges %>%
      group_by(Library, bin) %>%
      summarize(Coverage=sum(Coverage)) %>%
      group_by(Library) %>%
      mutate(Percentage = (Coverage/sum(Coverage)*100))
    
    if(opt$normdepth==FALSE){
      illumina.cover$Percentage<-illumina.cover$Percentage/len_for_norm # Normalize to length - experimental
    }
    
    illumina.cover$bin<-as.numeric(as.character(illumina.cover$bin))

    all_coverage<-rbind(all_coverage, illumina.cover)
  }
  
  all_coverage$Library<-factor(all_coverage$Library, levels = unique(all_coverage$Library))  
  
  # Get min and max coords for the plots
  if(nchar(coords_in)>0){
    min_coord<-mround(coords$V2, bin_by)
    max_coord<-mround(coords$V3, bin_by)
  }else{
    min_coord<-bin_by
    max_coord<-mround(max(coverage_main$end), bin_by)
  }
  
  if(opt$ifile=="stdin"){
    if(nchar(opt$suffix)>0){
      pdf(paste0("results/compare", "/coverage-", coords$V4, "-", opt$suffix, ".pdf"))
    }else{
      pdf(paste0("results/compare", "/coverage-", coords$V4, ".pdf"))
    }
  }else{
    if(nchar(opt$suffix)>0){
      pdf(paste0(dirname(opt$ifile), "/coverage-", coords$V4, "-", opt$suffix, ".pdf"))
    }else{
      pdf(paste0(dirname(opt$ifile), "/coverage-", coords$V4, ".pdf"))
    }
  }  
    p<-ggplot(all_coverage, aes(x=bin, y=Percentage, color=Library)) +
      geom_line() +
      geom_vline(xintercept = c(min_coord, max_coord), color="grey", size=0.25, alpha=0.5)  + 
      geom_vline(xintercept = unique(c(mround(add_line$V2, bin_by), mround(add_line$V3, bin_by))), color="orange", linetype = "dashed", size=0.25, alpha=0.5)  + 
      theme_classic() +
      ggtitle(paste0("Coverage of ", coords$V4)) +
      xlab(paste0("Position (binned by ", opt$bin, " bp)")) +
      ylab(paste0("Percentage coverage per bin (binned by ", opt$bin, " bp)")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position="bottom") +
      guides(color=guide_legend(ncol=1))

    if(!is.null(add_point)){    
      if(nrow(add_point)>0){
        p<-p + 
          geom_point(data = add_point, aes(x=mround(V2, bin_by), y=-max(all_coverage$Percentage)/35), color="black", size=1.4, alpha=0.5)
      }  
    }
    
    if(!is.null(add_point2)){    
      if(nrow(add_point2)>0){
        p<-p + 
          geom_point(data = add_point2, aes(x=mround(V2, bin_by), y=-max(all_coverage$Percentage)/65), color="violet", size=1.4, alpha=0.5)
      }  
    }
    
    if(!is.null(add_region)){
      if(nrow(add_region)>0){
          p<-p + 
            geom_segment(data = add_region, aes(x=mround(V2, bin_by), xend=mround(V3, bin_by), y=-max(all_coverage$Percentage)/25, yend=-max(all_coverage$Percentage)/25), 
                         color="blue", size=1.4, alpha=0.5)
      }
    }
    
    print(p)
    
    # Whole plot
    p<-ggplot(all_forplot, aes(x=bin, y=Percentage, color=End)) +
      geom_point(alpha=0.5) +
      geom_vline(xintercept = c(min_coord, max_coord), color="grey", size=0.25, alpha=0.5)  +
      geom_vline(xintercept = unique(c(mround(add_line$V2, bin_by), mround(add_line$V3, bin_by))), color="orange", linetype = "dashed", size=0.25, alpha=0.5)  + 
      theme_classic() +
      ggtitle(paste0("Coverage of ", coords$V4)) +
      xlab(paste0("Position (binned by ", opt$bin, " bp)")) + 
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      facet_grid(Library ~ .) # original start:end is 451326 and 452640 but we are binning by 5 so we have to go to the correct bins

    if(!is.null(add_point)){
      if(nrow(add_point)>0){
        p<-p +
          geom_point(data = add_point, aes(x=mround(V2, bin_by), y=-max(all_forplot$Percentage)/25), color="black", alpha=0.5)
      }
    }
    
    if(!is.null(add_point2)){
      if(nrow(add_point2)>0){
        p<-p +
          geom_point(data = add_point2, aes(x=mround(V2, bin_by), y=-max(all_forplot$Percentage)/40), color="violet", alpha=0.5)
      }
    }
    
    if(!is.null(add_region)){
      if(nrow(add_region)>0){
        p<-p +
          geom_segment(data = add_region, aes(x=mround(V2, bin_by), xend=mround(V3, bin_by), y=-max(all_forplot$Percentage)/12, yend=-max(all_forplot$Percentage)/12), 
                       size=2, color="blue", alpha=0.5)
      }
    }
    
    print(p)
    
    # Zoom in to low-coverage points
    p<-ggplot(all_forplot, aes(x=bin, y=Percentage, color=End)) +
      geom_point(alpha=0.5) +
      geom_vline(xintercept = c(min_coord, max_coord), color="grey", size=0.25, alpha=0.5) + # original start:end is 451326 and 452640 but we are binning by 5 so we have to go to the correct bins
      geom_vline(xintercept = unique(c(mround(add_line$V2, bin_by), mround(add_line$V3, bin_by))), color="orange", linetype = "dashed", size=0.25, alpha=0.5)  + 
      coord_cartesian(ylim=c(-0.5, 10))  +
      ggtitle(paste0("Coverage of ", coords$V4, " - low-coverage positions (<=10)")) +
      xlab(paste0("Position (binned by ", opt$bin, " bp)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      facet_grid(Library ~ .)

    if(!is.null(add_point)){
      if(nrow(add_point)>0){
        p<-p +
          geom_point(data = add_point, aes(x=mround(V2, bin_by), y=-(10.5/30)), color="black", alpha=0.5) 
      }
    }
    
    if(!is.null(add_point2)){
      if(nrow(add_point2)>0){
        p<-p +
          geom_point(data = add_point2, aes(x=mround(V2, bin_by), y=-(10.5/50)), color="violet", alpha=0.5) 
      }
    }
    
    if(!is.null(add_region)){
      if(nrow(add_region)>0){
 
        # If start and end of the resiong after rounding up ends up in the same bin (510-510) the geom_segment doesn't plot it
        # To avoid this we will 'cheat' and modify one of the boundaries (end) to fall to the next bin
        if(sum(mround(add_region$V2, bin_by) == mround(add_region$V3, bin_by))>0){
          add_region$V3[mround(add_region$V2, bin_by) == mround(add_region$V3, bin_by)]<-add_region$V3[mround(add_region$V2, bin_by) == mround(add_region$V3, bin_by)]+5
        }
        
        p<-p +
          geom_segment(data = add_region, aes(x=mround(V2, bin_by), xend=mround(V3, bin_by), y=-(10.5/15), yend=-(10.5/15)), size=2, color="blue", alpha=0.5)
      }
    }
    
    print(p)
    
  dev.off()
}  