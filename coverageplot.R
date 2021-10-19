#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-r", "--region"), type="character", default=NULL, 
              help="Region of the genome that you want to visualise. Can either be a region copied from Ensembl or the Ensembl gene id", metavar="character"),
  make_option(c("-b", "--bedbam"), type="character", default=NULL, 
              help="BED or BAM file with the correct extension", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="GTF file", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="Whole genome FASTA file", metavar="character"), 
  make_option(c("-o", "--oligo"), type="character", default=NULL, 
              help="If you want to visualise where an oligonucleotide will bind, use this argument. If multiple sequences, separate by ','.", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(tidyverse)
library(plotly)
library(ggpubr)

if (grepl(".bed", opt$bedbam)) {
bed <- read_table(opt$bedbam, col_names = FALSE, col_types = "cddccc") 
} else {
  system(paste0("bedtools bamtobed -i ", opt$bedbam, " > tmp.bed"))
  bed <- read_table("tmp.bed", col_names = FALSE, col_types = "cddccc")
  system("rm tmp.bed")
}
gtf <- read_table(opt$gtf, skip = 5, col_names = FALSE, col_types = "cccddccccccccccccc") 
gtf_clean <- gtf %>%
  filter(X3 == "gene") %>%
  mutate(gene = str_remove_all(X10, "[\";]")) 


ggcoverage <- function(region, OLIGO = NULL) {
  if (!(region %in% gtf_clean$gene) & grepl("ENSG", region)) {
    stop("proposed gene region not found")
  } else if (region %in% gtf_clean$gene) {
    chr <- gtf_clean %>% filter(gene == region) %>% .$X1
    maximum <- gtf_clean %>% filter(gene == region) %>% .$X5
    minimum <- gtf_clean %>% filter(gene == region) %>% .$X4
  } else {
    region <- str_remove_all(region, "[Chromosome ,]")
    chr <- str_remove(region, ":[0-9]*-[0-9]*")
    maximum <- as.double(str_remove(region, ".*:[0-9]*-"))
    minimum <- as.double(str_remove(str_remove(region, paste0(chr, ":")), paste0("-", maximum)))
  }
  seq <- unlist(strsplit(paste(system(paste0("samtools faidx ", opt$fasta, " ",
                             chr, ":", minimum, "-", maximum), intern = TRUE)[-1], collapse = ""), split = ""))
  genome_region <- tibble(chr = chr, 
                          sequence = seq, 
                          pos = minimum-1+1:length(seq))
  start = seq(minimum, maximum-2, by = 2)
  data <- bed %>%
    filter(X1 == chr, 
           X2 <= maximum, 
           X3 >= minimum) %>%
    rowwise() %>%
    mutate(position = list(start[(start >= X2) & (start <= X3)])) %>%
    unnest(cols = position) %>%
    group_by(position) %>%
    summarise(counts = n()) %>%
    ungroup() 
  
  coverage <- data %>%
    ggplot(aes(position, counts)) +
    geom_bar(stat = "identity", col = "black", fill = "black", width = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.ticks.y = element_line(),
          axis.text.x = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    xlab("") +
    ylab("") +
    ggtitle(paste0("\n", region)) +
    xlim(minimum, maximum)
  
  sequence_plot <- genome_region %>%
    ggplot(aes(pos, fill = sequence)) +
    geom_histogram(binwidth = 1) +
    theme_minimal() +
    theme(axis.text.y = element_blank(), 
          panel.grid = element_blank(), 
          legend.position = "none") +
    xlab("") +
    ylab("") +
    xlim(minimum, maximum) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#000000", 
                                 "#F0E442")) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#000000", 
                                 "#F0E442"))
  
  p1 <- ggplotly(sequence_plot, tooltip = c("sequence"), dynamicTicks = TRUE)
  p2 <- ggplotly(coverage, dynamicTicks = TRUE)
  
   if (!is.null(OLIGO)) {
     oligos <- tibble(sequence = c(OLIGO, toupper(spgs::reverseComplement(OLIGO))), 
                    direction = rep(c("original", "revcomp"), each = length(OLIGO)))
     OLIGO_plot <- oligos %>%
       rowwise() %>%
       mutate(start = list(str_locate(paste(seq, collapse =""), sequence)[,"start"]), 
              end = list(str_locate(paste(seq, collapse =""), sequence)[,"end"])) %>%
       unnest(cols = c(start, end)) %>%
       mutate(id = 1:nrow(.),
              start = start + minimum,
              end = end + minimum) %>%
       ggplot() +
       geom_segment(aes(x = start, xend = end, y= id, yend = id, col = direction, label = sequence), size = 2) +
       theme_minimal() +
       theme(axis.text.y = element_blank(), 
             panel.grid = element_blank(), 
             legend.position = "none") +
       xlab("") +
       ylab("") +
       xlim(minimum, maximum) +
       scale_color_manual(values = c("#E69F00", "#56B4E9"))
     p3 <- ggplotly(OLIGO_plot, dynamicTicks = TRUE, tooltip = c("sequence", "start", "end", "direction"))
     subplot(p2, p1, p3, nrows = 3, heights = c(0.90, 0.05, 0.05), shareX = TRUE)
   } else {
     subplot(p2, p1, nrows = 2, heights = c(0.95, 0.05), shareX = TRUE)
   }
}

htmlwidgets::saveWidget(ggcoverage(opt$region, unlist(strsplit(opt$oligo, split = ","))), "region_export.html")

