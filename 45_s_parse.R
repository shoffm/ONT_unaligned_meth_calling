#parsing bam file into 45s units

library(Rsamtools)
library(tidyverse)
library(slider)
library(ggplot2)
library(grid)
library(gridExtra)
source("GC_p_meth_45s.R")

#read in the data
#create pointer to bam file
bamPath <- "path/to/all.bam"
bamFile <- BamFile(bamPath)
#read in bam file with name, sequence, and methylC tags
bam <- scanBam(bamFile, param=ScanBamParam(what=c("qname", "seq"), tag="MC"))

#read in 45s coordinates
coords <- read_tsv("....")
#colnames must be ID, start, end
#coltypes must be char, int, int

#process each 45s region
df <- data.frame(t(sapply(1:nrow(coords), GC_p_meth_45s, bam=bam, coords=coords))) %>% 
  mutate(mean = as.numeric(mean), p_CG = as.numeric(p_CG), 
         read_ID = as.character(samp), length=as.numeric(length), 
         sd = as.numeric(sd), unit_45s_ID = as.character(unit_45s_ID))

#read in reference file for which read belongs to which region
read_strats <- read_tsv("/path/to/read_ref.txt", col_names = FALSE) %>% 
  setNames(c("ID", "chr", "reg")) %>% mutate(chr = as.factor(chr), reg = as.factor(reg)) %>% distinct()

df_merged <- df %>% left_join(read_strats, by = c("read_ID" = "ID")) %>% mutate(chr = as.factor(chr))
#write out this dataframe to use later 
write.table(df_merged, file = "/path/to/read_summaries.txt", sep = '\t', col.names = TRUE, quote = FALSE, row.names = FALSE)
