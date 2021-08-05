#parsing bam file into 45s units

library(Rsamtools)
library(tidyverse)
library(ggplot2)
source("GC_p_meth_45s.R")

#read in the data
#create pointer to bam file
bamPath <- "path/to/all.bam"
bamFile <- BamFile(bamPath)
#read in bam file with name, sequence, and methylC tags
bam <- scanBam(bamFile, param=ScanBamParam(what=c("qname", "seq"), tag="MC"))

#read in 45s coordinates
coords <- read_tsv("path/to/all_merged-45s.paf", col_names = FALSE) %>% 
  select(1,3:5) %>% `colnames<-`(c("ID", "start", "end", "dir"))
#colnames must be ID, start, end
#coltypes must be char, int, int

#process each 45s region
df_sum <- data.frame(t(sapply(1:nrow(coords), GC_p_meth_45s, bam=bam, coords=coords))) #%>% 
df <- df_sum %>% mutate(mean = as.numeric(mean), 
                        p_CG = as.numeric(p_CG),
                        read_ID = as.character(read_ID),
                        length=as.numeric(length),
                        sd = as.numeric(sd),
                        unit_45s_ID = as.character(unit_45s_ID))

#read in reference file for which read belongs to which region
read_ref <- read_tsv(file = "/path/to/read_ref.txt", col_names = FALSE) %>% 
            `colnames<-`(c("ID", "chr", "reg")) %>% unique()
  
df_merged <- df %>% left_join(read_ref, by = c("read_ID" = "ID")) %>% mutate(chr = as.factor(chr))
#write out this dataframe to use later 
write.table(df_merged, file = "/path/to/read_summaries.txt", sep = '\t', col.names = TRUE, quote = FALSE, row.names = FALSE)

df %>% ggplot(aes(x = chr, y = mean, fill = chr)) + 
  geom_violin() + geom_boxplot(width=0.1) + 
  theme(legend.position = "none") + 
  xlab("Acrocentric Chromosome") + ylab("Mean CpG site methylation \nprobability for 45s units") + ggtitle("Mean CpG Site Methylation Probability by Chromosome")

