#finding the locations of the reads we want

#load packages
library(tidyverse)
library(stringr)

#read in reads
reads <- read_tsv("/path/to/read_IDs.txt", col_names = FALSE)

#read in read locations from nanopolish index
read_locs <- read_tsv("/path/to/combined_nodups.fastq.index.readdb", col_names = FALSE) #this might take awhile to run

reads_j <- reads %>% left_join(read_locs)
colnames(reads_j) <- c("ID", "PATH")

#write this out as an rdata file to use later if I need to 
saveRDS(reads_j2, file = "/path/to/reads_j.rds")
write.table(reads_j4, file = "/path/to/rDNA_read_paths.txt", sep = '\t', col.names = FALSE, quote = FALSE, row.names = FALSE)

#generate a list of fast5 file paths we will want to extract from
paths <- as.data.frame(table(reads_j$PATH)) %>% mutate(Var2 = paste0("/data/PoreTenders/nanopolish/", Var1)) %>% select(Var2)
write.table(paths, file = "/path/to/fast5_paths.txt", sep = '\t', col.names = FALSE, quote = FALSE, row.names = FALSE)