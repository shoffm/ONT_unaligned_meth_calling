GC_p_meth_45s <- function(index, bam, coords){
  #find read name
  ID <- coords$ID[index]
  #find start position
  start <- coords$start[index] #+ 1 #to correct for zero indexing
  #find end position
  end <- coords$end[index] #+ 1 #to correct for zero indexing
  
  #reach into bam file and find which index the read name has
  read_index = try({
    which(bam[[1]]$qname == ID)
  })
  if(length(read_index) == 0){
    print(ID)
    print(index)
    print("DNE")
    returns <- list(mean = NA, 
                    p_CG = NA,
                    length = NA,
                    sd = NA,
                    read_ID=ID,
                    unit_45s_ID=paste(ID,index, sep="."))
    return(returns)
  }
  
  #extract that read
  str <- strsplit(as.character((bam[[1]]$seq[read_index])), "")[[1]]
  #get the list of the tags
  tag <- bam[[1]]$tag$MC[[read_index]]
  #make a fake position matrix
  pos <- 0:length(str)
  
  #make df
  df <- data.frame(cbind(str,tag,pos)) %>% 
    mutate(m5c_tag = ifelse(str == "C", tag, NA), 
           m5c_tag=as.character(m5c_tag),
           str=as.character(str),
           tag=as.numeric(tag), 
           pos=as.numeric(pos))
  
  #subset to the region of interest
  test_df <- df %>% filter(pos %in% (start:end)) %>%
    mutate(pos_45s = 0:(end-start))
  #test_df <- test_df %>% mutate(pos_45s = 1:nrow(test_df))
  #add a 45s position
  
  #look for CGs next to each other, filter to these 
  CGs <- str_locate_all(as.character(paste(test_df$str, collapse="")), "CG")[[1]]
  #make this into one list
  CGs_list <- c(CGs[,1], CGs[,2])
  #order the list correctly
  CGs_list_ord <- sort(CGs_list)
  
  
  
  #test_df_filtered <- data.frame(test_df[pos %in% CGs_list_ord,])
  test_df_filtered <-test_df %>% filter(pos_45s %in% CGs_list_ord, str == "C")
  print(bam[[1]]$qname[read_index])
  print(index)
  returns <- list(mean = mean(as.numeric(test_df_filtered$tag))/255, 
                  p_CG = (nrow(test_df_filtered)*2)/nrow(test_df),
                  length = nrow(test_df),
                  sd = sd(as.numeric(test_df_filtered$tag)/255),
                  read_ID=bam[[1]]$qname[read_index],
                  unit_45s_ID=paste(bam[[1]]$qname[read_index],index, sep="."))
  return(returns)
  
}
