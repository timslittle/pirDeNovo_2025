read_HMMertxtgz <- function(hmmer_gz.path) {
  hmmer_gz_all.df <- readLines(hmmer_gz.path)
  header.index <- first(which(grepl(hmmer_gz_all.df, 
                                    pattern = 'E-value')))
  end.index <- first(which(grepl(hmmer_gz_all.df, 
                                 pattern = 'inclusion threshold|Domain annotation for each sequence')))
  if(any(grepl(hmmer_gz_all.df, pattern = 'Domain annotation for each sequence'))) end.index <- end.index - 2
  hmmer_gz.df <- hmmer_gz_all.df[c(header.index,
                                   (header.index+2):(end.index - 1))]
  hmmer_gz.df <- sapply(hmmer_gz.df, function(line){
    #Split by whitespace and remove the first element of each line (which is blank)
    str_split(line, pattern = "\\s+")[[1]][-1]
  })
  #Get the number of columns as this informs where the sequence information column is 
  # (this column can contain blank space and be difficult to parse otherwise)
  num_cols <- length(hmmer_gz.df[[1]])
  #Save column names for later
  col_nam <- hmmer_gz.df[[1]]
  #Edit column names to remove duplicates and clarify what they mean.
  col_nam <- c(paste0(col_nam[1:3], "_fullseq"), paste0(col_nam[4:6], "_bestdom"), paste0(col_nam[7:8], "_numdom"), col_nam[9:10])
  hmmer_gz.df <- lapply(hmmer_gz.df[-1], 
                        function(line){
                          data.frame(t(unlist(c(line[1:num_cols - 1], 
                                                paste(line[num_cols:length(line)], 
                                                      collapse = "")))))
                        }) %>% bind_rows
  colnames(hmmer_gz.df) <- col_nam
  hmmer_gz.df <-  mutate(hmmer_gz.df, Sequence = as.character(Sequence))
  
  #I can't get this working in 'mutate' for whatever reason, so have to do it manually and then add it in.
  hmm_cvrg.vec <- unlist(lapply(hmmer_gz.df$Sequence, function(Sequence) {
    as.numeric(strsplit(hmmer_gz_all.df[which(str_detect(hmmer_gz_all.df, 
                                                         pattern = Sequence))[2]+3],
                        split = '\\s+')[[1]][9]) - as.numeric(strsplit(hmmer_gz_all.df[which(str_detect(hmmer_gz_all.df, 
                                                                                                        pattern = Sequence))[2]+3],
                                                                       split = '\\s+')[[1]][8])
  }))
  
  hmmer_gz.df$hmm_cvrg <- hmm_cvrg.vec
  
  return(hmmer_gz.df)
}