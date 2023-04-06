
# We need to write a function to calculate the TPM from the count data we have generated. This function uses the _lengths.transcript_ object and should generate errors if it detects discrepancies between the Geneids so that the wrong transcript length is not used.

tpm = function(count_data,lengths) {
  
  #Print an error if the id columns do not match
  if(any(grepl(colnames(lengths), 
               pattern = 'Geneid', 
               ignore.case = T)) & any(grepl(colnames(count_data), 
                                                  pattern = 'Geneid', 
                                                  ignore.case = T))){
    if(all(count_data[grep(colnames(count_data), 
                           pattern = 'Geneid', ignore.case = T)] == lengths[grep(colnames(lengths), 
                                                                             pattern = 'Geneid', ignore.case = T)])) {
      print('Count data and length gene ids both match')
    } else {
      stop('Warning:Count data and length gene ids DO NOT MATCH')
    }
  } else {
    stop('Count data and/or length table supplied do not have a Geneid column, make sure they are ordered correctly')
  }
  
  #Remove the id column from supplied count_data
  # if(any(grepl(colnames(count_data), pattern = 'Geneid',ignore.case = T))){
    counts_tpm = count_data[,-grep(colnames(count_data), pattern = 'Geneid', ignore.case = T)]
  # } else{
  #   counts_tpm <- count_data}
  
  #Remove the id column from supplied lengths
  if(!any(grepl(colnames(lengths), pattern = 'Transcript.Length', ignore.case = F))){
    stop('Length table needs a column labelled \"Transcript.Length\" (with capitals)')
  }
    
    lengths.table = lengths[,'Transcript.Length', drop = TRUE]
    # }else{
  #   lengths.table <- lengths}
  
  #Now for the actual calculation
  rpm <- counts_tpm
  for (i in 1:dim(counts_tpm)[2]) {
    for (j in 1:dim(counts_tpm)[1]) {
      rpm[j,i] <- (counts_tpm[j,i]/(lengths.table[j]/1000))
    }
    scale.factor <- sum(rpm[,i])/(1E6)
    counts_tpm[,i] <- apply(rpm[,i,drop=F],1,function(x){x/scale.factor})
  }
  # if(any(grepl(colnames(count_data), pattern = 'Geneid', ignore.case = T))){
    counts_tpm = bind_cols(count_data[,grep(colnames(count_data), pattern = 'Geneid', ignore.case = T), drop = FALSE], 
                           counts_tpm)
  # } 
  return(counts_tpm)
}