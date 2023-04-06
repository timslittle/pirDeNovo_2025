#Make unique identifiers from a vector of names by appending '_X' where X is the number of the 'replicate'.
uniquify <- function(vector_names, first_one = TRUE){
  new_vector <- c()
  df_num <- data.frame(names = unique(vector_names), rep = 1)
  if(first_one){
    duplicated_names <- unique(vector_names[duplicated(vector_names)])
  }else{
    duplicated_names <- c()
  }
  for(i in 1:length(vector_names)){
    if(vector_names[i] %in% duplicated_names|!first_one){
      name <- paste(vector_names[i],
                    df_num$rep[match(vector_names[i], df_num$names)],
                    sep = '_')
      #Increase rep counter
      df_num$rep[match(vector_names[i], df_num$names)] <- df_num$rep[match(vector_names[i], df_num$names)] + 1
      new_vector[i] <- name
    }else{
      new_vector[i] <- vector_names[i]
    }
  }
  new_vector
}