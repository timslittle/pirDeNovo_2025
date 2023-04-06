ecdf_percentile <- function(df_tpm){
  Fn <- ecdf(
    unlist( #To convert into an atomic vector that ecdf will like.
      select_at(
        df_tpm,
        vars(
          matches('med_fpkm|med_tpm|^tpm$') #Default matches behaviour is ignore.case = TRUE
        ) #Could be med_fpkm or med_tpm.
      )
    )
  )
  #Using a lapply to make a list in case there are multiple ancestral genes.
  list_percentile <- lapply(
    df_tpm$Geneid[
      df_tpm$ancestral_ortholog == 1
    ],
    function(ancestor_gene){
      #For debugging
      # ancestor_gene <- "PCHCB_0101200"
      list(
        Fn(
          filter(df_tpm, Geneid %in% ancestor_gene) %>% 
            select_at(
              vars(
                matches('med_fpkm|med_tpm|^tpm$')
              ) #Could be med_fpkm or med_tpm
            ) 
        )
      )
    }
  )
  #Create data.frame of the output
  gene_percentile_df <- data.frame(
    Geneid = df_tpm$Geneid[
      df_tpm$ancestral_ortholog ==1
    ],
    percentile = unlist(
      lapply(
        list_percentile, `[`, 1
      )
    )
  )
  plot(Fn)
  return(gene_percentile_df)
}