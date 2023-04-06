draw_heatmap_zscore <- function(zscore_list,
                                title = '',
                                gene_ordering_vector = NULL,
                                filter_genes_vector = NULL,
                                extra_title = NULL,
                                annotation_clade = NULL,
                                annotation_col,
                                threshold_tpm = NULL,
                                use_times = FALSE,
                                use_order_vector = TRUE,
                                tpm_mat = TRUE,
                                tpm_breaks = NULL,
                                column_title_name = NULL,
                                ...){
  if(is.null(tpm_breaks) & tpm_mat){
    #Need the generate the breaks in the data for the heatmap legend if not provided.
    max_tpm <- lapply(zscore_list, 
                      function(list_element){
                        #list_element <- zscore_list[[1]]
                        heatmap_data <- as.data.frame(list_element$zscores)
                        title <- list_element$experiment
                        
                        if(!'Geneid' %in% colnames(heatmap_data)){
                          stop('No Geneid column in the heatmap data.frame')
                        }
                        
                        if(!is.null(filter_genes_vector)){
                          heatmap_data <- as.data.frame(
                            filter(heatmap_data, 
                                   Geneid %in% filter_genes_vector)
                          )
                        }
                        
                        rownames(heatmap_data) <- as.character(heatmap_data$Geneid)
                        
                        heatmap_data <- heatmap_data[, !colnames(heatmap_data) %in% 'Geneid', drop = TRUE]
                        #This won't work if heatmap_data is a tibble, must be a data.frame, so annoying!!
                        
                        if(!is.null(gene_ordering_vector)){
                          #Remove genes not in the main data frame but in the ordering vector, 
                          # and then order based on the vector
                          remove_index <- rownames(heatmap_data) %in% gene_ordering_vector
                          heatmap_data <- heatmap_data[remove_index,,drop = FALSE]
                          heatmap_data <- heatmap_data[na.omit(
                            match(gene_ordering_vector, 
                                  rownames(heatmap_data))
                          ),,
                          drop=FALSE]
                        }
                        
                        geneids <- rownames(heatmap_data)
                        
                        #Get max/min TPM scores for the side annotation heatmap                       
                        tpm_mat <- list_element$tpm[match(rownames(heatmap_data),
                                                          list_element$tpm$Geneid),] %>% 
                          select(-Geneid) %>% 
                          as.matrix %>% 
                          log1p
                        rownames(tpm_mat) <- geneids
                        
                        tpm_min <- rowMins(tpm_mat)
                        names(tpm_min) <- geneids
                        tpm_max <- rowMaxs(tpm_mat)
                        names(tpm_max) <- geneids
                        return(list(tpm_max = tpm_max,
                                    tpm_min = tpm_min))
                      })
    
    tpm_max_all <- mround(
      max(
        unlist(
          lapply(max_tpm, 
                 getElement,
                 'tpm_max')
        )
      ), 
      base = 1)
    
    tpm_breaks <- round(
      c(0,
        tpm_max_all/4,
        tpm_max_all/2,
        tpm_max_all*(3/4),
        tpm_max_all),
      digits = 2
    )
  }
  
  #Generate the colour function from the TPM legend breaks for the TPM annotation
  if(tpm_mat){
    tpm_col_fun <- colorRamp2(
      breaks = tpm_breaks,
      colors = viridis(length(tpm_breaks))
    )
  }
  
  if(is.numeric(threshold_tpm)){
    #Apply a threshold to the TPM in order to remove low expressed genes.
    # Remember to add log1p here! 
    # Threshold uses the max value, so it ony needs to be above thresh 
    #   in one stage avg in one transmission method.
    #Will alter the gene_filter vector in order to only select the genes above threshold
    
    gene_maxs <- bind_rows(
      lapply(max_tpm, 
             getElement,
             'tpm_max')
    )
    
    remove_genes <- names(gene_maxs)[apply(
      gene_maxs, 
      2, 
      function(maxs){
        all(maxs < log2(1+threshold_tpm))
      }
    )]
    
    #Set the filter_genes_vector as the genes to remove if it doesn't already exist.
    if(is.null(filter_genes_vector)){
      filter_genes_vector <- remove_genes
    }else{
      filter_genes_vector <- filter_genes_vector[!filter_genes_vector %in% remove_genes]
    }
  }
  
  if(is.character(annotation_clade)){
    order_cir <- gene_ordering_vector[gene_ordering_vector %in% filter_genes_vector]
    clade_info <- getElement(cir_info,
                             annotation_clade)[match(order_cir, 
                                                     cir_info$Geneid)]
    names(clade_info) <- str_split_fixed(order_cir, 
                                         pattern = '_', 
                                         n = 2)[,2]
    annotation_col <- annotation_col
  }
  
  heatmaps <- lapply(zscore_list, 
                     function(list_element){
                       #list_element <- zscore_list[[1]]
                       heatmap_data <- as.data.frame(list_element$zscores)
                       title <- list_element$experiment
                       
                       if(!'Geneid' %in% colnames(heatmap_data)){
                         stop('No Geneid column in the heatmap data.frame')
                       }
                       
                       if(!is.null(filter_genes_vector)){
                         heatmap_data <- as.data.frame(
                           filter(heatmap_data, 
                                  Geneid %in% filter_genes_vector)
                         )
                       }
                       
                       rownames(heatmap_data) <- as.character(heatmap_data$Geneid)
                       
                       heatmap_data <- heatmap_data[, !colnames(heatmap_data) %in% 'Geneid', drop = TRUE]
                       #This won't work if heatmap_data is a tibble, must be a data.frame, so annoying!!
                       
                       if(!is.null(gene_ordering_vector)){
                         #Remove genes not in the main data frame but in the ordering vector,
                         #  and then order based on the vector
                         remove_index <- rownames(heatmap_data) %in% gene_ordering_vector
                         heatmap_data <- heatmap_data[remove_index,,drop = FALSE]
                         heatmap_data <- heatmap_data[na.omit(
                           match(gene_ordering_vector, 
                                 rownames(heatmap_data))
                         ),,
                         drop=FALSE]
                       }
                       
                       if(tpm_mat){
                         #Get max/min TPM scores for the side annotation heatmap                       
                         tpm_mat.mat <- list_element$tpm[match(rownames(heatmap_data),
                                                               list_element$tpm$Geneid),] %>%
                           select(-Geneid) %>% 
                           as.matrix %>% 
                           log1p
                         
                         tpm_min <- rowMins(tpm_mat.mat)
                         tpm_max <- rowMaxs(tpm_mat.mat)
                         
                         names(tpm_min) = rownames(heatmap_data)
                         names(tpm_max) = rownames(heatmap_data)
                         
                         tpm_ha <- rowAnnotation(
                           min.max = anno_simple(cbind(tpm_min, 
                                                       tpm_max), 
                                                 col = tpm_col_fun,
                                                 # gp = gpar(col = 'black',
                                                 #           lwd = 0.1),
                                                 # border = TRUE,
                                                 simple_anno_size = unit(0.2, "cm")),
                           show_annotation_name = FALSE)
                       }
                       
                       #Draw the left annotation only for the first element of the heatmap list.
                       
                       if(
                         identical(list_element,
                                   dplyr::first(zscore_list)) & is.character(annotation_clade)
                       ){
                         clade_ha <- rowAnnotation(
                           'subfamily' = anno_simple(
                             clade_info,
                             col = annotation_col
                           )
                         )
                       }
                       
                       #Change rownames to remove PCHAS_
                       rownames(heatmap_data) = str_split_fixed(rownames(heatmap_data), 
                                                                pattern = '_', 
                                                                n = 2)[,2]
                       
                       if(use_times){
                         colnames(heatmap_data) <- str_extract(colnames(heatmap_data), 
                                                               pattern = '[[:digit:]]{2}h')
                       }
                       
                       if(use_order_vector){
                         #Order the columns by ordering_vector
                         heatmap_data <- heatmap_data[,order(match(colnames(heatmap_data), 
                                                                   order_vector)),
                                                      drop = FALSE]
                       }else{
                         heatmap_data <- heatmap_data[,order(colnames(heatmap_data)),
                                                      drop = FALSE] 
                       }
                       
                       h = Heatmap(heatmap_data,
                                   name = paste0('z-score',
                                                 extra_title),
                                   col = c(cbPalette[6],
                                           'white', 
                                           cbPalette[7]),
                                   row_names_gp = gpar(fontsize = ifelse(
                                     nrow(heatmap_data) > 13, 
                                     12*(13/nrow(heatmap_data)), 
                                     12)
                                   ),
                                   column_names_gp = gpar(fontsize = 10),
                                   row_dend_width = unit(3, "cm"),
                                   right_annotation = if(tpm_mat){
                                     tpm_ha
                                   }else{
                                     NULL
                                   },
                                   left_annotation = 
                                     if(identical(list_element,
                                                  dplyr::first(zscore_list)) & is.character(annotation_clade)){
                                       clade_ha
                                     }else{
                                       NULL
                                     }, 
                                   #Only want the clade annotation when it is first heatmap
                                   width = unit(4.5, 'cm'),
                                   cluster_columns = FALSE,
                                   height = unit(7,'cm'),
                                   row_title = NULL,
                                   cluster_rows = FALSE,
                                   column_title = paste(title),
                                   heatmap_legend_param = list(
                                     at = rev(c(3,2,1,0,-1,-2,-3)),
                                     labels_gp = gpar(fontsize = 10),
                                     title_gp = gpar(fontsize = 10)
                                   ),
                                   ...
                       )
                       return(h)
                     })
  
  ht_list <- Reduce('+', heatmaps) #To add the list elements together
  
  if(tpm_mat){
    #Draw the legend
    tpm_lgd <- Legend(
      title = "Min/Max log(TPM + 1)", 
      col_fun = tpm_col_fun, 
      at = tpm_breaks, 
      labels = tpm_breaks,
      title_gp = gpar(fontsize = 7)
    )
  }
  
  if(is.character(annotation_clade)){
    clade_lgd <- Legend(
      title = 'sub-family',
      at = names(annotation_col),
      legend_gp = gpar(fill = annotation_col)
    )
  }
  
  draw(ht_list, 
       row_title_side = 'left',
       column_title_side = 'top',
       heatmap_legend_side = 'left',
       column_title = column_title_name,
       #Only include the annotation legends for the annotations we actually have.
       annotation_legend_list = if(is.character(annotation_clade) & tpm_mat){
         list(tpm_lgd,
              clade_lgd)
       }else if(tpm_mat & !is.character(annotation_clade) ){
         list(tpm_lgd)
       } else if(!tpm_mat & is.character(annotation_clade) ){
         list(clade_lgd)
       })
}