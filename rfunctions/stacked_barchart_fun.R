stacked_barchart <- function(melted_tpm,
                             threshold = 1,
                             stack_group_name = 'Geneid',
                             sd_bars = FALSE,
                             inset = TRUE){
  #Have to use sym and !! in order to get mutate to work 
  # (converting the stack_group_name to a symbol and then unquoting in the mutate function)
  stack_group_name <- rlang::sym(paste(stack_group_name))
  
  tpm_data <- mutate(melted_tpm, 
                     stack_group = (!!stack_group_name)) %>% 
    mutate(tpm = ifelse(tpm > threshold, tpm, 0)) %>% #Filtering step for low value pir genes.
    group_by(stage, SorL, stack_group, sample) %>% 
    summarise(sum_tpm = sum(tpm)) %>% 
    group_by(SorL, stack_group, stage) %>%
    summarise(mean_sum_tpm = mean(sum_tpm),
              sd_sum_tpm = sd(sum_tpm)) %>% 
    arrange(desc(SorL),
            mean_sum_tpm) %>% 
    ungroup() %>% 
    mutate(order = row_number())
  
  tpm_data_SorL_stackgenes <- 
    #Order by SorL and then by mean_sum_tpm, then convert stack_group to a factor so ggplot orders correctly
    arrange(tpm_data,
            stage,
            desc(SorL), 
            mean_sum_tpm) %>% 
    mutate(stage = factor(stage,
                          levels = ordering(unique(tpm_data$stage),
                                            order_vector = r_order_vector)),
           stage_num = as.integer(
             factor(stage,
                    levels = ordering(unique(tpm_data$stage),
                                      order_vector = r_order_vector))
           )) %>% 
    group_by(stage) %>% 
    mutate(stacked_mean = cumsum(mean_sum_tpm)) %>% 
    mutate(sd_upper = stacked_mean + sd_sum_tpm,
           sd_lower = ifelse(sd_sum_tpm > mean_sum_tpm,
                             lag(stacked_mean, default = 0), 
                             stacked_mean - sd_sum_tpm)
    )
  
  tpm_data_SorL_stackgenes$stack_group_order <- seq_len(length.out = nrow(tpm_data_SorL_stackgenes))
  
  stacked_max <- tpm_data %>% 
    ungroup %>% 
    group_by(stage) %>% 
    dplyr::summarise(max_tpm = sum(mean_sum_tpm)) %>% 
    .$max_tpm %>% 
    max
  
  max_tpm <- c(
    seq(
      0,
      mround(
        stacked_max,
        500),
      length.out = 11)
  )
  
  stack_label <- 'sub-family'
  
  stage_vec <- tpm_data_SorL_stackgenes$stage
  col_group <- case_when(str_detect(stage_vec, pattern = 'Asex') ~ 'zzAsex' ,
                         str_detect(stage_vec, pattern = 'Liv') ~ 'zzLiv',
                         str_detect(stage_vec, pattern = 'Gam') ~ 'zzGam',
                         str_detect(stage_vec, pattern  = 'Bld|Sporo|Ook') ~ 'zzMos')
  tpm_data_SorL_stackgenes$stage_col <- col_group
  fills <- c('#F6B9AD','#96d3ff','#f0d4e3','#a7e9af') 
  #Note that this order is influenced by sort(col_group)
  
  # Note that this code does not work with ggplot 3.3.5, but it does work with 3.3.3. 
  #  Some issue with the colours for the geom_rect background.
  
  y <- ggplot(data = tpm_data_SorL_stackgenes, 
              aes(
                x = stage_num,
                y = mean_sum_tpm, 
                fill = SorL
              )) +
    geom_rect(
      aes(
        xmin = stage_num - .7,
        xmax = stage_num + .7,
        ymin = -Inf,
        ymax = Inf,
        fill = stage_col
      ),
      alpha = .5
    ) +
    geom_bar(position = position_stack(),
             stat = 'identity',
             colour = 'black',
             size = 0.1) +
    {if(sd_bars)list(
      geom_errorbar(aes(ymax = sd_upper,
                        ymin = sd_lower),
                    position = position_dodge(width = .85),
                    size = .25),
      geom_point(aes(y = stacked_mean), 
                 position = position_dodge(width = .85),
                 show.legend = FALSE,
                 size = .25)
    )}+ 
    scale_fill_manual(
      breaks = unique(tpm_data_SorL_stackgenes$SorL),
      values = c(
        SorL_col_num,
        fills
      ),
      name = 'Sub-family'
    ) +
    scale_color_manual('black') + 
    scale_y_continuous(breaks = max_tpm) +
    theme_classic() +
    ylab(
      bquote('total '~italic(.('pir'))~' TPM')
    ) +
    scale_x_continuous(
      breaks = tpm_data_SorL_stackgenes$stage_num %>% unique(),
      labels = tpm_data_SorL_stackgenes$stage %>% unique(),
      expand = c(0, 0)
    ) +
    theme(
      axis.text.x = element_text(angle = 90, 
                                 vjust = .5,
                                 size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.key.size = unit(3,'mm'),
      legend.key.width = unit(2,'mm'),
    ) +
    xlab('stage') + 
    facet_grid(.~stage, 
               scales = "free_x",
               space = "free_x") +
    scale_x_discrete(expand=c(0,0)) +
    theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.key.size = unit(3,'mm'),
      legend.key.width = unit(2,'mm'),
      strip.text.x = element_text(size=8, 
                                  face = 'bold', 
                                  angle = 90),
      strip.background = element_rect(
        color="black", 
        size=0.1, 
        linetype="solid"
      ),
      panel.spacing.x = unit(0, "null")
    ) 
  
  g <- ggplot_gtable(ggplot_build(y))
  strip_t <- which(grepl('strip-t', g$layout$name))
  stage_vec_strip <- sort(unique(tpm_data_SorL_stackgenes$stage))
  fills <- case_when(str_detect(stage_vec_strip, pattern = 'Asex') ~ 'coral1' ,
                     str_detect(stage_vec_strip, pattern = 'Liv') ~ 'plum2',
                     str_detect(stage_vec_strip, pattern = 'Gam') ~ 'lightskyblue1',
                     str_detect(stage_vec_strip, pattern  = 'Bld|Sporo|Ook') ~ 'lightgreen')
  k <- 1
  for (i in strip_t) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  grid.newpage()
  grid.draw(g)
  
  #Draw inset of the smaller value samples
  if(inset){
    #Include all stages expressed lower than a certain amount
    
    low_express <- tpm_data_SorL_stackgenes %>% 
      group_by(stage) %>% 
      summarise(sum_tpm = sum(mean_sum_tpm),
                max_stan_dev = max(sd_upper)) %>% 
      filter(sum_tpm < 350)
    
    #Calculate the range for the y value to make it look nicer
    
    stacked_max <- max(low_express$max_stan_dev)
    
    max_tpm <- c(
      seq(
        0,
        mround(
          stacked_max,
          100),
        length.out = 11)
    )
    
    tpm_data_SorL_stackgenes_low <- filter(tpm_data_SorL_stackgenes,
                                           stage %in% low_express$stage)
    
    y <- ggplot(data = tpm_data_SorL_stackgenes_low, 
                aes(
                  x = stage_num,
                  y = mean_sum_tpm, 
                  fill = SorL
                )) +
      geom_bar(position = position_stack(),
               stat = 'identity',
               colour = 'black',
               size = 0.1) +
      {if(sd_bars)list(
        geom_errorbar(aes(ymax = sd_upper,
                          ymin = sd_lower),
                      position = position_dodge(width = .85),
                      size = .25),
        geom_point(aes(y = stacked_mean), 
                   position = position_dodge(width = .85),
                   show.legend = FALSE,
                   size = .25)
      )}+ 
      scale_fill_manual(
        breaks = unique(tpm_data_SorL_stackgenes$SorL),
        values = c(
          SorL_col_num,
          fills
        )
      ) +
      scale_color_manual('black') + 
      scale_y_continuous(breaks = max_tpm) +
      theme_classic() +
      ylab(
        bquote('total '~italic(.('pir'))~' TPM')
      ) +
      scale_x_continuous(
        breaks = tpm_data_SorL_stackgenes$stage_num %>% unique(),
        labels = tpm_data_SorL_stackgenes$stage %>% unique(),
        expand = c(0, 0)
      ) +
      theme(
        axis.text.x = element_text(angle = 90, 
                                   vjust = .5,
                                   size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(3,'mm'),
        legend.key.width = unit(2,'mm'),
      ) +
      xlab('stage') + 
      facet_grid(.~stage, 
                 scales = "free_x",
                 space = "free_x") +
      scale_x_discrete(expand=expand_scale(add=1)) +
      theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(3,'mm'),
        legend.key.width = unit(2,'mm'),
        strip.text.x = element_text(size=8, 
                                    face = 'bold', 
                                    angle = 90),
        strip.background = element_rect(
          color="black", 
          size=0.1, 
          linetype="solid"
        ),
        panel.spacing.x = unit(0, "null")
      ) 
    
    #This bit is to change the colour of the facet boxes
    
    g_inset <- ggplot_gtable(ggplot_build(y))
    strip_t <- which(grepl('strip-t', g_inset$layout$name))
    stage_vec <- unique(
      levels(tpm_data_SorL_stackgenes_low$stage)[levels(tpm_data_SorL_stackgenes_low$stage) %in% as.character(tpm_data_SorL_stackgenes_low$stage)]
    )
    fills <- case_when(str_detect(stage_vec, pattern = 'Asex') ~ 'coral1' ,
                       str_detect(stage_vec, pattern = 'Liv') ~ 'plum2',
                       str_detect(stage_vec, pattern = 'Gam') ~ 'lightskyblue1',
                       str_detect(stage_vec, pattern  = 'Bld|Sporo|Ook') ~ 'lightgreen')
    k <- 1
    for (i in strip_t) {
      j <- which(grepl('rect', g_inset$grobs[[i]]$grobs[[1]]$childrenOrder))
      g_inset$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    grid.newpage()
    grid.draw(g_inset)
    return(list("main" = g,
                "inset" = g_inset,
                "data" = tpm_data_SorL_stackgenes))
  } else {
    return(list("main" = g,
                "data" = tpm_data_SorL_stackgenes))
  }
  
}