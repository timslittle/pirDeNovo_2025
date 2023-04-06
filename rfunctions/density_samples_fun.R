density_samples <- function(input, 
                            save = FALSE,
                            line_val = 1){
  #Get some colours for the lines
  colours = rep_len(cbPalette, length.out = length(colnames(input$tpm)[-1]))
  list_dens = lapply(colnames(input$tpm)[-1], function(sample){
    data = as.numeric(unlist(select(input$tpm, paste(sample))))
    dens = density(log1p(data), bw = bw_mean)
    #Use the index of the column to get the colour
    colouring = colours[grep(colnames(input$tpm)[-1], pattern = paste(sample))]
    list_return = list(dens, colouring)
    names(list_return) = c('dens','colouring')
    list_return
  })
  max_dens = max(unlist(lapply(list_dens, function(x){max(x$dens$y)})))
  max_x = max(unlist(lapply(list_dens, function(x){max(x$dens$x)})))
  line_width = 1
  par(new=FALSE)
  if(save){
    png(file = paste0('plots/thesis_chapter3_figure_',
                      input$experiment, 
                      '_densityplot_fixedbw.png'), 
        width = 15/2.54, 
        height = 11.5/2.54,
        units = "in",
        res = 330)
  }
  lapply(seq_along(list_dens), function(index){
    x = list_dens[[index]]
    #Change lty value every 9 samples (i.e. when the colours reset)
    # %/% integer division gives the integer result, e.g. 8%/%9 = 0, 9%/%9 = 1, 18%/%9 = 2.
    lty <- ((index-1L) %/% 9L) + 1L
    plot(x$dens, 
         xlim = c(0, max_x), 
         ylim = c(0,max_dens), 
         type = 'l', 
         col = x$colouring,
         main = paste(input$experiment),
         lwd = line_width,
         lty = lty,
         xlab = '')
    if(is.numeric(line_val)){
      abline(v = log1p(line_val))
    }
    legend(max_x-4, 
           max_dens, 
           legend = paste(colnames(input$tpm)[-1]), 
           col = colours,
           lty = c(rep(1, times = 9), 
                   rep(2, times = 9), 
                   rep(3, times = 9)),
           cex = 0.5,
           ncol = ifelse(length(list_dens) > 15, 2, 1),
           xjust = ifelse(length(list_dens) > 15, 0.25, 0))
    par(new=TRUE)})
  if(save){
    dev.off()
  }
}