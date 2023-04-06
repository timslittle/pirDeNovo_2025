#Use an order_vector of partial string matches to order.
ordering = function(to_order, order_vector){
  order = unlist(sapply(order_vector, 
                        function(x){unique(str_subset(to_order, 
                                                      pattern = paste(x)))}))
  order_missing = unique(to_order[!to_order %in% order])
  unique(c(order, order_missing), fromLast = TRUE)
}