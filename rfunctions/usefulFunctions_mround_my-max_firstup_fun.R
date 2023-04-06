#Make nicer numbers rounded up, useful for making scales for a plot.
mround <- function(x,base){ 
  base*ceiling(x/base) 
} 

#max robust to NAs
my_max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA) 

#Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
