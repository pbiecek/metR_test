
load('data/quantile.function.rda')

get.quantile <- function(x, n){
  quantile.function[[n]](x)
}
