test.ks <- function(data){
  data %<>% arrange(prob, poz) 
  n = nrow(data)/2
  tryCatch({p.value = ks.test(x = data$meth.rate[1:n], y = data$meth.rate[(n+1):(2*n)],
                              alternative ="two.sided")$p.value
  return(data.frame(p.value = p.value))}            
  ,error=function(cond) {return(data.frame(p.value = NA))})
}      