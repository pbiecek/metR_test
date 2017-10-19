reg.mixed <- function(data){
  require(lme4)
  tryCatch({
    mod.temp <- glmer(cbind(meth, unmeth) ~ prob + (1|poz), 
                      data = data, family = binomial)
    p.value <- summary(mod.temp)$coefficients[2,4]
    beta.coef<- summary(mod.temp)$coefficients[2,1]
    return(data.frame(p.value = p.value, beta.coef = beta.coef))
  }, error = function(cond){
    return(data.frame(p.value = NA, beta.coef = NA))})
}