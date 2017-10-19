reg.corr.mixed <- function(data, acf){

  require(nlme)
  require(MASS)

  data %<>% arrange(prob, poz)

  tryCatch({
    n <- nrow(data)/2
    poz <- data$poz[1:n]

    m1 <- matrix(poz, nrow = length(poz), ncol = length(poz))
    m2 <- t(m1)
    m3 <- m1 - m2
    m3 <- m3[m3>0]
    cor <- acf[m3]

    data$poz2 <- c(1:n, 1:n)
    X <- diag(1, nrow = n)
    X[lower.tri(X)] <- cor
    X[upper.tri(X)] <- cor
    mat_corr <- nearPD(X, corr= TRUE)$mat
    cs <- Initialize(corSymm(value =  mat_corr[lower.tri(mat_corr)], form= ~ poz2, fixed = TRUE), data=data[1:n,])

    mod.temp <- glmmPQL(cbind(meth, unmeth) ~ prob , random = ~1|factor(poz),
                        data = data, family = binomial, correlation = cs, verbose = FALSE)
    p.value <- summary(mod.temp)$tTable[2,1]
    beta.coef <- summary(mod.temp)$tTable[2,1]
    return(data.frame(p.value = p.value, beta.coef = beta.coef))
  }, error = function(cond){
    return(data.frame(p.value = NA, beta.coef = NA))})
}
