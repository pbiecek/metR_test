#' Finding DMR
#'
#' Finding DMR by Wilcoxon, t-Student, Kolmogorov-Smirnow tests or logistic regression, logistic regression with mixed models,
#' logistic regression with mixed models with correlation matrix
#' @param data There are two options:  1. dataframe with specyfic columns: chr, poz, prob, no, meth, unmeth, meth.rate.
#' This dataframe is result of function preprocessing.
#' 2. dataframe with specyfic columns: chr, poz, prob, no, meth, unmeth, meth.rate, tiles and possible tiles.common columns. This dataframe is result of function create.tiles.min.gap or
#' create.tiles.fixed.length.
#' @param methods vectors with given methods. Possible values are: 'Wilcoxon', 'Ttest', 'KS', 'Reg.Log', 'Reg.Mixed',
#' 'Reg.Corr.Mixed'.
#' @param p.value.log.reg if not NULL regions with p.value of prob variable smaller than p.value.log.reg are returned and  decreasingly ordered by absolute value of beta coefficient
#' of prob variable otherwise regions ale increasingly ordered by p.value
#' @param p.value.reg.mixed if not NULL regions with p.value of prob variable smaller than p.value.log.reg are returned and  decreasingly ordered by absolute value of beta coefficient
#' of prob variable otherwise regions ale increasingly ordered by p.value
#' @param p.value.reg.corr.mixed if not NULL regions with p.value of prob variable smaller than p.value.log.reg are returned and  decreasingly ordered by absolute value of beta coefficient
#' of prob variable otherwise regions ale increasingly ordered by p.value
#' @return list object. Elements of list are results of given methods. The most interesting regions are on the top
#' @export
find.DMR <- function(data, methods, p.value.log.reg = NULL, p.value.reg.mixed= NULL, p.value.reg.corr.mixed= NULL){


  require(dplyr)
  require(data.table)
  require(magrittr)

  if ('tiles.common' %in% colnames(data)) {
    data.colnames <- c('chr', 'poz', 'prob', 'no', 'meth', 'unmeth', 'meth.rate', 'tiles', 'tiles.common')

    if(!all.equal(colnames(data), data.colnames))
      stop("Error: Incorrect colnames in data")

    data.types <- c("character","integer","character","integer","integer", "integer", 'double', 'integer', 'integer')
    names(data.types) <- data.colnames

    if(!all.equal(sapply(data, typeof), data.types))
      stop("Error: Incorrect datatypes in data")

    if (!all(data$poz > 0) | !all(data$chr %in% paste0('chr', c(1:22, "X", "Y", "M"))) |
        !all(data$meth + data$unmeth == data$no) | !all(round(data$meth /data$no,3) == round(data$meth.rate,3)) |
        !all(data$tiles >= 0) | !all(zapsmall(data$tiles, 15) == round(data$tiles)) |
        !all(data$tiles.common >= 0) | !all(zapsmall(data$tiles.common, 15) == round(data$tiles.common)))
      stop("Error: Incorrect values in data")

  }else{

    data.colnames <- c('chr', 'poz', 'prob', 'no', 'meth', 'unmeth', 'meth.rate', 'tiles')

    if(!all.equal(colnames(data), data.colnames))
      stop("Error: Incorrect colnames in data")

    data.types <- c("character","integer","character","integer","integer", "integer", 'double', 'integer')
    names(data.types) <- data.colnames

    if(!all.equal(sapply(data, typeof), data.types))
      stop("Error: Incorrect datatypes in data")

    if (!all(data$poz > 0) | !all(data$chr %in% paste0('chr', c(1:22, "X", "Y", "M"))) |
        !all(data$meth + data$unmeth == data$no) | !all(round(data$meth /data$no,3) == round(data$meth.rate,3)) |
        !all(data$tiles >= 0) | !all(zapsmall(data$tiles, 15) == round(data$tiles)))
      stop("Error: Incorrect values in data")

  }

result <- list()

data %>% group_by(chr, tiles) %>%
  summarize(start = min(poz), end = max(poz)) %>%
  dplyr::select(c(chr, start, end)) -> map

if ('tiles.common' %in% colnames(data)) {
  data %>% group_by(chr, tiles.common) %>%
  summarize(start = min(poz), end = max(poz)) %>%
    dplyr::select(c(chr, start, end)) -> map.common

rbind(map, map.common) %>% distinct() -> map
rm(map.common)
}

data %<>% mutate(poz.new = poz)

setDT(data)
setDT(map)

data <- data[map, on=.(poz.new >= start, poz.new <= end, chr = chr), nomatch = 0,
             .(chr, poz, prob, no, meth, unmeth, meth.rate, start, end), allow.cartesian = T] %>%
  group_by(chr, start, end)



if ('Wilcoxon' %in% methods){
  print("Started: Finding DMR by Wilcoxon test")
  result$Wilcoxon <- data %>% do(test.wilcoxon(data=.)) %>%
    filter(!is.nan(p.value)) %>% arrange(p.value) %>% ungroup()
}

if ('Ttest' %in% methods){
  print("Started: Finding DMR by t-test")
  result$Ttest <- data %>% do(test.t(data=.)) %>%
    filter(!is.nan(p.value)) %>% arrange(p.value) %>% ungroup()
}

if ('KS' %in% methods){
  print("Started: Finding DMR by KS test")
  result$KS <- data %>% do(test.ks(data=.)) %>%
    filter(!is.nan(p.value)) %>% arrange(p.value) %>% ungroup()
}


if ('Reg.Log' %in% methods){
  print("Started: Finding DMR by Logistic Regression")

  if (!is.null(p.value.log.reg)){
    result$Reg.Log <- data %>% do(reg.log(data=.)) %>%
      filter(p.value < p.value.log.reg) %>% arrange(-abs(beta.coef)) %>% ungroup()
  }else{
    result$Reg.Log <- data %>% do(reg.log(data=.)) %>%
      filter(!is.nan(p.value)) %>% arrange(p.value) %>% ungroup()}

  }


if ('Reg.Mixed' %in% methods){
  print("Started: Finding DMR by Logistic Regression with Mixed Effects")

  if (!is.null(p.value.reg.mixed)){
    result$Reg.Mixed <- data %>% do(reg.mixed(data=.)) %>%
      filter(p.value < p.value.reg.mixed) %>% arrange(-abs(beta.coef)) %>% ungroup()
  }else{
    result$Reg.Mixed <- data %>% do(reg.mixed(data=.)) %>%
      filter(!is.nan(p.value)) %>% arrange(p.value) %>% ungroup()}

}


if ('Reg.Corr.Mixed' %in% methods){
  print("Started: Finding DMR by Logistic Regression with Mixed Effects with Correlation Matrix")
  load('mean.acf.chr.rda')
  acf <- mean.acf.chr[-1]
  if (!is.null(p.value.reg.corr.mixed)){
    result$Reg.Corr.Mixed <- data %>% do(reg.corr.mixed(data=., acf = acf)) %>%
      filter(p.value < p.value.reg.corr.mixed) %>% arrange(-abs(beta.coef)) %>% ungroup()
  }else{
    result$Reg.Corr.Mixed <- data %>% do(reg.corr.mixed(data=., acf = acf)) %>%
      filter(!is.nan(p.value)) %>% arrange(p.value) %>% ungroup()}

}
return(result)
}

