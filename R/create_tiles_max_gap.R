#' Create regions for testing methylation difference
#'
#' The same regions are observation that maximum difference position is gaps.length argument
#' @param data dataframe with specyfic columns: chr, poz, prob, no, meth, unmeth, meth.rate. This dataframe is result of function preprocessing.
#' @param gaps.length integer number that specifes maximum difference position between the same methylation regions
#' @return data.frame from parameter data with extra column tiles that is region id number within chromosomes
#' @export
create.tiles.max.gap <- function(data, gaps.length){


  require(dplyr)
  require(data.table)

  data.colnames <- c('chr', 'poz', 'prob', 'no', 'meth', 'unmeth', 'meth.rate')

  if(!all.equal(colnames(data), data.colnames))
    stop("Error: Incorrect colnames in data")

  data.types <- c("character","integer","character","integer","integer", "integer", 'double')
  names(data.types) <- data.colnames

  if(!all.equal(sapply(data, typeof), data.types))
    stop("Error: Incorrect datatypes in data")

  if (!all(data$poz > 0) | !all(data$chr %in% paste0('chr', c(1:22, "X", "Y", "M"))) |
      !all(data$meth + data$unmeth == data$no) | !all(round(data$meth /data$no,5) == round(data$meth.rate,5)))
    stop("Error: Incorrect values in data")

  if(!zapsmall(gaps.length, 15) == round(gaps.length) | gaps.length <= 0)
    stop("Error: Incorrect gaps.length argument")

  data <- split(data, f = data$chr)
  data <- lapply(data, create.tiles.chr, gaps.length)
  do.call(rbind, data)
}



