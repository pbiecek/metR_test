#' Create regions for testing methylation difference
#'
#' The same regions are observation that maximum difference position is gaps.length argument
#' @param data dataframe with specyfic columns: chr, poz, prob, no, meth, unmeth, meth.rate. This dataframe is result of function preprocessing.
#' @param tiles.length integer number that specifes maximum difference between minimum and maximum position in the same methylation regions
#' @param common logi value. If TRUE this function creates second regions group that are min postion is (min postion + max position)/2 of k-region and
#' max position is (min position + max position) of k+1 region.
#' @return data.frame from parameter data with extra column tiles that is region id number within chromosomes and extra column
#' tiles.common if argument tiles.common is not null
#' @export
create.tiles.fixed.length <- function(data, tiles.length, common = NULL){

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
      !all(data$meth + data$unmeth == data$no) | !all(round(data$meth /data$no,3) == round(data$meth.rate,3)))
    stop("Error: Incorrect values in data")

  if(!zapsmall(tiles.length, 15) == round(tiles.length) | tiles.length <= 0)
    stop("Error: Incorrect tiles.length argument")

  data %<>% mutate(tiles = as.integer(poz %/% tiles.length))

  if (!is.null(common)){
    if(!is.logical(common))
      stop("Error: Incorrect common argument")

    if (!is.null(common) & common == T){

    common.length = tiles.length %/% 2
    data %<>% mutate(tiles.common = as.integer((poz + common.length) %/% tiles.length))
    }
  }

  data
}
