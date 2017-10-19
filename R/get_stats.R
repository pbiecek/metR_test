#' Summarize regions
#'
#' Summarize regions by minimum, maximmum, mean, standard deviation of methylation rate in two probes and methylation diff rate withib two probes and estimated quantile basen od methylation diff rate
#' @param data dataframe with specyfic columns: chr, poz, prob, no, meth, unmeth, meth.rate, tiles and possible tiles.common columns. This dataframe is result of function create.tiles.min.gap or
#' create.tiles.fixed.length
#' @return data.frame which is summing-up regions specifed by tiles and tiles.common columns in data
#' @export
get.stats <- function(data){

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
    group_by(chr, start, end, prob) %>%
    summarise(meth.mean = mean(meth.rate), meth.cov = n(), meth.sd = sd(meth.rate),
              meth.min = min(meth.rate), meth.max = max(meth.rate))  %>%
    gather(Var, val, starts_with("meth.")) %>%
    unite(Var1,Var, prob) %>%
    spread(Var1, val) %>% dplyr::select(-meth.cov_x) %>% rename(meth.cov = meth.cov_y) %>%
    mutate(meth.diff = abs(meth.mean_x - meth.mean_y)) -> data

  data$quantile <- mapply(get.quantile, x = data$meth.diff, n = data$meth.cov)
  data
  }




