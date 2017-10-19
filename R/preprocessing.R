#' Preprocessing data
#'
#' Preprocessing data for methylation analyses
#' @param path.1 a path to file.1 with methylation data
#' @param path.2 a path to file.2 with methylation data
#' @return data.frame with methylation values on common positions and chromosomes from data in path.1 and path.2
#' @export
preprocessing <- function(path.1, path.2){

  require(data.table)
  require(dplyr)
  require(readr)
  require(magrittr)
  require(tidyr)


  if(!file.exists(path.1))
    stop("Error: Path to sample.1 is incorrect")

  if(!file.exists(path.2))
    stop("Error: Path to sample.2 is incorrect")

  sample.1 <- suppressMessages(read_delim(path.1, delim = "\t"))
  sample.2 <- suppressMessages(read_delim(path.2, delim = "\t"))
  sample.colnames <- c("_CHROM","POS","ID","REF","ALT","DP","AF","meth","unmeth")

  if(!all.equal(colnames(sample.1), sample.colnames))
    stop("Error: Incorrect colnames in sample.1")

  if(!all.equal(colnames(sample.2), sample.colnames))
    stop("Error: Incorrect colnames in sample.2")

  sample.types <- c("character","integer","character","character","character","integer","double","integer","integer")
  names(sample.types) <-sample.colnames

  if(!all.equal(sapply(sample.1, typeof), sample.types))
    stop("Error: Incorrect datatypes in sample.1")

  if(!all.equal(sapply(sample.2, typeof), sample.types))
    stop("Error: Incorrect datatypes in sample.2")

  if (!all(sample.1$POS > 0) & all(sample.1$`_CHROM` %in% paste0('chr', c(1:22, "X", "Y", "M"))) &
      all(sample.1$meth + sample.1$unmeth == sample.1$DP) )
    stop("Error: Incorrect data in sample.1")

  if (!all(sample.2$POS > 0) & all(sample.2$`_CHROM` %in% paste0('chr', c(1:22, "X", "Y", "M"))) &
      all(sample.2$meth + sample.2$unmeth == sample.2$DP) )
    stop("Error: Incorrect data in sample.2")

    sample.1 %>% inner_join(sample.2, by = c("_CHROM"="_CHROM", "POS"="POS")) %>%
    dplyr::select(-c(ID.x, REF.x, ALT.x,ID.y, REF.y, ALT.y, AF.x, AF.y ))  %>%
    gather(key, value, -`_CHROM`, -POS) %>%
    extract(key, c('type', "prob"), "(.*)\\.(.)") %>%
    spread(type, value) %>%
    mutate(meth.rate = meth/DP) %>%
    rename(no = DP,  chr = `_CHROM`, poz = POS) -> data
  rm(sample.1, sample.2, sample.types, sample.colnames)
  data
  }


