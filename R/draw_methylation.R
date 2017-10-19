#' Plot DMR region
#'
#' Visualizing DMR region by plotting methylation rate within two probes
#' @param data There are two options:  1. dataframe with specyfic columns: chr, poz, prob, no, meth, unmeth, meth.rate.
#' This dataframe is result of function preprocessing.
#' 2. dataframe with specyfic columns: chr, poz, prob, no, meth, unmeth, meth.rate, tiles and possible tiles.common columns. This dataframe is result of function create.tiles.min.gap or
#' create.tiles.fixed.length.
#' @param chr chromosom name of region that are being plotted
#' @param start minimum position of region that are being plotted
#' @param end maximum position of region that are being plotted
#' @return ggplot object with visualization of regions
#' @export
draw.methylation <- function(data, chromosom, start, end, plot.title = 26, axis.title.x  = 23,
                             axis.title.y = 23, legend.position = 'right', axis.text.x = 20, axis.text.y = 20,
                             legend.text = 16, legend.title = 18){

  require(ggplot2)
  require(dplyr)
  require(tidyr)


  if(!zapsmall(start, 15) == round(start) | start <= 0)
    stop("Error: Incorrect start argument")

  if(!zapsmall(end, 15) == round(end) | end <= 0)
    stop("Error: Incorrect end argument")

  if(!chromosom %in% paste0('chr', c(1:22, "X", "Y", "M")))
    stop("Error: Incorrect chr argument")

  data.colnames <- c('chr', 'poz', 'prob', 'no', 'meth', 'unmeth', 'meth.rate')

  if ('tiles.common' %in% colnames(data)) data <- data[, -ncol(data)]
  if ('tiles' %in% colnames(data)) data <- data[, -ncol(data)]

  if(!all.equal(colnames(data), data.colnames))
    stop("Error: Incorrect colnames in data")

  data.types <- c("character","integer","character","integer","integer", "integer", 'double')
  names(data.types) <- data.colnames

  if(!all.equal(sapply(data, typeof), data.types))
    stop("Error: Incorrect datatypes in data")

  if (!all(data$poz > 0) | !all(data$chr %in% paste0('chr', c(1:22, "X", "Y", "M"))) |
      !all(data$meth + data$unmeth == data$no) | !all(round(data$meth /data$no,3) == round(data$meth.rate,3)))
    stop("Error: Incorrect values in data")


  DT <- data %>% filter(chr == chromosom, poz >=start, poz <= end)

  DT1 <- DT %>% mutate(cov = round(log(no),1)) %>%  dplyr::select(poz, chr, prob, cov) %>% spread(prob, cov) %>%
    rename(cov_I = x, cov_II = y)

  DT  %>% dplyr::select(poz, chr, meth.rate, prob) %>% spread(prob, meth.rate) %>%
    rename(met_I = x, met_II = y) %>% left_join(DT1, by = c("poz", "chr")) -> DT

  p <- ggplot(DT)  + xlim(start,end)  +
    geom_point(DT,map=aes(x=poz, y=met_I, colour=cov_I), size = 15, shape = 20) +
    scale_color_gradient(low="pink1", high="red", name = "log(covI)") +
    geom_point(DT,map=aes(x=poz, y=met_II, fill=cov_II), size = 9, shape = 21) +
    scale_fill_gradient(low="cadetblue1", high="royalblue4", name = "log(covII)")

  x = as.matrix(DT[,1])
  y1 = as.matrix(DT[,'met_I'])
  y2 = as.matrix(DT[,'met_II'])

  p <- p + geom_segment(aes(x = x, y = y1, xend = x, yend = y2) )  +
    ggtitle(paste0("Methylation rate within two probes in ", chromosom)) +
    xlab("DNA position") + ylab("Methylation rate") + theme_minimal() +
    theme(
      plot.title = element_text(size=plot.title),
      axis.title.x = element_text(size=axis.title.x),
      axis.title.y = element_text(size=axis.title.y),
      legend.position=legend.position,
      axis.text.x = element_text(size=axis.text.x),
      axis.text.y = element_text(size=axis.text.y),
      legend.text = element_text(size = legend.text),
      legend.title = element_text(size = legend.title),
      plot.margin = unit(c(1,1,1,1), "cm"))
  p

}

