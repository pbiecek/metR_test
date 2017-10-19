create.tiles.chr <- function(data.chr, gaps.length){

  require(data.table)
  require(dplyr)

  data.chr %>% distinct(poz) %>% arrange(poz) %>%
    mutate(diff = poz - lag(poz), id = seq_along(poz))  -> poz

  s1 <- c(1, which(poz$diff >= gaps.length))
  s2 <- c(s1[-1] - 1, nrow(poz))
  tiles <- 1:length(s1)
  to_bind <- data.table(s1,s2,tiles)
  setDT(poz)
  poz <- poz[to_bind, on=.(id >= s1, id <= s2), nomatch = 0, .(poz, tiles)]

  data.chr %>% left_join(poz, by = c("poz"= "poz"))
}
