% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tiles_max_gap.R
\name{create.tiles.max.gap}
\alias{create.tiles.max.gap}
\title{Create regions for testing methylation difference}
\usage{
create.tiles.max.gap(data, gaps.length)
}
\arguments{
\item{data}{dataframe with specyfic columns: chr, poz, prob, no, meth, unmeth, meth.rate. This dataframe is result of function preprocessing.}

\item{gaps.length}{integer number that specifes maximum difference position between the same methylation regions}
}
\value{
data.frame from parameter data with extra column tiles that is region id number within chromosomes
}
\description{
The same regions are observation that maximum difference position is gaps.length argument
}
