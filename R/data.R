#' Genetics datasets
#'
#' Eight different genetics datasets \code{b1, b2, b3, g1, g2, w1, w2, w3} are used in this paper
#'
#' @docType data
#' @keywords datasets
#' @usage data(b1)
#' @name genes
#' @format Each of the eight datasets contains 40 columns. Columns 2 through 40 contain the 20 nucleotides to the left and 19 nucleotides to the right of a particular nucleotide. The first column, \code{count}, contains the number of times those sequences were found in different RNA sequences.
#'
#' The datasets beginning with \code{b} are human and originated with the Burge group \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2593745/}.
#'
#' The datasets beginning with \code{g} are from mice and due to the Grimmond group \url{http://www.nature.com/nmeth/journal/v5/n7/full/nmeth.1223.html}.
#'
#' The datasets beginning with \code{w} are also from mice and due to the Wold group \url{http://www.nature.com/nmeth/journal/v5/n7/full/nmeth.1226.html}.
#'
#' @source These data were downloaded from Prof. Jun Li and then processed using his package \code{\link[=https://cran.r-project.org/web/packages/mseq/index.html]{mseq}} \url{http://www3.nd.edu/~jli9/mseq/data_top100.zip}.
NULL

#' @rdname genes
'b1'

#' @rdname genes
'b2'

#' @rdname genes
'b3'

#' @rdname genes
'g1'

#' @rdname genes
'g2'

#' @rdname genes
'w1'

#' @rdname genes
'w2'

#' @rdname genes
'w3'
