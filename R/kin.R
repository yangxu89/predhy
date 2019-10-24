#' @title Calculate Kinship Matrix
#' @description Calculate the additive and dominance kinship matrix (Xu et al. 2014).
#' @param gen a matrix for genotypes, coded as 1, 0, -1 for AA, Aa, aa. Each row represents an individual and each column represents a marker.
#'
#' @return a kinship matrix
#' @export
#'
#' @examples
#' ## random population with 100 lines and 1000 markers
#' gen <- matrix(rep(0,100*1000),100,1000)
#' gen <- apply(gen,2,function(x){x <- sample(c(-1,0,1), 100, replace = TRUE)})
#'
#' ## generate 100 × 100 kinship matrix
#' k <- kin(gen)
#' @references
#' Xu, S., Zhu, D. and Zhang, Q. (2014) Predicting hybrid performance in rice using genomic best linear unbiased prediction. Proc. Natl. Acad. Sci. USA 111, 12456-12461.
#'
kin <- function(gen) {
  gen <- as.matrix(gen)
  m <- ncol(gen)
  n <- nrow(gen)
  kk <- matrix(0, n, n)
  for (k in 1:m) {
    z <- gen[, k, drop = F]
    kk <- kk + z %*% t(z)
  }
  kk <- kk/m
  return(kk)
}

