#' @title Infer Genotype of Hybrids
#' @description Infer additive and dominance genotypes of hybrids based on their parental genotypes.
#' @param inbred_gen a matrix for genotypes of parental lines in numeric format, coded as 1, 0 and -1. The row.names of inbred_gen must be provied. It can be obtained from the original genotype using \code{\link{convertgen}} function.
#' @param hybrid_phe a data frame with three columns. The first column and the second column are the names of male and female parents of the corresponding hybrids, respectively; the third column is the phenotypic values of hybrids.
#'     The names of male and female parents must match the rownames of inbred_gen. Missing (NA) values are not allowed.
#'
#' @return A list with following information is returned:
#'
#'     $add  additive genotypes of hybrids
#'
#'     $dom  dominance genotypes of hybrids
#'
#' @export
#'
#' @examples
#'
#' ## load example data from hypred package
#' # data(hybrid_phe)  ## remove one '#' mark to run the examples
#' # head(hybrid_phe)
#' # data(input_geno)
#'
#' ## convert original genotype, and the type is selected based on the format of genotype
#' # inbred_gen <- convertgen(input_geno, type = "hmp2")
#'
#' # gena <- infergen(inbred_gen, hybrid_phe)$add
#' # gend <- infergen(inbred_gen, hybrid_phe)$dom
#'
infergen <- function(inbred_gen, hybrid_phe) {
  p_name <- row.names(inbred_gen)
  p1 <- as.character(hybrid_phe[, 1])
  p2 <- as.character(hybrid_phe[, 2])
  inbred_gen <- as.matrix(inbred_gen)
  p1_gen <- inbred_gen[match(p1, p_name), ]
  p2_gen <- as.matrix(inbred_gen[match(p2, p_name), ])
  f1_gena <- (p1_gen + p2_gen)/2
  f1_gend <- abs(p1_gen - p2_gen)/2
  return(list(add = f1_gena, dom = f1_gend))
}
