#' @title Solve Mixed Model
#' @description Solve linear mixed model using restricted maximum likelihood (REML). Multiple variance components can be estimated.
#' @param fix a design matrix of the fixed effects. If not passed, a vector of ones is added for the intercept.
#' @param y a vector of the phenotypic values.
#' @param kk a list of one or multiple kinship matrices.
#'
#' @return A list with following information is returned:
#'     $v_i  the inverse of the phenotypic variance-covariance matrix.
#'     $var  estimated variance components of genetic effects
#'     $ve   estimated residual variance
#'     $beta estimated fixed effects
#' @export
#'
#' @examples
#' ## load example data from hypred package
#' # data(hybrid_phe)  ## remove one '#' mark to run the examples
#' # data(input_geno)
#'
#' ## convert original genotype
#' # inbred_gen <- convertgen(input_geno, type = "hmp2")
#'
#' ## infer the additive and dominance genotypes of hybrids
#' # gena <- infergen(inbred_gen, hybrid_phe)$add
#' # gend <- infergen(inbred_gen, hybrid_phe)$dom
#'
#' ## calculate the additive and dominance kinship matrix
#' # ka <- kin(gena)
#' # kd <- kin(gend)
#'
#' ## for the additive model
#' # parm <- mixed(y = hybrid_phe[,3], kk = list(ka))
#'
#' ## for the additive-dominance model
#' # parm <- mixed(y = hybrid_phe[,3], kk = list(ka, kd))
#'
#' @references
#' Xu, S., Zhu, D. and Zhang, Q. (2014) Predicting hybrid performance in rice using genomic best linear unbiased prediction. Proc. Natl. Acad. Sci. USA 111, 12456-12461.
#'
mixed <- function(fix = NULL, y, kk) {
  n <- length(y)
  y <- as.matrix(y)
  if (is.null(fix)) {
    fix <- matrix(1, n, 1)
  } else {
    fix <- as.matrix(fix)
  }
  ##Negative nloglikelihood Function
  g <- length(kk)
  nloglik_REML <- function(pm) {
    v_phi <- 0
    for (p in 1:g) {
      v_phi <- v_phi + kk[[p]] * pm[p]
    }
    v_sigma <- diag(n) * pm[g + 1]
    v <- v_phi + v_sigma + diag(1e-09, n)
    v_i <- solve(v, tol = -50)
    beta <- solve(t(fix) %*% v_i %*% fix, t(fix) %*% v_i %*% y)
    nloglik <- 0.5 * (unlist(determinant(v))[[1]] + unlist(determinant(t(fix) %*%
                                                                         v_i %*% fix))[[1]] + t(y - fix %*% beta) %*% v_i %*% (y - fix %*% beta))
    return(nloglik)
  }
  parm0 <- rep(1, g + 1)
  parm <- optim(par = parm0, fn = nloglik_REML, method = "L-BFGS-B", hessian = FALSE,
                lower = 0)
  v_phi <- 0
  for (p in 1:g) {
    v_phi <- v_phi + kk[[p]] * parm$par[p]
  }
  v_sigma <- diag(n) * parm$par[g + 1]
  v <- v_phi + v_sigma
  v_i <- solve(v, tol = -50)
  beta <- solve(t(fix) %*% v_i %*% fix, t(fix) %*% v_i %*% y)
  res <- list(v_i = v_i, var = parm$par[-(g + 1)], ve = parm$par[g + 1], beta = beta)
  return(res)
}
