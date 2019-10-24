#' @title Evaluate Trait Predictability via the HAT Method
#' @description The HAT method is a fast algorithm for the ordinary cross validation. It is highly recommended for large dataset (Xu et al. 2017).
#' @param fix a design matrix of the fixed effects. If not passed, a vector of ones is added for the intercept.
#' @param y a vector of the phenotypic values.
#' @param kk a list of one or multiple kinship matrices.
#' @param nfold the number of folds, default is 5. For the HAT Method, nfold can be set as the sample size (leave-one-out CV) to avoid
#' variation caused by random partitioning of the samples, but it is not recommended for \code{\link{cv}}.
#' @param seed the random number, default is 123.
#'
#' @return  Trait predictability
#' @export
#'
#' @examples
#'
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
#' ##for the additive model
#' # predictability <- cv_fast(y = hybrid_phe[,3], kk = list(ka))
#'
#' ##for the additive-dominance model
#' # predictability <- cv_fast(y = hybrid_phe[,3], kk = list(ka,kd))
#'
#' @references
#' Xu, S. (2017) Predicted residual error sum of squares of mixed models: an application for genomic prediction. G3 (Bethesda) 7, 895-909.
cv_fast <- function(fix = NULL, y, kk, nfold = 5, seed = 123) {
  n <- length(y)
  y <- as.matrix(y)
  if (is.null(fix)) {
    fix <- matrix(1, n, 1)
  }
  ### Negative nloglikelihood Function
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
  if (parm$convergence != 0)
    stop("optim() failed to converge")
  v_phi <- 0
  for (p in 1:g) {
    v_phi <- v_phi + kk[[p]] * parm$par[p]
  }
  v_sigma <- diag(n) * parm$par[g + 1]
  v <- v_phi + v_sigma
  v_i <- solve(v)
  beta <- solve(t(fix) %*% v_i %*% fix, t(fix) %*% v_i %*% y)

  H <- v_phi %*% v_i
  r <- y - fix %*% beta
  SSP <- var(r) * (n - 1)
  rhat <- H %*% r
  ehat <- r - rhat
  if (nfold == n) {
    foldid <- c(1:n)
  } else {
    set.seed(seed)
    foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
  }
  PRESS <- 0
  for (i in 1:nfold) {
    indx <- which(foldid == i)
    nk <- length(indx)
    Hkk <- H[indx, indx]
    e <- solve(diag(nk) - Hkk) %*% ehat[indx]
    PRESS <- PRESS + sum(e^2)
  }
  R2 <- 1 - PRESS/SSP
  return(R2)
}
