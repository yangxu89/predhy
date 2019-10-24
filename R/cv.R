#' @title Evaluate Trait Predictability via Cross Validation
#' @description Evaluate trait predictability of the GBLUP method via k-fold cross validation. For k-fold cross validation, the sample is randomly divided into k equal
#' sized parts and each part is predicted once using parameters estimated based on the other k – 1 parts. The trait predictability is defined as the squared Pearson correlation coefficient between the observed and the predicted trait values.
#' @param fix a design matrix of the fixed effects. If not passed, a vector of ones is added for the intercept.
#' @param y a vector of the phenotypic values.
#' @param kk a list of one or multiple kinship matrices.
#' @param nfold the number of folds. Default is 5.
#' @param seed the random number. Default is 123.
#'
#' @return  Trait predictability
#' @export
#'
#' @examples
#'
#' ## load example data from hypred package
#' # data(hybrid_phe) ## remove one '#' mark to run the examples
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
#' # predictability <- cv(y = hybrid_phe[,3], kk = list(ka))
#'
#' ## for the additive-dominance model
#' # predictability <- cv(y = hybrid_phe[,3], kk = list(ka,kd))
#'
cv <- function(fix = NULL, y, kk, nfold = 5, seed = 123) {
  n <- length(y)
  y <- as.matrix(y)
  if (is.null(fix)) {
    fix <- matrix(1, n, 1)
  } else {
    fix <- as.matrix()
  }
  set.seed(seed)
  foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
  yobs <- NULL
  yhat <- NULL
  fold <- NULL
  id <- NULL
  k11 <- list()
  k21 <- list()
  g <- length(kk)
  for (k in 1:nfold) {
    i1 <- which(foldid != k)
    i2 <- which(foldid == k)
    x1 <- fix[i1, , drop = F]
    y1 <- y[i1, , drop = F]

    x2 <- fix[i2, , drop = F]
    y2 <- y[i2, , drop = F]
    # n1 <- length(y1)
    for (i in 1:g) {
      k11[[i]] <- kk[[i]][i1, i1]
      k21[[i]] <- kk[[i]][i2, i1]
    }
    parm <- mixed(fix = x1, y = y1, kk = k11)
    G21 <- 0
    for (i in 1:g) {
      G21 <- G21 + k21[[i]] * parm$var[i]
    }
    v_i <- parm$v_i
    beta <- parm$beta
    y3 <- x2 %*% beta + G21 %*% v_i %*% (y1 - x1 %*% beta)
    fold <- c(fold, rep(k, length(y2)))
    yobs <- c(yobs, y2)
    yhat <- c(yhat, y3)
    id <- c(id, i2)
  }
  pred <- data.frame(fold, id, yobs, yhat)
  R2 <- cor(pred$yobs, pred$yhat)^2
  return(R2)
}

