#' @title Predict the Performance of Hybrids
#' @description Predict all potential crosses of a given set of parents using a subset of crosses as the training sample.
#' @param inbred_gen a matrix for genotypes of parental lines in numeric format, coded as 1, 0 and -1. The row.names of inbred_gen must be provied. It can be obtained from the original genotype using  \code{\link{convertgen}} function.
#' @param hybrid_phe a data frame with three columns. The first column and the second column are the names of male and female parents of the corresponding hybrids, respectively; the third column is the phenotypic values of hybrids.
#' The names of male and female parents must match the rownames of inbred_gen. Missing (NA) values are not allowed.
#' @param predparent_gen a matrix for genotypes of a given parental lines. All potential crosses derived from these parental lines will be predicted. Default is the same as inbred_gen.
#' @param fix a design matrix of the fixed effects for the parental lines. If not passed, a vector of ones is added for the intercept.
#' @param fixnew a design matrix of the fixed effects for all potential crosses. If not passed, a vector of ones is added for the intercept.
#' @param model the model of prediction. There are two options: model = "AD" for the additive-dominance model, model = "A" for the additive model. Default is model = "AD".
#' @param select the selection of hybrids based on the prediction results. There are three options: select = "all", which selects all potential crosses.
#' select = "top", which selects the top n crosses. select = "top", which selects the bottom n crosses. The n is determined by the param number.
#' @param number the number of selected top or bottom hybrids, only when select = "top" or select = "bottom".
#'
#' @return A data frame of prediction result with two columns. The first column denotes the names of male and female parents of the predicted hybrids,
#' and the second column denotes the phenotypic values of the predicted hybrids.
#' @export
#'
#' @examples
#' ## load example data from hypred package
#' # data(hybrid_phe) ## remove one '#' mark to run the examples
#' # data(input_geno)
#' # inbred_gen <- convertgen(input_geno, type = "hmp2")
#'
#' ## to save time, only predict 45 crosses derived from the first 10 lines of parental lines
#' ## select all hybrids with additive-dominance model
#' # pred1 <- predhybrid(inbred_gen, hybrid_phe, predparent_gen = inbred_gen[c(1:10),], model = "AD")
#'
#' ## select top 20 hybrids with additive model
#' # pred2 <- predhybrid(inbred_gen, hybrid_phe, predparent_gen = inbred_gen[c(1:10),],
#' # select = "top", number = 10, model = "A")
#'
predhybrid <- function(inbred_gen, hybrid_phe, predparent_gen = inbred_gen, fix = NULL,
                    fixnew = NULL, model = "AD", select = "all", number = NULL) {
  n1 <- nrow(hybrid_phe)
  y <- hybrid_phe[, 3]

  if (model == "AD") {
    gena <- infergen(inbred_gen,hybrid_phe)[[1]]
    ka <- kin(gena)
    gend <- infergen(inbred_gen,hybrid_phe)[[2]]
    kd <- kin(gend)
    parm <- mixed(fix = fix, y = y, kk = list(ka, kd))
    v_i <- parm$v_i
    beta <- parm$beta
    va <- parm$var[1]
    vd <- parm$var[2]
    ve <- parm$ve

    ka21 <- NULL
    kd21 <- NULL
    phe_name <- NULL
    pred_name <- row.names(predparent_gen)
    for (i in 1:(nrow(predparent_gen) - 1)) {
      ha1 <- (predparent_gen[i, ] + predparent_gen[-(1:i), ])/2
      ka2 <- tcrossprod(ha1, gena)/ncol(gena)
      ka21 <- rbind(ka21, ka2)
      hd1 <- abs(predparent_gen[i, ] - predparent_gen[-(1:i), ])/2
      kd2 <- tcrossprod(hd1, gend)/ncol(gend)
      kd21 <- rbind(kd21, kd2)
      test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
      phe_name <- c(phe_name, test_name)
    }
    n2 <- nrow(kd21)
    if (is.null(fixnew)) {
      fixnew <- matrix(1, n2, 1)
    } else {
      fixnew <- as.matrix(fixnew)
    }
    if (is.null(fix)) {
      fix <- matrix(1, n1, 1)
    } else {
      fix <- as.matrix(fix)
    }
    G21 <- ka21 * va + kd21 * vd
    pred_phe <- fixnew %*% beta + G21 %*% v_i %*% (y - fix %*% beta)
    row.names(pred_phe) <- phe_name
  }

  if (model == "A") {
    gena <- infergen(inbred_gen,hybrid_phe)[[1]]
    ka <- kin(gena)
    parm <- mixed(fix = fix, y = y, kk = list(ka))
    v_i <- parm$v_i
    beta <- parm$beta
    va <- parm$var[1]
    ve <- parm$ve
    ka21 <- NULL
    phe_name <- NULL
    pred_name <- row.names(predparent_gen)
    for (i in 1:(nrow(predparent_gen) - 1)) {
      ha1 <- (predparent_gen[i, ] + predparent_gen[-(1:i), ])/2
      ka2 <- tcrossprod(ha1, gena)/ncol(gena)
      ka21 <- rbind(ka21, ka2)
      test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
      phe_name <- c(phe_name, test_name)
    }
    n2 <- nrow(ka21)
    if (is.null(fixnew)) {
      fixnew <- matrix(1, n2, 1)
    } else {
      fixnew <- as.matrix(fixnew)
    }
    if (is.null(fix)) {
      fix <- matrix(1, n1, 1)
    } else {
      fix <- as.matrix(fix)
    }
    G21 <- ka21 * va
    pred_phe <- fixnew %*% beta + G21 %*% v_i %*% (y - fix %*% beta)
    row.names(pred_phe) <- phe_name
  }

  if (select == "all") {
    pred_phe_select <- as.data.frame(pred_phe)
    colnames(pred_phe_select) <- paste("all_", nrow(pred_phe), sep = "")
  } else if (select == "top") {
    pred_phe_select <- as.data.frame(sort(pred_phe[, 1], decreasing = T)[c(1:number)])
    names(pred_phe_select) <- paste("top_", number, sep = "")
  } else if (select == "bottom") {
    pred_phe_select <- as.data.frame(sort(pred_phe[, 1], decreasing = F)[c(1:number)])
    colnames(pred_phe_select) <- paste("bottom_", number, sep = "")
  }
  return(pred_phe_select)
}

