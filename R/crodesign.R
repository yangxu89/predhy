#' @title Generate Mating Design
#' @description Generate a mating design for a subset of crosses based on a balanced random partial rectangle cross-design (BRPRCD) (Xu et al. 2016).
#' @param d an integer denoting 1/d percentage of crosses to be evaluated in the field.
#' @param male_name a character string for the names of male parents.
#' @param female_name a character string for the names of male parents.
#' @param seed the random number, default is 123.
#'
#' @return  A data frame of mating design result with three columns. The first column is "crossID", the second column is the "male_Name" and the third column is the "female_Name".
#' @export
#'
#' @examples
#'
#' ## generate a mating design with 100 male parents and 150 female parents
#' ## for 1/d = 1/50 percentage of crosses to be evaluated in the field.
#' ## the total number of potential crosses is 100 × 150 = 15000.
#' ## The number of crosses to be field evaluated is 15000 × (1/50) = 300.
#'
#' male_name <- paste("m", 1:100, sep = "")
#' female_name <- paste("f", 1:150, sep = "")
#' design <- crodesign(d = 50, male_name, female_name)
#'
#' @references
#' Xu, S., Xu, Y., Gong, L. and Zhang, Q. (2016) Metabolomic prediction of yield in hybrid rice. Plant J. 88, 219-227.
#'
crodesign <- function(d, male_name, female_name, seed = 123) {
  m <- length(male_name)
  male_id <- c(1:m)
  n <- length(female_name)
  n1 <- ifelse(m < n, n, m)
  temp <- ceiling(n1/d)
  n1 <- temp * d + 1
  a <- matrix(data = 1:(n1 * n1) * 0, nrow = n1)
  for (i in 1:n1) {
    j <- i
    while (j <= n1) {
      a[i, j] <- 1
      a[j, i] <- 1
      j <- j + d
    }
  }
  b <- a[1:n, 1:m]

  crossID <- NULL
  parent1 <- NULL
  parent2 <- NULL
  t <- 0
  for (p in 1:n) {
    for (q in 1:m) {
      if (b[p, q] == 1) {
        t <- t + 1
        crossID <- c(crossID, t)
        parent1 <- c(parent1, p)
        parent2 <- c(parent2, q)
      }
    }
  }
  set.seed(seed)
  nfemale <- sample(1:n)
  nmale <- sample(1:m)
  female_Name <- female_name[nfemale][parent1]
  male_Name <- male_name[nmale][parent2]
  crossdesign <- data.frame(crossID, male_Name, female_Name)
  return(crossdesign)
}

