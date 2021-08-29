#' Compute the sGBJ statistic along with its p-value quantifying the association
#' between a gene set and survival outcome
#'
#' @param surv a \code{\link[survival]{Surv}} object of length \code{n}
#' @param factor_matrix a \code{n x p} \code{data.frame} of the expression for the
#' particular gene set of interest being tested
#' @param covariates a \code{n x l} matrix of the covariates to adjust upon. Default is \code{NULL}
#' @param nperm number of permutations performed to estimate the \code{epsilon} matrix.
#' Default is \code{300}.
#'
#' @importFrom stats cor
#' @export
#'
#' @return a list containing the sGBJ statistic estimation and its associated p-value
#' @examples
#'  n <- 100
#'  surv_data <- data.frame(Time = runif(n = n, min = 0, max = 100),
#'                          event = rbinom(n = n, size = 1, prob = 0.5))
#'  surv <- survival::Surv(time = surv_data$Time, event = surv_data$event)
#'
#'  factor_matrix <- data.frame(P1 = rnorm(n = n),
#'                               P2 = rnorm(n = n))
#'
#'  sGBJ::sGBJ_scores(surv,factor_matrix, nperm = 2)
#'
#'  # with covariates
#'
#'  covariates <- data.frame(age = runif(n = n, 60, 90))
#'
#'  sGBJ_scores(surv,factor_matrix, nperm = 2, covariates = covariates)
sGBJ_scores=function(surv,factor_matrix,covariates=NULL,nperm=300){

  # computation of the score vector
  lsScores <- .survival_scores(factor_matrix = factor_matrix,
                               covariates = covariates,
                               surv = surv)

  # computation of the epsilon matrix by permutation
  epsilon <- .epsilon_matrix(Z = lsScores$Z,
                             nperm = nperm,
                             surv = surv,
                             factor_matrix = lsScores$updatedFactor_matrix,
                             covariates = covariates,
                             datas = lsScores$datas)

  scores_GBJ=list(test_stats=lsScores$Z, cor_mat=epsilon)

  return(scores_GBJ)
}

