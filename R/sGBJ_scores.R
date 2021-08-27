#' Compute the pvalues and GBJ value associated with a pathway and survival
#'
#' @param surv a surv object of size n
#' @param factor_matrix a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)
#' @param nperm number of permutations to perform to estimate the matrix epsilon (default=300)
#'
#' @importFrom stats cor
#' @export
#'
#' @return The GBJ value and it's pvalue associated
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

