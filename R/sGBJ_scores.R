#' Compute the pvalues and GBJ value associated with a pathway and survival
#'
#' @param surv a surv object of size n
#' @param counts_pathway a data frame of the counts for the particular pathway of interest of size nxp
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
#'  counts_pathway <- data.frame(P1 = rnorm(n = n),
#'                               P2 = rnorm(n = n))
#'
#'  sGBJ::sGBJ_scores(surv,counts_pathway, nperm = 2)
#'
#'  # with covariates
#'
#'  covariates <- data.frame(age = runif(n = n, 60, 90))
#'
#'  sGBJ_scores(surv,counts_pathway, nperm = 2, covariates = covariates)
sGBJ_scores=function(surv,counts_pathway,covariates=NULL,nperm=300){

  # computation of the score vector
  lsScores <- .survival_scores(counts_pathway = counts_pathway,
                               covariates = covariates,
                               surv = surv)

  # computation of the epsilon matrix by permutation
  epsilon <- .epsilon_matrix(Z = lsScores$Z,
                             nperm = nperm,
                             surv = surv,
                             counts_pathway = lsScores$updatedCount_pathway,
                             covariates = covariates,
                             datas = lsScores$datas)

  scores_GBJ=list(test_stats=lsScores$Z, cor_mat=epsilon)

  return(scores_GBJ)
}
