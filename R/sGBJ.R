#' Compute the pvalues and GBJ value associated with a pathway and survival
#'
#' @param surv a surv object of size n
#' @param counts_pathway a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)
#' @param nperm number of permutations to perform to estimate the matrix epsilon (default=300)
#'
#'
#' @export
#' @return The GBJ value and it's pvalue associated
#' @examples
#' n <- 100
#' surv_data <- data.frame(Time = runif(n = n, min = 0, max = 100),
#'                         event = rbinom(n = n, size = 1, prob = 0.5))
#' surv <- survival::Surv(time = surv_data$Time, event = surv_data$event)
#'
#' counts_pathway <- data.frame(P1 = rnorm(n = n),
#'                              P2 = rnorm(n = n))
#'
#' sGBJ::sGBJ(surv,counts_pathway, nperm = 2)
sGBJ=function(surv,counts_pathway,covariates=NULL,nperm=300){

  scores_GBJ=sGBJ_scores(surv,counts_pathway,covariates,nperm)
  GBJOut <- GBJ::GBJ(test_stats=scores_GBJ$test_stats, cor_mat=scores_GBJ$cor_mat)
  return(GBJOut)
}
