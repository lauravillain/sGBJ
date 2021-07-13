#' Compute the pvalues and GBJ value associated with a pathway and survival
#'
#' @param surv a surv object of size n
#' @param counts_pathway a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)
#' @param nperm number of permutations to perform to estimate the matrix epsilon (default=300)


#' @return The GBJ value and it's pvalue associated
#' @examples
#' sGBJ(surv,counts_pathway)

sGBJ=function(surv,counts_pathway,covariates=NULL,nperm=300){

  scores_GBJ=sGBJ_scores(surv,counts_pathway,covariates=NULL,nperm=300)
  GBJOut <- GBJ(test_stats=Z, cor_mat=epsilon)
  return(GBJOut)
}
