#' .epsilon_matrix
#'
#' @description Compute the epsilon matrix by permutation for the sGBJ_scores() function.
#'
#' @param Z The score vector returned by .survival_scores() function.
#' @param surv a surv object of size n
#' @param counts_pathway a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)
#' @param nperm number of permutations to perform to estimate the matrix epsilon (default=300)
#'
#' @return The epsilon matrix.
.epsilon_matrix <- function(Z, nperm, surv, counts_pathway, covariates = NULL){

  # pre compute usefull parameter/data
  perm=sample(length(surv))
  surv_perm=surv[perm]

  if(is.null(covariates)){
    datas_perm=counts_pathway
  } else {
    covariates_perm=covariates[perm,]
    datas_perm=cbind(covariates_perm,counts_pathway)
    if (is.null(dim(covariates))){
      size_covariates=1
    }else{
      size_covariates=ncol(covariates)
    }
  }

  # build Z_matrix
  Z_matrix<- matrix(nrow=(length(Z)), ncol=nperm)
  Z_matrix[,1]=Z
  i=2

  # fill Z_matrix
  while(i<=nperm){
    perm_OK=TRUE
    for (j in 1:(length(Z))){
      if(is.null(covariates)){
        model=try(survival::coxph(surv_perm~datas_perm[,j]))
      } else {
        model=try(survival::coxph(surv_perm~datas[,j+size_covariates]+datas_perm[,2:size_covariates], data = datas_perm))
      }

      boolLengthModel <- length(model)>10
      if(boolLengthModel){
        if(is.null(covariates)){
          Z_matrix[j,i]=model$coefficients/summary(model)$coefficients[3]
        } else {
          Z_matrix[j,i]=model$coefficients[1]/summary(model)$coefficients[1,3]
        }
      }

      boolZna <- !is.na(Z_matrix[j,i])
      boolZinf <- abs(Z_matrix[j,i]) != Inf

      perm_OK <- boolLengthModel & boolZna & boolZinf
    }
    i <- i + perm_OK
  }

  epsilon=cor(t(Z_matrix))

  return(epsilon)
}
