#' .survival_scores
#'
#' @description Compute survival score.
#'
#' @param surv a surv object of size n
#' @param factor_matrix a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)
#'
#' @importFrom stats as.formula
#'
#' @return A list of length 3 with the updated factor_matrix (same as factor_matrix but removing columns for which survival model failed to converge), the Z matrix and the datas used to fit survival model.
.survival_scores <- function(factor_matrix, covariates = NULL, surv){
  remove_Z=NULL
  Z=numeric(ncol(factor_matrix))

  if(is.null(covariates)){
    datas=factor_matrix
  } else {
    datas=cbind(covariates,factor_matrix)
    if (is.null(dim(covariates))){
      size_covariates=1
    }else{
      size_covariates=ncol(covariates)
    }
  }

  for (i in 1:ncol(factor_matrix)){

    if(is.null(covariates)){
      model=try(survival::coxph(surv~datas[,i]))
    } else {
      vecCovariates <- "datas[,1]"
      if (size_covariates>1){
        for (j in 2:size_covariates){
        vecCovariates=paste(vecCovariates,"+datas[,",j,"]")
      }}
   #   vecCovariates <- colnames(datas)[1:size_covariates]
      vecPathway <- paste("datas[,",i+size_covariates,"]",sep="")
      formX <- paste(c(vecPathway, vecCovariates), collapse = " + ")
      form <- as.formula(paste0("surv ~ ", formX))

      model=try(survival::coxph(form, data = datas))
    }

    boolLengthModel <- length(model)>10
    if(boolLengthModel){
      if(is.null(covariates)){
        Z[i]=model$coefficients/summary(model)$coefficients[3]
      } else {
        Z[i]=model$coefficients[1]/summary(model)$coefficients[1,3]
      }
    }

    boolZna <- is.na(Z[i])
    boolZinf <- abs(Z[i]) == Inf
    boolRemove <- (!boolLengthModel) | boolZna | boolZinf

    if(boolRemove){
      remove_Z=c(remove_Z,i)
    }
  }

  if (length(remove_Z)!=0){
    Z=Z[-remove_Z]
    updatedFactor_matrix = factor_matrix[,-remove_Z]
  } else {
    updatedFactor_matrix = factor_matrix
  }

  return(list(Z = Z,
              updatedFactor_matrix = updatedFactor_matrix,
              datas = datas))
}
