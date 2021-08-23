#' .survival_scores
#'
#' @description Compute survival score.
#'
#' @param surv a surv object of size n
#' @param counts_pathway a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)
#'
#' @return A list of length 2 with the updated counts_pathway (same as counts_pathway but removing pathways for which survival model failed to converge) and the Z matrix.
.survival_scores <- function(counts_pathway, covariates = NULL, surv){
  remove_Z=NULL
  Z=numeric(ncol(counts_pathway))

  if(is.null(covariates)){
    datas=counts_pathway
    loopLength <- ncol(datas)
  } else {
    datas=cbind(covariates,counts_pathway)
    if (is.null(dim(covariates))){
      size_covariates=1
    }else{
      size_covariates=ncol(covariates)
    }
    loopLength <- ncol(datas)-size_covariates
  }

  for (i in 1:loopLength){

    if(is.null(covariates)){
      model=try(survival::coxph(surv~datas[,i]))
    } else {
      model=try(survival::coxph(surv~datas[,i+size_covariates]+datas[,2:size_covariates]))
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
    updatedCount_pathway = counts_pathway[,-remove_Z]
  } else {
    updatedCount_pathway = counts_pathway
  }

  return(list(Z = Z,
              updatedCount_pathway = updatedCount_pathway))
}
