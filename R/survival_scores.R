#' .survival_scores
#'
#' @description Compute survival score.
#'
#' @param surv a surv object of size n
#' @param counts_pathway a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)
#'
#' @importFrom stats as.formula
#'
#' @return A list of length 3 with the updated counts_pathway (same as counts_pathway but removing pathways for which survival model failed to converge), the Z matrix and the datas used to fit survival model.
.survival_scores <- function(counts_pathway, covariates = NULL, surv){
  remove_Z=NULL
  Z=numeric(ncol(counts_pathway))
  
  if(is.null(covariates)){
    datas=counts_pathway
  } else {
    datas=cbind(covariates,counts_pathway)
    if (is.null(dim(covariates))){
      size_covariates=1
    }else{
      size_covariates=ncol(covariates)
    }
  }
  
  for (i in 1:ncol(counts_pathway)){
    
    if(is.null(covariates)){
      model=try(survival::coxph(surv~datas[,i]))
    } else {
      vecCovariates <- colnames(datas)[1:size_covariates]
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
    updatedCount_pathway = counts_pathway[,-remove_Z]
  } else {
    updatedCount_pathway = counts_pathway
  }
  
  return(list(Z = Z,
              updatedCount_pathway = updatedCount_pathway,
              datas = datas))
}