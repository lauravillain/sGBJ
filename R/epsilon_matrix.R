epsilon_matrix <- function(Z, perm, surv, counts_pathway, covariates = NULL){

  Z_matrix<- matrix(nrow=(length(Z)), ncol=nperm)
  Z_matrix[,1]=Z
  i=2
  perm=sample(length(surv))
  surv_perm=surv[perm]

  if(is.null(covariates)){
    datas_perm=counts_pathway
  } else {
    covariates_perm=covariates[perm,]
    datas_perm=cbind(covariates_perm,counts_pathway)
  }

  while(i<=nperm){
    perm_OK=TRUE
    for (j in 1:(length(Z))){
      if(is.null(covariates)){
        model=try(survival::coxph(surv_perm~datas_perm[,j]))
      } else {
        model=try(survival::coxph(surv_perm~datas[,j+(size_covariates)]+datas_perm[,2:(size_covariates)], data = datas_perm))
      }
      if (length(model)>10){
        if(is.null(covariates)){
          Z_matrix[j,i]=model$coefficients/summary(model)$coefficients[3]
        } else {
          Z_matrix[j,i]=model$coefficients[1]/summary(model)$coefficients[1,3]
        }
        if(is.na(Z_matrix[j,i])){
          perm_OK=FALSE
        } else{
          if (abs(Z_matrix[j,i])==Inf){
            perm_OK=FALSE
          }
        }
      } else {
        perm_OK=FALSE
      }
    }
    if (perm_OK==TRUE){
      i=i+1
    }
  }
}
