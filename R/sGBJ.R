#' Compute the pvalues and GBJ value associated with a pathway and survival
#'
#' @param surv a surv object of size n
#' @param counts_pathway a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)

#' @return The GBJ value and it's pvalue associated
#' @examples
#' sGBJ(surv,counts_pathway)

sGBJ=function(surv,counts_pathway,covariates=NULL){


  remove_Z=NULL
  Z=numeric(dim(counts_pathway)[2])

  if (is.null(covariates)){
    datas=cbind(surv,counts_pathway)
    for (i in 1:(dim(datas)[2]-1)){
      model=try(coxph(surv~datas[,i+1]))
      if (length(model)>10){
        Z[i]=model$coefficients/summary(model)$coefficients[3]
        if (is.na(Z[i])){
          remove_Z=c(remove_Z,i)
        }else{
          if (abs(Z[i])==Inf){
            remove_Z=c(remove_Z,i)
          }}
      }else{
        remove_Z=c(remove_Z,i)
      }
    }
    if (length(remove_Z)!=0){
      Z=Z[-remove_Z]
      counts_pathway=counts_pathway[,-remove_Z]
      datas=cbind(surv,covariates,counts_pathway)
      datas=as.data.frame(datas)
    }

    Z_matrix<- matrix(nrow=(length(Z)), ncol=300)
    Z_matrix[,1]=Z
    i=2
    while(i<=300){
      perm=sample(length(surv))
      surv_perm=surv[perm]
      datas_perm=cbind(surv_perm,counts_pathway)
      perm_OK=TRUE
      for (j in 1:(length(Z))){
        model=try(coxph(surv_perm~datas_perm[,j+1]))
        if (length(model)>10){
          Z_matrix[j,i]=model$coefficients/summary(model)$coefficients[3]
          if( is.na(Z_matrix[j,i])){
            perm_OK=FALSE
          }else{
            if (abs(Z_matrix[j,i])==Inf){
              perm_OK=FALSE
            }}


        }else{perm_OK=FALSE}

      }
      if (perm_OK==TRUE){
        i=i+1
      }
    }

  }else{
    datas=cbind(surv,covariates,counts_pathway)
    size_covariates=dim(covariates)[2]
    for (i in 1:(dim(datas)[2]-(size_covariates+1))){
      model=try(coxph(surv~datas[,i+(size_covariates+1)]+datas[,2:(1+size_covariates)]))
      if (length(model)>10){
        Z[i]=model$coefficients[1]/summary(model)$coefficients[1,3]
        if (is.na(Z[i])){
          remove_Z=c(remove_Z,i)
        }else{
          if (abs(Z[i])==Inf){
            remove_Z=c(remove_Z,i)
          }}
      }else{
        remove_Z=c(remove_Z,i)
      }
    }
    if (length(remove_Z)!=0){
      Z=Z[-remove_Z]
      counts_pathway=counts_pathway[,-remove_Z]
      datas=cbind(surv,covariates,counts_pathway)
      datas=as.data.frame(datas)
    }

    Z_matrix<- matrix(nrow=(length(Z)), ncol=300)
    Z_matrix[,1]=Z
    i=2
    while(i<=300){
      perm=sample(length(surv))
      surv_perm=surv[perm]
      covariates_perm=covariates[perm,]
      datas_perm=cbind(surv_perm,covariates_perm,counts_pathway)
      perm_OK=TRUE
      for (j in 1:(length(Z))){
        model=try(coxph(surv_perm~datas[,j+(size_covariates+1)]+datas_perm[,2:(1+size_covariates)], data = datas_perm))
        if (length(model)>10){
          Z_matrix[j,i]=model$coefficients[1]/summary(model)$coefficients[1,3]
          if( is.na(Z_matrix[j,i])){
            print(c("Z matrix na",i,j))
            print(c("coeffs",model$coefficients,sqrt(model$var)))
            perm_OK=FALSE
          }else{
            if (abs(Z_matrix[j,i])==Inf){
              perm_OK=FALSE
            }}


        }else{perm_OK=FALSE}

      }
      if (perm_OK==TRUE){
        i=i+1
      }
    }

  }

  epsilon=cor(t(Z_matrix))

  scores_GBJ=list(test_stats=Z,cor_mat=epsilon)
  GBJOut <- GBJ(test_stats=Z, cor_mat=epsilon)
  return(GBJOut)
}
