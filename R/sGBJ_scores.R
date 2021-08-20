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
#' n <- 5
#' surv_data <- data.frame(Time = runif(n = n, min = 0, max = 100),
#'                         event = rbinom(n = n, size = 1, prob = 0.5))
#' surv <- survival::Surv(time = surv_data$Time, event = surv_data$event)
#'
#' counts_pathway <- data.frame(P1 = rnorm(n = n),
#'                              P2 = rnorm(n = n))
#'
#' sGBJ::sGBJ_scores(surv,counts_pathway, nperm = 2)
sGBJ_scores=function(surv,counts_pathway,covariates=NULL,nperm=300){


  remove_Z=NULL
  Z=numeric(dim(counts_pathway)[2])

  # Computation of scores if there is no covariates

  if (is.null(covariates)){
    datas=counts_pathway
    for (i in 1:(dim(datas)[2])){
      model=try(survival::coxph(surv~datas[,i]))
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
      datas=cbind(covariates,counts_pathway)
      datas=as.data.frame(datas)
    }

    #Computation of the epsilon matrix by permutation
    Z_matrix<- matrix(nrow=(length(Z)), ncol=nperm)
    Z_matrix[,1]=Z
    i=2
    while(i<=nperm){
      perm=sample(length(surv))
      surv_perm=surv[perm]
      datas_perm=counts_pathway
      perm_OK=TRUE
      for (j in 1:(length(Z))){
        model=try(survival::coxph(surv_perm~datas_perm[,j]))
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
    #Computation of score is there is covariates
    datas=cbind(covariates,counts_pathway)
    if (dim(covariates)==NULL){
      size_covariates=1
    }else{
      size_covariates=dim(covariates)[2]
    }

    for (i in 1:(dim(datas)[2]-(size_covariates))){
      model=try(survival::coxph(surv~datas[,i+(size_covariates)]+datas[,2:(size_covariates)]))
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
      datas=cbind(covariates,counts_pathway)
      datas=as.data.frame(datas)
    }

    #Computation of the epsilon matrix

    Z_matrix<- matrix(nrow=(length(Z)), ncol=nperm)
    Z_matrix[,1]=Z
    i=2
    while(i<=nperm){
      perm=sample(length(surv))
      surv_perm=surv[perm]
      covariates_perm=covariates[perm,]
      datas_perm=cbind(covariates_perm,counts_pathway)
      perm_OK=TRUE
      for (j in 1:(length(Z))){
        model=try(survival::coxph(surv_perm~datas[,j+(size_covariates)]+datas_perm[,2:(size_covariates)], data = datas_perm))
        if (length(model)>10){
          Z_matrix[j,i]=model$coefficients[1]/summary(model)$coefficients[1,3]
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

  }

  epsilon=cor(t(Z_matrix))

  scores_GBJ=list(test_stats=Z,cor_mat=epsilon)
  return(scores_GBJ)
}
