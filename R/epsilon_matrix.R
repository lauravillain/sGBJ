#' .epsilon_matrix
#'
#' @description Compute the epsilon matrix by permutation for the sGBJ_scores() function.
#'
#' @param Z The score vector returned by .survival_scores() function.
#' @param surv a surv object of size n
#' @param counts_pathway a data frame of the counts for the particular pathway of interest of size nxp
#' @param covariates a matrix nxl of the covariates to adjust (default=NULL)
#' @param nperm number of permutations to perform to estimate the matrix epsilon (default=300)
#' @param datas Data used to fit survival model returned by .survival_scores() function.
#'
#' @return The epsilon matrix.
.epsilon_matrix <- function(Z, nperm, surv, counts_pathway, covariates = NULL, datas){
  
  if(!is.null(covariates)){
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
    perm=sample(length(surv))
    surv_perm=surv[perm]
    if(is.null(covariates)){
      datas_perm=counts_pathway
    } else {
      covariates_perm=covariates[perm,]
      datas_perm=cbind(covariates_perm,counts_pathway)
    }
    for (j in 1:(length(Z))){
      if(is.null(covariates)){
        model=try(survival::coxph(surv_perm~datas_perm[,j]))
      } else {
        vecCovariates <- colnames(datas_perm)[1:size_covariates]
        vecPathway <- paste("datas_perm[,",j+size_covariates,"]",sep="")
        formX <- paste(c(vecPathway, vecCovariates), collapse = " + ")
        form <- as.formula(paste0("surv ~ ", formX))
        model=try(survival::coxph(form, data = datas_perm))
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