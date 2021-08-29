#' .epsilon_matrix
#'
#' @description Compute the epsilon matrix by permutation for the \code{sGBJ_scores()} function.
#'
#' @param Z the score vector returned by \code{.survival_scores()} function.
#' @param surv a \code{\link[survival]{Surv}} object of length \code{n}
#' @param factor_matrix a \code{n x p} \code{data.frame} of the expression for the
#' particular gene set of interest being tested
#' @param covariates a \code{n x l} matrix of the covariates to adjust upon. Default is \code{NULL}
#' @param nperm number of permutations performed to estimate the \code{epsilon} matrix.
#' Default is \code{300}.
#' @param dat data used to fit survival model returned by \code{.survival_scores()} function.
#'
#' @return The epsilon matrix.
.epsilon_matrix <- function(Z, nperm, surv, factor_matrix, covariates = NULL, dat){

  if(!is.null(covariates)){
    if (is.null(dim(covariates))){
      size_covariates <- 1
    }else{
      size_covariates <- ncol(covariates)
    }
  }

  # build Z_matrix
  Z_matrix <- matrix(nrow = (length(Z)), ncol = nperm)
  Z_matrix[,1] <-  Z
  i <- 2

  # fill Z_matrix
  while(i <= nperm){
    perm_OK <- TRUE
    perm <- sample(length(surv))
    surv_perm <- surv[perm]
    if(is.null(covariates)){
      dat_perm <- factor_matrix
    } else {
      if (size_covariates>1){
        covariates_perm <- covariates[perm,]
      }else{
        covariates_perm <- unlist(covariates)[perm]
      }
      dat_perm=cbind(covariates_perm, factor_matrix)
    }
    for (j in 1:(length(Z))){
      if(is.null(covariates)){
        model=try(survival::coxph(surv_perm ~ dat_perm[, j]))
      } else {
        vecCovariates <- paste("dat_perm[,", size_covariates, "]", collapse=" + ")
        vecPathway <- paste("dat_perm[,", j + size_covariates,"]", sep="")
        formX <- paste(c(vecPathway, vecCovariates), collapse = " + ")
        form <- as.formula(paste0("surv_perm ~ ", formX))
        model <- try(survival::coxph(form, data = dat_perm))
      }

      boolLengthModel <- length(model)>10
      if(boolLengthModel){
        if(is.null(covariates)){
          Z_matrix[j, i] <- model$coefficients / summary(model)$coefficients[3]
        } else {
          Z_matrix[j, i] <- model$coefficients[1] / summary(model)$coefficients[1,3]
        }
      }

      boolZna <- !is.na(Z_matrix[j,i])
      boolZinf <- abs(Z_matrix[j,i]) != Inf

      perm_OK <- boolLengthModel & boolZna & boolZinf
    }
    i <- i + perm_OK
  }

  epsilon <- cor(t(Z_matrix))

  return(epsilon)
}
