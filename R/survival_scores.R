#' .survival_scores
#'
#' @description Compute the survival score
#'
#' @param surv a \code{\link[survival]{Surv}} object of length \code{n}
#' @param factor_matrix a \code{n x p} \code{data.frame} of the expression for
#' the particular gene set of interest being tested
#' @param covariates a matrix nxl of the covariates to adjust. Default is \code{NULL}
#'
#' @importFrom stats as.formula
#'
#' @return A list of length 3 with the updated factor_matrix (same as factor_matrix but removing columns for which survival model failed to converge), the Z matrix and the data used to fit survival model.
.survival_scores <- function(factor_matrix, covariates = NULL, surv){
  remove_Z <- NULL
  Z <- numeric(ncol(factor_matrix))

  if(is.null(covariates)){
    dat <- factor_matrix
  } else {
    dat <- cbind(covariates, factor_matrix)
    if (is.null(dim(covariates))){
      size_covariates <- 1
    }else{
      size_covariates <- ncol(covariates)
    }
  }

  remove_Z <- rep(FALSE, ncol(factor_matrix))
  for (i in 1:ncol(factor_matrix)){

    if(is.null(covariates)){
      model <- try(survival::coxph(surv ~ dat[,i]))
    } else {
      vecCovariates <- paste("dat[,", size_covariates, "]", collapse=" + ")

      vecPathway <- paste("dat[,", i+size_covariates,"]", sep = "")
      formX <- paste(c(vecPathway, vecCovariates), collapse = " + ")
      form <- as.formula(paste0("surv ~ ", formX))

      model <- try(survival::coxph(form, data = dat))
    }

    boolLengthModel <- length(model) > 10
    if(boolLengthModel){
      if(is.null(covariates)){
        Z[i] <- model$coefficients / summary(model)$coefficients[3]
      } else {
        Z[i] <- model$coefficients[1] / summary(model)$coefficients[1, 3]
      }
    }

    boolZna <- is.na(Z[i])
    boolZinf <- abs(Z[i]) == Inf
    boolRemove <- (!boolLengthModel) | boolZna | boolZinf

    if(boolRemove){
      remove_Z[i] <- TRUE
    }
  }

  Z <- Z[!remove_Z]
  updatedFactor_matrix <-  factor_matrix[, !remove_Z]

  return(list(Z = Z,
              updatedFactor_matrix = updatedFactor_matrix,
              dat = dat))
}
