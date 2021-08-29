#' Compute the sGBJ statistic and its p-value quantifying a gene set expression association
#' with survival
#'
#' @param surv a \code{\link[survival]{Surv}} object of length \code{n}
#' @param factor_matrix a \code{n x p} \code{data.frame} of the expression for the
#' particular gene set of interest being tested
#' @param covariates a \code{n x l} matrix of the covariates to adjust upon. Default is \code{NULL}
#' @param nperm number of permutations performed to estimate the \code{epsilon} matrix.
#' Default is \code{300}.
#'
#' @return The sGBJ statistic and its associated p-value associated
#' @examples
#' n <- 100
#' surv_data <- data.frame(Time = runif(n = n, min = 0, max = 100),
#'                         event = rbinom(n = n, size = 1, prob = 0.5))
#' surv <- survival::Surv(time = surv_data$Time, event = surv_data$event)
#'
#' factor_matrix <- data.frame(P1 = rnorm(n = n),
#'                              P2 = rnorm(n = n))
#'
#' sGBJ::sGBJ(surv,factor_matrix, nperm = 2)
#'
#' @export
#'
sGBJ <- function(surv, factor_matrix, covariates = NULL, nperm = 300){

  scores_GBJ <- sGBJ_scores(surv, factor_matrix, covariates, nperm)

  GBJOut <- GBJ::GBJ(test_stats = scores_GBJ$test_stats,
                     cor_mat = scores_GBJ$cor_mat)

  names(GBJOut) <- c("sGBJ_stat", "sGBJ_pvalue", "err_code")

  return(GBJOut)
}
