#' surv_calc_scores_stats
#'
#' @description An adaptation of \code{GBJ::calc_scores_stats()} to survival context.
#' Wrapper of sGBJ_scores() function.
#'
#' @param null_model An R cox model fitted with \code{survival::coxph()}.
#' @param factor_matrix An \code{n x p} matrix with each factor as one column. There should be no missing data.
#' @param nperm Number of permutations (default is 300)
#'
#' @return A list with the elements:
#' \item{test_stats}{The \code{p} score test statistics.}
#' \item{cor_mat}{The \code{p x p} matrix giving the pairwise correlation of every test statistic pairs.}
#'
#' @export
#'
#' @examples
#' n <- 100
#' surv_data <- data.frame(Time = runif(n = n, min = 0, max = 100),
#'                         event = rbinom(n = n, size = 1, prob = 0.5))
#' surv <- survival::Surv(time = surv_data$Time, event = surv_data$event)
#'
#' factor_matrix <- data.frame(P1 = rnorm(n = n),
#'                              P2 = rnorm(n = n))
#'
#' covariates <- data.frame(age = runif(n = n, 60, 90))
#'
#' null_model <- survival::coxph(surv ~ age, data = covariates, x = TRUE)
#' surv_reg_stats <- surv_calc_scores_stats(null_model = null_model,
#'                                          factor_matrix = factor_matrix,
#'                                          nperm = 2)#nperm = 300)
#'
#' GBJ::GBJ(test_stats=surv_reg_stats$test_stats, cor_mat=surv_reg_stats$cor_mat)
#'
surv_calc_scores_stats <- function(null_model,
                                   factor_matrix,
                                   nperm = 300){

  # check null_model
  if(class(null_model) != "coxph") {
    stop("null_model must be a 'coxph' object from survival::coxph()")
  }
  if(!any(names(null_model) == "x")) {
    stop("Please be sure to provide the 'x = TRUE' option to survival::coxph() function")
  }

  res <- sGBJ_scores(surv = null_model$y,
                     factor_matrix = factor_matrix,
                     nperm = nperm,
                     covariates = null_model$x)

  return(res)
}
