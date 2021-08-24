test_that("sGBJ_scores_consistency", {
  set.seed(1)
  n <- 10
  surv_data <- data.frame(Time = runif(n = n, min = 0, max = 100),
                          event = rbinom(n = n, size = 1, prob = 0.5))
  surv <- survival::Surv(time = surv_data$Time, event = surv_data$event)

  counts_pathway <- data.frame(P1 = rnorm(n = n),
                               P2 = rnorm(n = n))

  ls_results <- sGBJ::sGBJ_scores(surv,counts_pathway, nperm = 2)

  expect_equal(ls_results$test_stats, sGBJ::ls_test_results$test_stats)
  expect_equal(ls_results$cor_mat, sGBJ::ls_test_results$cor_mat)

})
