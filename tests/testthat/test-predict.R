# predict() (spec 03, "predict()"): the training texts give the fitted y_hat;
# the other text form and fits with a formula raise.

test_that("predict() on the training texts equals y_hat: fit_regression()", with_empty_cache({
  r <- pca_ols_fit()
  expect_equal(predict(r, golden_texts()$text), r$y_hat, tolerance = 1e-10)
  expect_equal(predict(r, golden_texts()$text[c(5, 1)]), r$y_hat[c(5, 1)], tolerance = 1e-10)
}))

test_that("predict() on the training texts equals y_hat: fit_pls()", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_fit()
  expect_equal(predict(p, golden_texts()$text), p$y_hat, tolerance = 1e-10)
}))

test_that("predict() takes the fit's text form and gives NA without a vector", with_empty_cache({
  t <- results_text()
  r <- pca_ols_fit()
  expect_error(predict(r, as.list(t$d$text)), "This result was fit on strings; pass strings, not token lists.",
               fixed = TRUE)
  expect_identical(predict(r, c("", "no lexicon words here")), c(NA_real_, NA_real_))

  s <- ssd(t$s$docs, t$emb, lexicon = t$d$lexicon)
  tok <- fit_regression(s, t$d$y, fixed_k = 5)
  expect_equal(predict(tok, s$docs[1:10]), tok$y_hat[1:10], tolerance = 1e-10)
  expect_error(predict(tok, t$d$text[1:3]),
               "This result was fit on token lists; pass token lists, one list of tokens per document.", fixed = TRUE)
}))

test_that("predict() raises for a formula and for models other than Gaussian", with_empty_cache({
  msg <- "predict() works only for fits without a formula, and for fit_regression() only with the Gaussian family."
  t <- results_text()
  expect_error(predict(interaction_fit(), t$d$text[1:2]), msg, fixed = TRUE)
  expect_error(predict(fit_regression(t$s, y ~ ssd(), data = t$data, fixed_k = 3), t$d$text[1:2]), msg, fixed = TRUE)
  yb <- as.numeric(t$d$y2 > 0)
  expect_error(predict(fit_regression(t$s, yb, family = binomial(), fixed_k = 3), t$d$text[1:2]), msg, fixed = TRUE)
}))
