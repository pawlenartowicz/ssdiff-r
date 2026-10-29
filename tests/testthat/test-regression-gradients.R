# fit_regression(ssd_gradients = N): N rotated gradients from the model's
# own component weights (spec 03 "Gradients (ssd_gradients)"), after
# Python's test_fit_regression_gradients.py.

span_projector <- function(a) {
  q <- qr.Q(qr(a))
  q %*% t(q)
}

test_that("ssd_gradients = 1 is the standard fit", with_empty_cache({
  r <- regression_text()
  one <- fit_regression(r$s, "y", data = r$data, fixed_k = 10)
  expect_null(one$units)
  expect_null(one$pls_info)
  expect_identical(one$fit_info$ssd_gradients, 1L)
  expect_true(is.na(one$fit_info$rotate))
  expect_identical(one$fit_args$ssd_gradients, 1)
}))

test_that("ssd_gradients = 3 keeps the standard fit; the block pieces sum to beta and the units span them",
          with_empty_cache({
  skip_if_not_installed("plskit")
  r <- regression_text()
  one <- fit_regression(r$s, "y", data = r$data, fixed_k = 10)
  three <- fit_regression(r$s, "y", data = r$data, fixed_k = 10, ssd_gradients = 3)
  expect_identical(three$beta, one$beta)
  expect_identical(three$stats, one$stats)
  expect_identical(three$test, one$test)
  expect_identical(three$y_hat, one$y_hat)
  expect_identical(c(three$fit_info$ssd_gradients, three$fit_info$rotate), c("3", "words"))
  expect_identical(names(three$units), c("dim-1", "dim-2", "dim-3"))
  expect_identical(three$pls_info$ssd_gradients, 3L)
  expect_true(is.na(three$pls_info$random_state))

  # K = 10 in three consecutive blocks: 4 / 3 / 3, highest variance first.
  comps <- three$texts[[1]]$pca_components
  w <- three$pca_weights
  pieces <- sapply(list(1:4, 5:7, 8:10), function(b) drop(crossprod(comps[b, ], w[b])))
  x_scale <- .standardize(three$texts[[1]]$x)$scale
  expect_equal(rowSums(pieces) / x_scale, three$beta, tolerance = 1e-12)
  W <- sweep(pieces, 2L, sqrt(colSums(pieces^2)), "/")
  expect_equal(crossprod(W), diag(3), tolerance = 1e-12)
  B <- sapply(three$units, function(u) u$beta * x_scale)
  expect_equal(span_projector(B), span_projector(W), tolerance = 1e-8)

  corr <- vapply(three$units, function(u) u$stats$score_corr, 0)
  expect_true(all(corr > 0))
  expect_identical(corr, sort(corr, decreasing = TRUE))
  expect_equal(sum(vapply(three$units, function(u) u$stats$cov_share, 0)), 1)
}))

test_that("gradients with a formula or binomial() leave the standard fit", with_empty_cache({
  skip_if_not_installed("plskit")
  r <- regression_text()
  for (a in list(list(y ~ cov), list(yb ~ cov, family = binomial()))) {
    plain <- do.call(fit_regression, c(list(r$s), a, list(data = r$data, fixed_k = 6)))
    res <- do.call(fit_regression, c(list(r$s), a, list(data = r$data, fixed_k = 6, ssd_gradients = 2)))
    expect_identical(names(res$units), c("dim-1", "dim-2"))
    expect_identical(res$stats$pvalue, plain$stats$pvalue)
    expect_identical(res$beta, plain$beta)
  }
}))

test_that("ssd_gradients guards raise with Python's messages", with_empty_cache({
  skip_if_not_installed("plskit")
  r <- regression_text()
  expect_error(fit_regression(r$s, "y", data = r$data, fixed_k = 4, ssd_gradients = 5),
               "ssd_gradients=5 is more than the K=4 PCA components of the fit; use ssd_gradients <= K.", fixed = TRUE)
  expect_error(fit_regression(r$s, "y", data = r$data, fixed_k = 4, ssd_gradients = 0), "whole number >= 1",
               fixed = TRUE)
  expect_error(fit_regression(r$s, "y", data = r$data, fixed_k = 4, ssd_gradients = 1.5), "whole number >= 1",
               fixed = TRUE)
  for (f in list(y ~ group * ssd(), y ~ cov:ssd())) {
    expect_error(fit_regression(r$s, f, data = r$data, fixed_k = 4, ssd_gradients = 2),
                 "needs the text as one plain ssd() term", fixed = TRUE)
  }
  expect_error(fit_regression(list(a = r$s, b = r$s), y ~ ssd(a) + ssd(b), data = r$data, fixed_k = 4,
                              ssd_gradients = 2),
               "needs the text as one plain ssd() term, but the formula is 'y ~ ssd(a) + ssd(b)'", fixed = TRUE)
}))

test_that("a block with an aliased or zero piece raises", {
  p <- list(components = diag(4))
  expect_error(.regression_gradients(p, c(1, NA, 1, 1), 0, 2),
               "ssd_gradients: gradient 1 (components 1-2) has a coefficient the model could not estimate (aliased).",
               fixed = TRUE)
  expect_error(.regression_gradients(p, c(1, 1, 1, 0), 0, 3),
               "ssd_gradients: gradient 3 (component 4) is zero: the model gives those components no weight.",
               fixed = TRUE)
})
