# fit_pls(): covariates, guards, the test switch and several gradients
# (spec 03 "fit_pls()").

# The golden text (240 documents, D = 16) and its data columns.
pls_text <- function() {
  d <- golden_texts()
  list(s = ssd(d$text, load_wordpack(test_pack()), lexicon = d$lexicon),
       data = data.frame(y = d$y, cov = d$cov, group = d$group, stringsAsFactors = FALSE))
}

FIT <- list(k = 2, n_splits = 10, random_state = 0)
fit <- function(...) do.call(fit_pls, c(list(...), FIT))

span_projector <- function(a) {
  q <- qr.Q(qr(a))
  q %*% t(q)
}

test_that("covariates come out by qr.resid, as lm() residuals column by column", with_empty_cache({
  p <- pls_text()
  x <- p$s$x
  res <- .pls_residualize(y ~ cov + group, p$data, p$data$y, x)
  expect_equal(res$y, unname(lm(y ~ cov + group, p$data)$residuals), tolerance = 1e-10)
  for (j in c(1, 7, 16)) {
    expect_equal(res$x[, j], unname(lm(x[, j] ~ cov + group, p$data)$residuals), tolerance = 1e-10)
  }
}))

test_that("a formula fits PLS on the residuals; r2 is the residualized y's", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  r <- fit(p$s, y ~ cov + group, data = p$data)
  res <- .pls_residualize(y ~ cov + group, p$data, p$data$y, p$s$x)
  manual <- with_mocked_bindings(fit(p$s, res$y), .x_rows = function(s, keep) res$x)
  expect_equal(r$stats$r2, manual$stats$r2, tolerance = 1e-12)
  expect_equal(r$beta, manual$beta, tolerance = 1e-12)
  expect_identical(r$fit_info$formula, "y ~ cov + group")
  expect_identical(r$fit_args$formula, "y ~ cov + group")
  expect_equal(r$residuals, manual$residuals, tolerance = 1e-10)
  # The raw y is described, the residualized y fitted.
  expect_equal(r$stats$y_mean, mean(p$data$y))
  expect_equal(r$y, p$data$y)
}))

test_that("values, a column name and a formula ~ 1 give the same fit", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  a <- fit(p$s, p$data$y)
  b <- fit(p$s, "y", data = p$data)
  f <- fit(p$s, y ~ 1, data = p$data)
  expect_identical(a$beta, b$beta)
  expect_equal(f$beta, a$beta, tolerance = 1e-10)
  expect_identical(c(a$outcome_name, b$outcome_name), c("y", "y"))
}))

test_that("formula, k, ssd_gradients, rotate and constant y raise with Python's messages", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  expect_error(fit_pls(p$s, y ~ cov + ssd(), data = p$data), "no ssd() term", fixed = TRUE)
  expect_error(fit_pls(p$s, y ~ cov + (1 | group), data = p$data),
               "fit_pls() takes fixed covariates only, but the formula has random effects", fixed = TRUE)
  p$data$y2 <- 2 * p$data$y + 1
  expect_error(fit_pls(p$s, y ~ y2, data = p$data),
               "The formula's covariates explain y fully, so no variance is left for the text.", fixed = TRUE)
  expect_error(fit_pls(p$s, rep(3, 240)),
               "fit_pls() requires y with nonzero variance (std=0). PLS1 is undefined when y is constant.", fixed = TRUE)
  expect_error(fit_pls(p$s, p$data$y, k = 0), "k must be >= 1, got 0", fixed = TRUE)
  expect_error(fit_pls(p$s, p$data$y, k = "two"), "k must be a whole number or \"auto\", got \"two\"", fixed = TRUE)
  expect_error(fit_pls(p$s, p$data$y, ssd_gradients = 0),
               "ssd_gradients must be a whole number in [1, 16] (min(n - 1, D) = min(239, 16)); got 0", fixed = TRUE)
  expect_error(fit_pls(p$s, p$data$y, ssd_gradients = 17), "[1, 16]", fixed = TRUE)
  expect_error(fit_pls(p$s, p$data$y, ssd_gradients = 2.5), "[1, 16]", fixed = TRUE)
  expect_error(fit_pls(p$s, p$data$y, ssd_gradients = 3, rotate = "outcomes"),
               "fit_pls() rotate must be 'words' or NULL; got 'outcomes'", fixed = TRUE)
}))

test_that("k = 'auto' with no component found raises", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  fake <- function(...) list(k_star = 0L)
  expect_error(with_mocked_bindings(fit_pls(p$s, p$data$y), pls1_find_k_optimal = fake, .package = "plskit"),
               "plskit found no component (k_star = 0)", fixed = TRUE)
}))

test_that("test = FALSE gives an NA p-value and the same fit", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  a <- fit(p$s, p$data$y)
  b <- fit(p$s, p$data$y, test = FALSE)
  expect_true(is.na(b$stats$pvalue))
  expect_true(all(is.na(unlist(b$test))))
  expect_false(is.na(a$stats$pvalue))
  expect_identical(a$test$pvalue, a$stats$pvalue)
  expect_identical(b$beta, a$beta)
  expect_null(a$find_k)
}))

test_that("ssd_gradients = 3 leaves the standard fit and adds three gradients spanning PLS's W", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  one <- fit(p$s, p$data$y)
  three <- fit(p$s, p$data$y, ssd_gradients = 3)
  expect_null(one$units)
  expect_null(one$pls_info)
  expect_identical(three$beta, one$beta)
  expect_identical(three$stats, one$stats)
  expect_identical(three$fit_info$n_components, one$fit_info$n_components)
  expect_identical(c(three$fit_info$ssd_gradients, three$fit_info$rotate), c("3", "words"))
  expect_identical(names(three$units), c("dim-1", "dim-2", "dim-3"))

  st <- .standardize(p$s$x)
  ys <- drop(.standardize(p$data$y)$x)
  W <- plskit::pls1_fit(st$x, ys, k = 3L, pre_standardized = TRUE, seed = 0)$W
  B <- sapply(three$units, function(u) u$beta * st$scale)
  expect_equal(span_projector(B), span_projector(W), tolerance = 1e-8)

  corr <- vapply(three$units, function(u) u$stats$score_corr, 0)
  expect_true(all(corr > 0))
  expect_identical(corr, sort(corr, decreasing = TRUE))
  expect_equal(sum(vapply(three$units, function(u) u$stats$cov_share, 0)), 1)
  expect_identical(three$pls_info$ssd_gradients, 3L)
  expect_identical(sort(three$pls_info$order[[1]]), 1:3)
}))

test_that("rotate = NULL keeps the NIPALS weights in order", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  res <- fit(p$s, p$data$y, ssd_gradients = 3, rotate = NULL)
  st <- .standardize(p$s$x)
  ys <- drop(.standardize(p$data$y)$x)
  W <- plskit::pls1_fit(st$x, ys, k = 3L, pre_standardized = TRUE, seed = 0)$W
  for (i in 1:3) expect_equal(res$units[[i]]$beta * st$scale, W[, i], tolerance = 1e-12)
  expect_equal(res$units[["dim-2"]]$stats$cov_share, 0, tolerance = 1e-12)
  expect_true(is.na(res$fit_info$rotate))
  expect_identical(res$pls_info$order[[1]], 1:3)
}))

test_that("the gradients do not depend on k: the main fit's weights or a new fit give the same", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  wide <- fit_pls(p$s, p$data$y, k = 3, ssd_gradients = 2, test = FALSE)    # W of the k = 3 fit
  narrow <- fit_pls(p$s, p$data$y, k = 1, ssd_gradients = 2, test = FALSE)  # W of a new k = 2 fit
  expect_equal(wide$units, narrow$units, tolerance = 1e-10)
  expect_equal(wide$pls_info, narrow$pls_info, tolerance = 1e-10)
}))

test_that("a formula with ssd_gradients keeps the standard fit", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  plain <- fit(p$s, y ~ cov, data = p$data)
  res <- fit(p$s, y ~ cov, data = p$data, ssd_gradients = 3)
  expect_identical(res$beta, plain$beta)
  expect_identical(res$stats$pvalue, plain$stats$pvalue)
  expect_identical(names(res$units), c("dim-1", "dim-2", "dim-3"))
}))

test_that("fit_args holds Python's keys", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_text()
  r <- fit_pls(p$s, "y", data = p$data, k = 1, n_splits = 5)
  expect_identical(names(r$fit_args),
                   c("formula", "k", "k_max", "ssd_gradients", "rotate", "n_splits", "random_state", "test"))
  expect_null(r$fit_args$formula)
  expect_s3_class(r, c("ssd_pls", "ssd_result"))
  expect_identical(r$texts[[1]]$keep, p$s$keep)
}))
