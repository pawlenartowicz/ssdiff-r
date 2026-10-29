# print() and summary() (spec 04, "print() and summary()").

test_that("summary() of a regression: covariates from the model generics, terms, the F test", with_empty_cache({
  r <- interaction_fit()
  out <- capture.output(s <- summary(r))
  expect_identical(names(s), c("coefficients", "terms", "model_test"))
  expect_identical(s$coefficients[, "Estimate"], coef(r))
  expect_equal(s$coefficients[, "Std. Error"], sqrt(diag(vcov(r))))
  expect_identical(colnames(s$coefficients), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  expect_identical(s$terms, r$terms)
  expect_identical(out[1L], "Formula: y ~ cov + group * ssd()")
  expect_true(any(grepl("^grouplow:ssd\\(\\) +10 ", out)))
  expect_true(any(grepl("^all text terms ", out)))
  f <- summary(r$model)$fstatistic
  expect_true(any(startsWith(out, sprintf("F-statistic: %s on %d and %d DF", formatC(f[[1]], digits = 4), f[[2]], f[[3]]))))
}))

test_that("summary() shows the formula as Python does, and the LR test for other models", with_empty_cache({
  out <- capture.output(summary(pca_ols_fit()))
  expect_identical(out[1L], "Formula: y ~ ssd()")
  t <- results_text()
  yb <- as.numeric(t$d$y2 > 0)
  d <- cbind(t$data, yb = yb)
  out <- capture.output(s <- summary(fit_regression(t$s, yb ~ cov, data = d, family = binomial(), fixed_k = 3)))
  expect_identical(out[1L], "Formula: yb ~ cov + ssd()")
  expect_identical(colnames(s$coefficients)[3:4], c("z value", "Pr(>|z|)"))
  expect_true(any(startsWith(out, "Likelihood-ratio test against the intercept-only model:")))
  out <- capture.output(summary(fit_regression(t$s, y ~ 0 + cov + ssd(), data = t$data, fixed_k = 3)))
  expect_true(any(out == "No test of the whole model: the formula has no intercept."))
}))

test_that("summary() prints the same rows on the lm and the fastglmm paths", with_empty_cache({
  skip_if_not_installed("fastglmm")
  t <- results_text()
  lm_fit <- fit_regression(t$s, y ~ cov + group * ssd(), data = t$data, fixed_k = 5)
  fg_fit <- with_mocked_bindings(fit_regression(t$s, y ~ cov + group * ssd(), data = t$data, fixed_k = 5),
                                 .regression_path = function(full, family) "fastglmm")
  expect_s3_class(fg_fit$model, "fastglmm")
  lm_out <- capture.output(lm_s <- summary(lm_fit))
  fg_out <- capture.output(fg_s <- summary(fg_fit))
  # The statistic and p-value columns are read by name and the whole model
  # is tested by the model (OLS: F), so every line, the footer too, is the same.
  expect_identical(fg_out, lm_out)
  expect_identical(dimnames(fg_s$coefficients), dimnames(lm_s$coefficients))
  expect_equal(fg_s$coefficients, lm_s$coefficients, tolerance = 1e-6)
  expect_equal(fg_s$terms, lm_s$terms, tolerance = 1e-6)
  expect_equal(fg_s$model_test, lm_s$model_test, tolerance = 1e-6)
}))

test_that("print() names the accessors; summary() of the other results", with_empty_cache({
  g <- groups_fit()
  expect_output(print(g), "ssd_groups  G=3  n=240  omnibus p=")
  expect_output(print(g), "words(res)", fixed = TRUE)
  out <- capture.output(s <- summary(g))
  expect_identical(s$pairs, g$pairs)
  expect_true(any(grepl("Groups: g1 = high, g2 = low, g3 = mid", out, fixed = TRUE)))
  r <- pca_ols_fit()
  expect_output(print(r), "sweep_table(res)", fixed = TRUE)
  cc <- contrast(r$gradient, g$units[[1]]$gradient, load_wordpack(test_pack()), labels = c("a", "b"))
  expect_output(print(cc), "ssd_contrast  a − b  cos(g1, g2)=", fixed = TRUE)
  expect_identical(capture.output(s <- summary(cc))[1L], "Contrast [ssdiff]: a − b")
  expect_identical(s$stats, cc$stats)
}))

test_that("print() and summary() of the PLS results", with_empty_cache({
  skip_if_not_installed("plskit")
  p <- pls_fit()
  expect_output(print(p), "ssd_pls  r²=")
  expect_output(print(p), "\"dim-1\"  \"dim-2\"", fixed = TRUE)
  out <- capture.output(summary(p))
  expect_identical(out[1:2], c("PLS fit [ssdiff]", "Number of obs: 240"))
  expect_identical(out[5], "       test not run")
  m <- multipls_fit()
  expect_output(print(m), "outcomes(res)", fixed = TRUE)
  out <- capture.output(s <- summary(m))
  expect_identical(s$dims$dim, c("dim-1", "dim-2"))
  expect_equal(sum(s$dims$cov_share), 1)
}))
