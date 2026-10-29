# fit_regression(): the model paths, the tests, the whole-model test and the
# model generics (spec 03 "fit_regression()", 04 "Model generics").

# The same fit on the fastglmm path.
on_fastglmm <- function(...) {
  with_mocked_bindings(fit_regression(...), .regression_path = function(full, family) "fastglmm")
}

expect_same_tests <- function(base, fg, tol = 1e-6) {
  expect_identical(base$test$label, fg$test$label)
  expect_identical(base$test$name, fg$test$name)
  expect_equal(base$test$pvalue, fg$test$pvalue, tolerance = tol)
  expect_equal(base$test$statistic, fg$test$statistic, tolerance = tol)
  expect_equal(base$test$df, fg$test$df)
  expect_equal(base$test$df_resid, fg$test$df_resid)
  expect_equal(base$test$r2_delta, fg$test$r2_delta, tolerance = tol)
  expect_identical(base$terms$term, fg$terms$term)
  expect_identical(base$terms$name, fg$terms$name)
  expect_equal(base$terms[c("statistic", "df", "df_resid", "pvalue")], fg$terms[c("statistic", "df", "df_resid", "pvalue")],
               tolerance = tol)
  expect_identical(base$covariates$term, fg$covariates$term)
  expect_equal(base$covariates[-1L], fg$covariates[-1L], tolerance = tol)
}

test_that("OLS on lm gives fastglmm's compare, contrast_test and covariates", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- regression_text()
  base <- fit_regression(r$s, y ~ cov + group * ssd(), data = r$data, fixed_k = 5)
  fg <- on_fastglmm(r$s, y ~ cov + group * ssd(), data = r$data, fixed_k = 5)
  expect_s3_class(base$model, "lm")
  expect_s3_class(fg$model, "fastglmm")
  expect_same_tests(base, fg)
  # The whole-model test follows the model (OLS: F), not the path.
  expect_equal(base$model_test, fg$model_test, tolerance = 1e-6)
  expect_true(is.finite(fg$model_test$fvalue))
  expect_equal(base$stats$r2, fg$stats$r2, tolerance = 1e-8)
  expect_equal(base$units[["grouplow:ssd()"]]$beta, fg$units[["grouplow:ssd()"]]$beta, tolerance = 1e-6)
}))

for (fam in c("binomial", "poisson")) {
  test_that(paste("a", fam, "model on glm gives fastglmm's numbers"), with_empty_cache({
    skip_if_not_installed("fastglmm")
    r <- regression_text()
    f <- if (fam == "binomial") yb ~ cov + group * ssd() else yc ~ cov + group * ssd()
    base <- fit_regression(r$s, f, data = r$data, family = fam, fixed_k = 3)
    fg <- on_fastglmm(r$s, f, data = r$data, family = fam, fixed_k = 3)
    expect_s3_class(base$model, "glm")
    expect_identical(base$test$name, c("lrt", "wald_chisq", "wald_chisq"))
    expect_same_tests(base, fg)
    expect_equal(base$model_test, fg$model_test, tolerance = 1e-6)
    expect_identical(base$stats$backend, "PCA+GLMM")
    expect_true(is.na(base$stats$r2))
    # Wald intervals, not profiled ones.
    expect_equal(confint(base), confint.default(base$model)[c("(Intercept)", "cov", "grouplow", "groupmid"), ],
                 ignore_attr = TRUE)
  }))
}

test_that("random effects and families base R leaves out need fastglmm, and say why", with_empty_cache({
  r <- regression_text()
  fake <- need_pkg
  environment(fake) <- list2env(list(requireNamespace = function(...) FALSE), parent = environment(need_pkg))
  local_mocked_bindings(need_pkg = fake)
  expect_error(fit_regression(r$s, y ~ cov + (1 | g), data = r$data, fixed_k = 3),
               "fit_regression() needs the 'fastglmm' package (for a formula with random effects)", fixed = TRUE)
  expect_error(fit_regression(r$s, y ~ cov, data = transform(r$data, y = exp(y)), family = Gamma, fixed_k = 3),
               "(for family = Gamma)", fixed = TRUE)
  expect_error(fit_regression(r$s, y ~ cov, data = r$data, family = gaussian("log"), fixed_k = 3),
               "(for family = gaussian(link = \"log\"))", fixed = TRUE)
  expect_error(fit_regression(r$s, y ~ cov, data = r$data, fixed_k = 3), NA)
}))

test_that("a mixed model: LR test of the text, Wald rows, ranef per level", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- regression_text()
  res <- suppressWarnings(fit_regression(r$s, y ~ cov + group * ssd() + (1 | g), data = r$data, fixed_k = 3))
  expect_s3_class(res$model, "fastglmm")
  expect_identical(res$stats$backend, "PCA+GLMM")
  expect_identical(res$test$name, c("lrt", "wald_f", "wald_f"))
  expect_true(all(is.na(res$model_test[c("fvalue", "f_pvalue")])))
  expect_true(is.finite(res$model_test$llr))
  re <- ranef(res)
  expect_identical(names(re), "g")
  expect_identical(nrow(re$g), 24L)
  expect_identical(ranef(fit_regression(r$s, y ~ cov, data = r$data, fixed_k = 3)), list())
  expect_identical(names(coef(res)), c("(Intercept)", "cov", "grouplow", "groupmid"))
  expect_identical(df.residual(res), 240L - 4L - 9L)
}))

test_that("the whole-model test is the F test against y ~ 1 on OLS and the LR test otherwise", with_empty_cache({
  r <- regression_text()
  res <- fit_regression(r$s, y ~ cov + ssd(), data = r$data, fixed_k = 4)
  m1 <- res$model
  a <- anova(lm(.ssd_y ~ 1, data = m1$model), m1)
  expect_equal(res$model_test$f_pvalue, a[2L, "Pr(>F)"], tolerance = 1e-12)
  expect_equal(res$model_test$fvalue, a[2L, "F"], tolerance = 1e-12)
  expect_true(all(is.na(res$model_test[c("llr", "llr_pvalue")])))
  # Every label a text label: the all-text test itself.
  res <- fit_regression(r$s, y ~ ssd(), data = r$data, fixed_k = 4)
  expect_identical(res$model_test$f_pvalue, res$test$pvalue[1L])

  res <- fit_regression(r$s, yb ~ cov + ssd(), data = r$data, family = binomial(), fixed_k = 4)
  m0 <- glm(.ssd_y ~ 1, family = binomial(), data = res$model$model)
  expect_equal(res$model_test$llr, 2 * (as.numeric(logLik(res$model)) - as.numeric(logLik(m0))), tolerance = 1e-10)
  expect_equal(res$model_test$llr_pvalue, pchisq(res$model_test$llr, 5, lower.tail = FALSE), tolerance = 1e-10)

  res <- fit_regression(r$s, y ~ 0 + cov + ssd(), data = r$data, fixed_k = 4)
  expect_true(all(is.na(res$model_test)))
}))

test_that("the model generics on the lm path are lm's without the text's rows", with_empty_cache({
  r <- regression_text()
  d <- r$data
  rownames(d) <- paste0("doc", seq_len(nrow(d)))
  d$y[c(3, 10)] <- NA
  res <- fit_regression(r$s, y ~ cov + group + ssd(), data = d, fixed_k = 4)
  m <- res$model
  cov_rows <- c("(Intercept)", "cov", "grouplow", "groupmid")
  expect_identical(coef(res), coef(m)[cov_rows])
  expect_identical(vcov(res), vcov(m)[cov_rows, cov_rows])
  expect_equal(confint(res), confint(m)[cov_rows, ])
  expect_equal(confint(res, "cov", level = 0.9), confint(m, "cov", level = 0.9))
  expect_identical(nobs(res), 238L)
  expect_identical(df.residual(res), df.residual(m))
  expect_identical(logLik(res), logLik(m))
  expect_identical(AIC(res), AIC(m))
  expect_identical(BIC(res), BIC(m))
  expect_identical(names(residuals(res)), rownames(d)[-c(3, 10)])
  expect_equal(unname(residuals(res)), unname(residuals(m)))
  expect_equal(unname(fitted(res)), unname(fitted(m)))
  expect_identical(res$covariates$term, cov_rows)
  expect_equal(res$covariates$p_value, unname(summary(m)$coefficients[cov_rows, 4]))
  expect_equal(res$covariates$ci_low, unname(confint(m)[cov_rows, 1]))
  expect_identical(res$stats$n_kept, 238L)
  expect_identical(res$stats$n_dropped, 2L)
  expect_identical(res$texts[[1]]$keep, !is.na(d$y))
}))

test_that("covariates are empty without a formula and with only the intercept", with_empty_cache({
  r <- regression_text()
  expect_identical(nrow(fit_regression(r$s, "y", data = r$data, fixed_k = 3)$covariates), 0L)
  expect_identical(nrow(fit_regression(r$s, y ~ ssd(), data = r$data, fixed_k = 3)$covariates), 0L)
  res <- fit_regression(r$s, r$data$y, fixed_k = 3)
  expect_identical(res$outcome_name, "y")
  expect_null(res$formula)
  expect_null(res$fit_args$formula)
  expect_identical(names(residuals(res)), as.character(1:240))
}))

test_that("the result records its settings and the text's PCA", with_empty_cache({
  r <- regression_text()
  res <- fit_regression(r$s, y ~ cov + ssd(), data = r$data, fixed_k = 7)
  expect_s3_class(res, c("ssd_pca_ols", "ssd_result"))
  expect_identical(res$fit_args, list(formula = "y ~ cov + ssd()", family = "gaussian", link = "identity",
                                      fixed_k = 7, k_min = 2, k_max = 120, k_step = 2, ssd_gradients = 1,
                                      sweep_args = NULL))
  expect_identical(res$fit_info$pca_k_source, "fixed")
  expect_identical(res$fit_info$formula, "y ~ cov + ssd()")
  expect_identical(res$texts[[1]]$k, 7L)
  expect_identical(dim(res$texts[[1]]$pca_components), c(7L, 16L))
  expect_null(res$texts[[1]]$sweep)
  expect_null(res$texts[[1]]$name)
  expect_identical(fit_regression(r$s, y ~ ssd(essay), data = r$data, fixed_k = 3)$texts[[1]]$name, "essay")
  expect_equal(res$alignment_scores, drop(res$texts[[1]]$x %*% res$gradient))
  expect_equal(res$pca_weights, unname(coef(res$model)[paste0(".ssd_pca_", 1:7)]) / res$stats$y_std)
  expect_identical(res$stats$backend, "PCA+OLS")
}))

test_that("bad arguments raise", with_empty_cache({
  r <- regression_text()
  expect_error(fit_regression(r$s, "y", data = r$data, ssd_gradients = 1.5), "ssd_gradients must be a whole number")
  expect_error(fit_regression(r$s, "y", data = r$data, sweep_args = list(save_tables = TRUE)),
               "sweep_args takes only the keys", fixed = TRUE)
  expect_error(fit_regression(y ~ ssd(essay), "y", data = r$data), "s is a formula, so leave out y", fixed = TRUE)
  expect_error(fit_regression(y ~ ssd(essay), r$data),
               "fit_regression(): pass data = df by name: fit_regression(y ~ ssd(essay), data = df).", fixed = TRUE)
  expect_error(fit_regression(r$s, "y", data = r$data, family = "nonesuch"), "is not a family function")
  expect_error(fit_regression(r$s, "group", data = r$data), "requires a numeric outcome")
  expect_error(fit_regression(r$s, y ~ age + ssd(), data = r$data), "the formula uses 'age', which data does not have")
  d <- r$data
  d$cov[5] <- NA
  expect_error(fit_regression(r$s, y ~ cov, data = d, fixed_k = 3), "missing values in 1 of the 240 kept documents")
  d <- r$data
  d$grp <- ifelse(seq_len(nrow(d)) <= 2L, "rare", "common")
  expect_error(fit_regression(r$s, y ~ grp * ssd(), data = d, fixed_k = 3), "cannot all be estimated", fixed = TRUE)
}))

test_that("sweep_args reach the sweep, and the sweep picks the best joint score", with_empty_cache({
  r <- regression_text()
  res <- fit_regression(r$s, "y", data = r$data, k_min = 4, k_max = 16, k_step = 4,
                        sweep_args = list(cluster_topn = 30, auck_radius = 1))
  sw <- res$texts[[1]]$sweep
  expect_identical(names(sw), .SWEEP_COLUMNS)
  expect_identical(sw$k, c(4L, 8L, 12L, 16L))
  expect_true(all(sw$total_size <= 60L))
  expect_identical(res$texts[[1]]$k, sw$k[which.max(sw$joint_score)])
  expect_identical(res$fit_info$best_k, res$texts[[1]]$k)
  expect_identical(res$fit_info$pca_k_source, "sweep")
  expect_equal(sw$interp_auck[2], mean(sw$interp_resid_z[1:3]))
  expect_true(is.na(sw$beta_delta_1_minus_cos[1]))
}))

test_that("the sweep's scoring helpers follow Python's _sweep_math", {
  expect_equal(.compute_auck(c(1, NA, 3, 5), 1), c(1, 2, 4, 4))
  expect_equal(.zscore_ignore_nan(c(1, NA, 3)), c(-1, NA, 1))
  expect_equal(.zscore_ignore_nan(c(2, 2)), c(0, 0))
  d <- .detrend_by_variance(c(10, 20, 40, NA), c(1, 2, 3, 4))
  expect_equal(d$hat[1:3], unname(fitted(lm(c(1, 2, 3) ~ log(c(10, 20, 40))))))
  expect_true(is.na(d$hat[4]))
  expect_true(all(is.na(.detrend_by_variance(c(10, 20), c(1, 2))$hat)))
  ov <- .overall_interpretability(list(list(size = 3, coherence = 0.5, centroid_cos_beta = -0.2),
                                       list(size = 1, coherence = 0.9, centroid_cos_beta = 0.6)))
  expect_equal(ov[["mean_coherence"]], 0.6)
  expect_equal(ov[["aggregate"]], 0.6 * 0.3)
  expect_identical(.cosine(c(0, 0), c(1, 0)), NaN)
})

test_that("a Gamma model goes to fastglmm and keeps the test fastglmm reports", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- regression_text()
  d <- transform(r$data, yg = exp(y / 2))
  res <- fit_regression(r$s, yg ~ cov + ssd(), data = d, family = Gamma(link = "log"), fixed_k = 3)
  expect_s3_class(res$model, "fastglmm")
  expect_identical(res$test$name, "f_test")
  expect_identical(res$terms$name, "wald_f")
  expect_identical(res$fit_args[c("family", "link")], list(family = "Gamma", link = "log"))
  expect_true(is.finite(res$model_test$llr))
  expect_true(is.na(res$stats$r2))
}))

test_that("verbose = TRUE reports the sweep's choice", with_empty_cache({
  r <- regression_text()
  expect_message(capture.output(res <- fit_regression(r$s, "y", data = r$data, k_min = 2, k_max = 6, verbose = TRUE)),
                 "[sweep] best K = ", fixed = TRUE)
}))

test_that("an offset or random term before - 1 fits as the 0 + form", with_empty_cache({
  r <- regression_text()
  d <- transform(r$data, t = 1 + seq_len(240) %% 3)
  a <- fit_regression(r$s, yc ~ cov + ssd() + offset(log(t)) - 1, data = d, family = poisson(), fixed_k = 3)
  b <- fit_regression(r$s, yc ~ 0 + cov + ssd() + offset(log(t)), data = d, family = poisson(), fixed_k = 3)
  expect_equal(a$test, b$test)
  expect_equal(a$terms, b$terms)
  expect_equal(a$model_test, b$model_test)
  skip_if_not_installed("fastglmm")
  a <- suppressWarnings(fit_regression(r$s, y ~ cov + ssd() + (1 | g) - 1, data = d, fixed_k = 3))
  b <- suppressWarnings(fit_regression(r$s, y ~ 0 + cov + ssd() + (1 | g), data = d, fixed_k = 3))
  expect_equal(a$test, b$test)
  expect_equal(a$terms, b$terms)
}))

test_that("fastglmm's naming rules raise before fitting; R's own names stay hidden", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- regression_text()
  d <- r$data
  d$age.years <- d$cov
  expect_error(fit_regression(r$s, y ~ age.years + (1 | g), data = d, fixed_k = 3),
               "takes column names of letters, digits and underscores only, starting with a letter. Rename 'age.years'",
               fixed = TRUE)
  expect_identical(fit_regression(r$s, y ~ age.years, data = d, fixed_k = 3)$covariates$term,
                   c("(Intercept)", "age.years"))
  d$ssd__foo <- d$cov
  expect_error(fit_regression(r$s, y ~ ssd__foo + (1 | g), data = d, fixed_k = 3),
               "data columns 'ssd__foo' start with 'ssd__'", fixed = TRUE)
}))

test_that("a missing fastglmm stops before the sweep", with_empty_cache({
  r <- regression_text()
  fake <- need_pkg
  environment(fake) <- list2env(list(requireNamespace = function(...) FALSE), parent = environment(need_pkg))
  local_mocked_bindings(need_pkg = fake, .pca_sweep = function(...) stop("the sweep ran"))
  expect_error(fit_regression(r$s, y ~ cov + (1 | g), data = r$data), "needs the 'fastglmm' package", fixed = TRUE)
}))

test_that("the fastglmm path works without data", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- regression_text()
  res <- fit_regression(r$s, exp(r$data$y / 2), family = Gamma(link = "log"), fixed_k = 3)
  expect_s3_class(res$model, "fastglmm")
  expect_identical(res$terms$term, "ssd()")
}))

test_that("data columns starting with .ssd_ are reserved", with_empty_cache({
  r <- regression_text()
  d <- r$data
  d$.ssd_y <- 0
  expect_error(fit_regression(r$s, y ~ cov, data = d, fixed_k = 3),
               "data columns \".ssd_y\" start with '.ssd_', which fit_regression() reserves", fixed = TRUE)
  expect_error(fit_regression(list(essay = r$s), y ~ ssd(essay), data = d, fixed_k = 3),
               "start with '.ssd_'", fixed = TRUE)
}))

test_that("a mixed model's covariates take fastglmm's t values and p-values by name", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- regression_text()
  res <- suppressWarnings(fit_regression(r$s, y ~ cov + group + ssd() + (1 | g), data = r$data, fixed_k = 3))
  sm <- summary(res$model)$coefficients
  rows <- c("(Intercept)", "cov", "grouplow", "groupmid")
  expect_true("df" %in% colnames(sm))
  expect_identical(res$covariates$term, rows)
  expect_equal(res$covariates$statistic, unname(sm[rows, "t value"]))
  expect_equal(res$covariates$statistic, res$covariates$estimate / res$covariates$se)
  expect_equal(res$covariates$p_value, unname(sm[rows, "Pr(>|t|)"]))
}))

test_that("a coefficient table without one statistic or p-value column raises", {
  tab <- function(cols) matrix(0, 1L, length(cols), dimnames = list("cov", cols))
  expect_identical(.coef_column(tab(c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")), "^[tz] value$",
                                "statistic"), 4L)
  expect_error(.coef_column(tab(c("Estimate", "Std. Error", "stat", "p")), "^[tz] value$", "statistic"),
               "fastglmm's coefficient table has no single statistic column (its columns: Estimate, Std. Error, stat, p)",
               fixed = TRUE)
  expect_error(.coef_column(tab(c("Estimate", "Pr(>|t|)", "Pr(>|z|)")), "^Pr\\(", "p-value"),
               "no single p-value column", fixed = TRUE)
})

test_that("a result fitted inside a function does not keep the function's variables", with_empty_cache({
  r <- regression_text()
  size <- function(x) length(serialize(x, NULL))
  # `call` evaluated in a frame that also holds a 16 MB local.
  inside <- function(call) {
    big <- numeric(2e6)
    eval(call)
  }
  calls <- list(quote(fit_regression(r$s, y ~ cov + ssd(), data = r$data, fixed_k = 4)),
                quote(fit_regression(r$s, "y", data = r$data, fixed_k = 4)),
                quote(fit_regression(r$s, yb ~ cov, data = r$data, family = binomial(), fixed_k = 4)))
  if (requireNamespace("fastglmm", quietly = TRUE)) {
    calls <- c(calls, quote(suppressWarnings(fit_regression(r$s, y ~ cov + (1 | g), data = r$data, fixed_k = 3))))
  }
  if (requireNamespace("plskit", quietly = TRUE)) {
    calls <- c(calls, quote(fit_pls(r$s, y ~ cov, data = r$data, k = 2, test = FALSE)))
  }
  for (call in calls) {
    expect_lt(size(inside(call)), size(eval(call)) + 1e5, label = deparse1(call))
  }
}))
