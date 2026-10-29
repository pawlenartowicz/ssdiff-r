# The ssd() term of a fit_regression formula: placeholders, null formulas,
# term units and their Wald tests (spec 03, "The ssd() term", "Null formula").

null_of <- function(f) {
  f_ph <- .to_placeholders(f)
  f_ph[[2L]] <- quote(.ssd_y)
  labels <- .formula_terms(f_ph)$labels
  deparse1(.null_formula(f_ph, .text_labels(labels, ".ssd_text")))
}

test_that("null formulas drop every text label and keep the rest as written", {
  expect_identical(null_of(y ~ age + sex * ssd()), ".ssd_y ~ age + sex")
  expect_identical(null_of(y ~ sex:ssd()), ".ssd_y ~ 1")
  expect_identical(null_of(y ~ 0 + age * ssd()), ".ssd_y ~ 0 + age")
  expect_identical(null_of(y ~ age + ssd() + (1 | g)), ".ssd_y ~ age + (1 | g)")
  expect_identical(null_of(y ~ (1 + age | g) + ssd() + offset(log(w))), ".ssd_y ~ 1 + (1 + age | g) + offset(log(w))")
  expect_error(null_of(y ~ 0 + sex:ssd()), "has no fixed term left once the text terms are dropped", fixed = TRUE)
})

test_that("offsets and random terms are kept wherever they stand, - 1 and + 0 included", {
  expect_identical(null_of(y ~ cov + ssd() + offset(log(t)) - 1), ".ssd_y ~ 0 + cov + offset(log(t))")
  expect_identical(null_of(y ~ age + ssd() + (1 | g) - 1), ".ssd_y ~ 0 + age + (1 | g)")
  expect_identical(null_of(y ~ age + ssd() + 0 + (1 | g)), ".ssd_y ~ 0 + age + (1 | g)")
  expect_identical(null_of(y ~ 0 + age + offset(w) + ssd()), ".ssd_y ~ 0 + age + offset(w)")
  expect_identical(null_of(y ~ age - 1 + (1 | g) + ssd()), ".ssd_y ~ 0 + age + (1 | g)")
  expect_identical(null_of(y ~ offset(w) + ssd() - 1 + age), ".ssd_y ~ 0 + age + offset(w)")
  expect_identical(null_of(y ~ (cov + offset(log(t))) + ssd()), ".ssd_y ~ cov + offset(log(t))")
})

test_that("placeholders expand to the component columns, the rest of the formula untouched", {
  f <- .to_placeholders(y ~ age + sex * ssd(essay) + (1 | g))
  expect_identical(deparse1(f), "y ~ age + sex * .ssd_essay + (1 | g)")
  expect_identical(deparse1(.expand(f, list(.ssd_essay = c(".ssd_pca_1", ".ssd_pca_2")))),
                   "y ~ age + sex * (.ssd_pca_1 + .ssd_pca_2) + (1 | g)")
  expect_identical(deparse1(.with_text_term(y ~ age)), "y ~ age + ssd()")
  expect_identical(.text_names(y ~ ssd(a) + x:ssd() + ssd(a)), c("a", NA))
})

test_that("coefficient names group into term units, and each unit into its label", {
  comps <- list(.ssd_pca_1 = list("ssd(e)", 1L), .ssd_pca_2 = list("ssd(e)", 2L))
  g <- .term_groups(c("(Intercept)", "sexM", ".ssd_pca_2", ".ssd_pca_1", "sexM:.ssd_pca_1", "sexM:.ssd_pca_2"), comps)
  expect_identical(names(g), c("ssd(e)", "sexM:ssd(e)"))
  expect_identical(g[["ssd(e)"]]$pos, c(4L, 3L))
  expect_identical(.label_of("sexM:ssd(e)", c("ssd(e)", "sex:ssd(e)")), "sex:ssd(e)")
  expect_identical(.label_of("ssd(e):sexM", c("ssd(e)", "sex:ssd(e)")), "sex:ssd(e)")
  # The longest variable wins: ageGroupB comes from ageGroup, not age.
  expect_identical(.label_of("ageGroupB:ssd(e)", c("age:ssd(e)", "ageGroup:ssd(e)")), "ageGroup:ssd(e)")
  expect_identical(.label_of("sexM:ssd(e)", "ssd(e)"), NA_character_)
})

test_that("a formula without ssd() is the formula with + ssd()", with_empty_cache({
  r <- regression_text()
  a <- fit_regression(r$s, y ~ cov, data = r$data, fixed_k = 6)
  b <- fit_regression(r$s, y ~ cov + ssd(), data = r$data, fixed_k = 6)
  expect_equal(a$beta, b$beta)
  expect_identical(a$test, b$test)
  expect_identical(a$terms, b$terms)
  expect_identical(a$terms$term, "ssd()")
  expect_null(a$units)
}))

test_that("the units of group*ssd() add up to the gradients of group + group:ssd() per level", with_empty_cache({
  r <- regression_text()
  a <- fit_regression(r$s, y ~ group * ssd(), data = r$data, fixed_k = 6)
  b <- fit_regression(r$s, y ~ group + group:ssd(), data = r$data, fixed_k = 6)
  expect_identical(names(a$units), c("ssd()", "grouplow:ssd()", "groupmid:ssd()"))
  expect_identical(names(b$units), c("grouphigh:ssd()", "grouplow:ssd()", "groupmid:ssd()"))
  expect_equal(a$units[["ssd()"]]$beta, b$units[["grouphigh:ssd()"]]$beta, tolerance = 1e-8)
  for (lv in c("low", "mid")) {
    key <- sprintf("group%s:ssd()", lv)
    expect_equal(a$units[["ssd()"]]$beta + a$units[[key]]$beta, b$units[[key]]$beta, tolerance = 1e-8)
  }
  # Several units: no top-level gradient, every unit's y_hat is the model's fit.
  expect_null(a$beta)
  expect_equal(a$y_hat, unname(fitted(a$model)))
  expect_true(is.na(a$stats$beta_norm))
  expect_identical(a$units[["grouplow:ssd()"]]$test, "group:ssd()")
  expect_identical(a$units[["ssd()"]]$test, "ssd()")
  expect_equal(a$units[["grouplow:ssd()"]]$stats$beta_norm, sqrt(sum(a$units[["grouplow:ssd()"]]$beta^2)))
}))

test_that("every label of a text with several labels gets a Wald row; on OLS it equals the refit F", with_empty_cache({
  r <- regression_text()
  res <- fit_regression(r$s, y ~ group * ssd(), data = r$data, fixed_k = 6)
  expect_identical(res$test$label, c(NA, "ssd()", "group:ssd()"))
  expect_identical(res$test$name, c("f_test", "wald_f", "wald_f"))
  expect_identical(res$test$df[2:3], c(6, 12))
  expect_identical(res$terms$term, names(res$units))

  res <- fit_regression(r$s, y ~ cov * ssd(), data = r$data, fixed_k = 6)
  row <- res$test[res$test$label %in% "cov:ssd()", ]
  pcs <- paste0(".ssd_pca_", 1:6)
  m0 <- lm(reformulate(c("cov", pcs), ".ssd_y"), data = res$model$model)
  a <- anova(m0, res$model)
  expect_equal(row$statistic, a[2L, "F"], tolerance = 1e-10)
  expect_equal(row$pvalue, a[2L, "Pr(>F)"], tolerance = 1e-8)
  expect_identical(row$df_resid, as.numeric(a[2L, "Res.Df"]))

  # One label: only the all-text row.
  res <- fit_regression(r$s, y ~ cov + ssd(), data = r$data, fixed_k = 6)
  expect_identical(nrow(res$test), 1L)
}))

test_that("a lone moderated term takes its y_hat from the model", with_empty_cache({
  r <- regression_text()
  res <- fit_regression(r$s, y ~ cov:ssd(), data = r$data, fixed_k = 4)
  expect_identical(res$terms$term, "cov:ssd()")
  expect_equal(res$y_hat, unname(fitted(res$model)))
  plain <- fit_regression(r$s, y ~ cov + ssd(), data = r$data, fixed_k = 4)
  x <- plain$texts[[1]]$x
  pred <- mean(plain$y) + plain$stats$y_std * drop(sweep(x, 2L, colMeans(x)) %*% plain$beta)
  expect_equal(plain$y_hat, pred)
}))

test_that("misplaced or ambiguous ssd() terms raise", with_empty_cache({
  r <- regression_text()
  expect_error(fit_regression(r$s, y ~ cov + (1 + ssd() | g), data = r$data),
               "ssd(...) cannot be used inside a random-effect term: (1 + ssd() | g)", fixed = TRUE)
  expect_error(fit_regression(r$s, y ~ ssd(a) + ssd(b), data = r$data),
               "fit_regression() got one text, but the formula has ssd(a), ssd(b)", fixed = TRUE)
  expect_error(fit_regression(r$s, y ~ ssd() + ssd(a), data = r$data),
               "the formula has ssd(), ssd(a)", fixed = TRUE)
  expect_error(fit_regression(list(a = r$s), y ~ ssd() + ssd(a), data = r$data),
               "bare ssd() needs one text from ssd() as s", fixed = TRUE)
  expect_error(fit_regression(list(a = r$s, b = r$s), y ~ ssd(a) * ssd(b), data = r$data, fixed_k = 3),
               "The term 'ssd(a):ssd(b)' has two texts; two texts cannot interact in one term.", fixed = TRUE)
  expect_error(fit_regression(r$s, y ~ ssd("essay"), data = r$data), "one bare column name", fixed = TRUE)
  expect_error(fit_regression(list(y = r$s), y ~ ssd(y), data = r$data),
               "the text name 'y' is reserved for the outcome", fixed = TRUE)
}))
