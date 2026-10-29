# Six documents; document 4 has no word with a vector.
outcome_text <- function() {
  emb <- load_wordpack(test_pack())
  ssd(list(c("dog", "cat"), c("cat", "duck"), c("dog", "swim"), c("zzz"), c("blue", "dog"), c("green", "cat")), emb)
}

test_that("values, a column name and a formula give the same rows", with_empty_cache({
  s <- outcome_text()
  df <- data.frame(score = c(1.5, NA, 3, 4, 5, 6), age = c(20, 30, 40, 50, 60, 70))
  by_values <- .outcome_rows(s, df$score, NULL, "pls")
  by_name <- .outcome_rows(s, "score", df, "pls")
  by_formula <- .outcome_rows(s, score ~ age, df, "pls")
  for (r in list(by_name, by_formula)) {
    expect_identical(r$y, by_values$y)
    expect_identical(r$keep, by_values$keep)
    expect_identical(r$x, by_values$x)
  }
  expect_equal(by_values$outcome_name, "y")
  expect_equal(by_name$outcome_name, "score")
  expect_equal(by_formula$outcome_name, "score")
  expect_identical(by_formula$formula, score ~ age)
  expect_equal(by_formula$data, data.frame(score = c(1.5, 3, 5, 6), age = c(20, 40, 60, 70)))
}))

test_that("rows with a missing outcome or no vector are dropped and counted", with_empty_cache({
  s <- outcome_text()
  r <- .outcome_rows(s, c(1, NA, 3, 4, Inf, 6), NULL, "regression")
  expect_equal(r$keep, c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
  expect_equal(r$y, c(1, 3, 6))
  expect_equal(c(r$n_raw, r$n_kept, r$n_dropped, r$n_no_vector), c(6, 3, 3, 1))
  # The kept rows' vectors are the text's: they are not rebuilt without the
  # dropped rows, and the SIF counts are those of all documents.
  expect_identical(r$x, s$x[c(1, 3, 5), ])
  expect_identical(s$sif_counts, .sif_counts(s$docs))

  labels <- .outcome_rows(s, factor(c("a", NA, "b", "a", "", "b")), NULL, "groups")
  expect_equal(labels$keep, c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
  expect_identical(labels$y, c("a", "b", "b"))

  m <- cbind(v = c(1, 2, NA, 4, 5, 6), a = c(1, 2, 3, 4, 5, NA))
  multi <- .outcome_rows(s, m, NULL, "multipls")
  expect_equal(multi$keep, c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_equal(multi$outcome_names, c("v", "a"))
  expect_equal(dim(multi$y), c(3L, 2L))
}))

test_that("each fit's outcome forms", {
  df <- data.frame(v = 1:6, a = c(2, 1, 2, 1, 2, 1), g = c(3, 1, 2, 3, 1, 2.5))
  expect_identical(.outcome_values(c(TRUE, FALSE, TRUE, NA, FALSE, TRUE), NULL, 6, "pls")$values,
                   c(1, 0, 1, NA, 0, 1))
  expect_equal(.outcome_values(c("v", "a"), df, 6, "multipls")$outcome_names, c("v", "a"))
  expect_equal(.outcome_values(df[c("v", "a")], NULL, 6, "multipls")$outcome_names, c("v", "a"))
  expect_equal(.outcome_values(unname(as.matrix(df[1:2])), NULL, 6, "multipls")$outcome_names, c("y1", "y2"))
  # Numbers as group labels read as Python writes them.
  expect_identical(.outcome_values("g", df, 6, "groups")$values, c("3", "1", "2", "3", "1", "2.5"))
  expect_identical(.outcome_values(df$g, NULL, 6, "groups_median")$values, df$g)
})

test_that("numeric labels and numbers in messages read as Python's str()", {
  # Python: str(v), and for labels str(int(v)) when v is whole.
  v <- c(5, 0.1 + 0.2, 1e-5, 1.5e-7, 1e16, 123456789012345678, 0.0001, 1e15, -2.5, 1 / 3, 1e20, 1234.5678, -0)
  expect_identical(.py_repr(v), c("5.0", "0.30000000000000004", "1e-05", "1.5e-07", "1e+16", "1.2345678901234568e+17",
                                  "0.0001", "1000000000000000.0", "-2.5", "0.3333333333333333", "1e+20", "1234.5678",
                                  "-0.0"))
  expect_identical(.number_labels(v), c("5", "0.30000000000000004", "1e-05", "1.5e-07", "10000000000000000",
                                        "123456789012345680", "0.0001", "1000000000000000", "-2.5",
                                        "0.3333333333333333", "100000000000000000000", "1234.5678", "0"))
})

test_that("an outcome a fit does not take raises", with_empty_cache({
  s <- outcome_text()
  df <- data.frame(v = 1:6, w = 6:1)
  expect_error(.outcome_rows(s, letters[1:6], NULL, "pls"),
               "fit_pls() requires a numeric outcome, but it holds categorical labels", fixed = TRUE)
  expect_error(.outcome_rows(s, letters[1:6], NULL, "groups_median"),
               "fit_groups(median_split = TRUE) requires a numeric outcome", fixed = TRUE)
  expect_error(.outcome_rows(s, as.matrix(df), NULL, "pls"),
               "fit_pls() takes one outcome, but it has 2 columns. Use fit_multipls()", fixed = TRUE)
  expect_error(.outcome_rows(s, df["v"], NULL, "regression"),
               "fit_regression() takes one outcome as a vector, but it is a 6 x 1 matrix", fixed = TRUE)
  expect_error(.outcome_rows(s, df$v, NULL, "multipls"),
               "fit_multipls() needs several outcomes (at least two columns)", fixed = TRUE)
  expect_error(.outcome_rows(s, v ~ w, df, "groups"), "fit_groups() takes no formula", fixed = TRUE)
  expect_error(.outcome_rows(s, ~ w, df, "pls"), "needs the outcome column on its left side", fixed = TRUE)
  expect_error(.outcome_rows(s, "v", NULL, "pls"),
               "fit_pls() got the column name \"v\" but no data. Pass data = df", fixed = TRUE)
  expect_error(.outcome_rows(s, "x", df, "pls"), "fit_pls(): data has no column \"x\".", fixed = TRUE)
  expect_error(.outcome_rows(s, 1:5, NULL, "pls"),
               "fit_pls(): the outcome has 5 values but the text has 6 documents.", fixed = TRUE)
}))

test_that("data must have one row per document and no reserved column", with_empty_cache({
  s <- outcome_text()
  expect_error(.outcome_rows(s, "v", data.frame(v = 1:5), "pls"),
               "nrow(data)=5 != s$n_docs=6. data must have one row per document.", fixed = TRUE)
  # The row check applies whatever form y takes.
  expect_error(.outcome_rows(s, 1:6, data.frame(v = 1:7), "pls"), "nrow(data)=7", fixed = TRUE)
  df <- data.frame(v = 1:6, .ssd_y = 1, .ssd_a = 2, check.names = FALSE)
  expect_error(.outcome_rows(s, "v", df, "regression"),
               "data columns \".ssd_a\", \".ssd_y\" start with '.ssd_', which fit_regression() reserves",
               fixed = TRUE)
}))

test_that("a formula column with a missing value in the kept rows raises", {
  data <- data.frame(age = c(1, NA, 3, NA), sex = c("f", "m", NA, "f"), essay = NA)
  expect_error(.check_nobs(y ~ age + sex, data),
               "The formula's columns have missing values in 3 of the 4 kept documents.", fixed = TRUE)
  # The column inside ssd() is a text, not a covariate.
  expect_silent(.check_nobs(y ~ ssd(essay), data))
  expect_error(.check_nobs(y ~ log(age) * ssd(essay), data), "missing values in 2 of the 4", fixed = TRUE)
  # A transform that makes NaN from a complete column counts; the outcome column does not.
  data <- data.frame(y = c(1, NA, 3, 4), cov = c(1, -1, 2, 3))
  expect_error(.check_nobs(y ~ log(cov) + ssd(), data), "missing values in 1 of the 4", fixed = TRUE)
  expect_silent(.check_nobs(y ~ cov + ssd(), data))
})
