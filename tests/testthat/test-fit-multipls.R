# fit_multipls(): rotation modes, outcome groups and guards (spec 03
# "fit_multipls()").

# The golden text with four outcomes y1..y4 (two pairs of related columns).
multi_text <- function() {
  d <- golden_texts()
  list(s = ssd(d$text, load_wordpack(test_pack()), lexicon = d$lexicon),
       data = data.frame(y1 = d$y, y2 = d$y + d$cov, y3 = d$y2, y4 = d$y2 - 0.5 * d$y * d$cov))
}

ys <- c("y1", "y2", "y3", "y4")
fit <- function(...) fit_multipls(..., n_perm = 19, n_splits = 5)

test_that("every rotation gives orthonormal saliences, the same p-value and cov_explained", with_empty_cache({
  skip_if_not_installed("plskit")
  m <- multi_text()
  fits <- list(none = fit(m$s, ys, data = m$data, ssd_gradients = 2),
               words = fit(m$s, ys, data = m$data, ssd_gradients = 2, rotate = "words"),
               outcomes = fit(m$s, ys, data = m$data, ssd_gradients = 2, rotate = "outcomes"))
  xs <- .standardize(m$s$x)$x
  Ys <- .standardize(as.matrix(m$data))$x
  total <- sum(crossprod(xs, Ys)^2)
  for (f in fits) {
    expect_equal(crossprod(f$U_rot), diag(2), tolerance = 1e-8)
    expect_equal(crossprod(f$V_rot), diag(2), tolerance = 1e-8)
    expect_equal(f$stats$cov_explained, fits$none$stats$cov_explained, tolerance = 1e-10)
    expect_identical(f$stats$pvalue, fits$none$stats$pvalue)
    shares <- vapply(f$units, function(u) u$stats$cov_share, 0)
    expect_equal(sum(shares), 1)
    for (i in 1:2) {
      c_i <- sum(crossprod(Ys, xs %*% f$U_rot[, i])^2)
      expect_equal(shares[[i]] * f$stats$cov_explained, c_i / total, tolerance = 1e-8)
    }
    # The largest outcome loading of each dimension is positive.
    expect_true(all(apply(f$V_rot, 2L, function(v) v[which.max(abs(v))] > 0)))
    expect_identical(names(f$units), c("dim-1", "dim-2"))
    expect_equal(f$units[["dim-1"]]$beta, f$U_rot[, 1] / .standardize(m$s$x)$scale)
  }
  expect_identical(fits$words$stats$rotate, "words")
  expect_true(is.na(fits$none$stats$rotate))
}))

test_that("the default ssd_gradients is q, or one per rotation target row with groups", with_empty_cache({
  skip_if_not_installed("plskit")
  m <- multi_text()
  expect_identical(names(fit(m$s, ys, data = m$data)$units), paste0("dim-", 1:4))
  g <- fit(m$s, ys, data = m$data, rotate = "outcomes", outcome_groups = list(g1 = c("y1", "y2")))
  expect_identical(g$stats$ssd_gradients, 3L)
  expect_identical(g$outcome_groups, list(g1 = c("y1", "y2")))
  expect_null(g$fit_args$ssd_gradients)
  expect_identical(names(g$fit_args), c("ssd_gradients", "rotate", "outcome_groups", "n_splits", "n_perm",
                                        "random_state"))
}))

test_that("rotate = 'outcomes' at q loads each dimension on one outcome", with_empty_cache({
  skip_if_not_installed("plskit")
  m <- multi_text()
  V_rot <- fit(m$s, ys, data = m$data, rotate = "outcomes")$V_rot
  expect_equal(apply(abs(V_rot), 2L, max), rep(1, 4), tolerance = 1e-10)
  expect_equal(colSums(abs(V_rot)), rep(1, 4), tolerance = 1e-10)
  expect_identical(sort(apply(abs(V_rot), 2L, which.max)), 1:4)
}))

test_that("bad outcomes, rotate, groups and ssd_gradients raise with Python's messages", with_empty_cache({
  skip_if_not_installed("plskit")
  m <- multi_text()
  expect_error(fit(m$s, m$data$y1), "For one outcome use fit_pls()", fixed = TRUE)
  expect_error(fit(m$s, "y1", data = m$data), "For one outcome use fit_pls()", fixed = TRUE)
  d <- m$data
  d$y3 <- 3
  expect_error(fit(m$s, ys, data = d), "requires outcomes with nonzero variance; constant: 'y3'.", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, ssd_gradients = 0),
               "ssd_gradients must be a whole number in [1, 4] (min(q, n - 1, D) = min(4, 239, 16)); got 0",
               fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, ssd_gradients = 5), "[1, 4]", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, ssd_gradients = 3, rotate = "outcomes",
                   outcome_groups = list(a = c("y1", "y2"), b = c("y3", "y4"))),
               "needs ssd_gradients <= the number of rotation target rows (2: one per group", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, outcome_groups = list(a = c("y1", "y2"))),
               "outcome_groups requires rotate='outcomes'; got rotate=NULL", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, rotate = "words", outcome_groups = list(a = c("y1", "y2"))),
               "outcome_groups requires rotate='outcomes'; got rotate='words'", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, rotate = "outcomes", outcome_groups = list(a = c("y1", "nope"))),
               "outcome_groups[['a']]: unknown outcome 'nope'; outcomes are 'y1', 'y2', 'y3', 'y4'.", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, rotate = "outcomes", outcome_groups = list(a = "y1", b = c("y1", "y2"))),
               "outcome 'y1' is in two groups: 'a' and 'b'.", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, rotate = "outcomes", outcome_groups = list(a = character(0))),
               "outcome_groups[['a']] is empty.", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, rotate = "outcomes", outcome_groups = list(c("y1", "y2"))),
               "outcome_groups must be a list of character vectors with unique names", fixed = TRUE)
  expect_error(fit(m$s, ys, data = m$data, rotate = "varimax"),
               "rotate must be NULL, 'words' or 'outcomes'; got 'varimax'", fixed = TRUE)
  d <- m$data
  d$y5 <- d$y1 + d$y3
  expect_error(fit(m$s, c(ys, "y5"), data = d, ssd_gradients = 5), "k_used", fixed = TRUE)
}))
