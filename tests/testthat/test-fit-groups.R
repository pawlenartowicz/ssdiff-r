# The golden text (240 documents, all with a vector) and its data columns.
groups_text <- function() {
  d <- golden_texts()
  list(s = ssd(d$text, load_wordpack(test_pack()), lexicon = d$lexicon),
       data = data.frame(group = d$group, y = d$y, stringsAsFactors = FALSE))
}

test_that("correction names map to p.adjust() with Python's correct_pvalues() numbers", {
  # Python: ssdiff.backends.group.correct_pvalues on the same p-values.
  p <- list(c(0.01, 0.04, 0.03), c(0.02, 0.02, 0.001, 0.3, 0.049, 0.2), c(0.5, 0.6, 0.4), c(0.8, 0.01))
  python <- list(
    holm = list(c(0.03, 0.06, 0.06), c(0.1, 0.1, 0.006, 0.4, 0.14700000000000002, 0.4), c(1, 1, 1),
                c(0.8, 0.02)),
    bonferroni = list(c(0.03, 0.12, 0.09), c(0.12, 0.12, 0.006, 1, 0.29400000000000004, 1), c(1, 1, 1),
                      c(1, 0.02)),
    fdr_bh = list(c(0.03, 0.04, 0.04), c(0.04, 0.04, 0.006, 0.3, 0.07350000000000001, 0.24000000000000005),
                  c(0.6, 0.6, 0.6), c(0.8, 0.02)),
    none = p
  )
  for (m in names(python)) {
    for (i in seq_along(p)) {
      # Spec 03: the last bit for holm, bonferroni and none; fdr_bh within 1 ulp
      # (p.adjust takes n / i * p, Python p * n / i).
      if (m == "fdr_bh") {
        expect_equal(.correct_pvalues(p[[i]], m), python[[m]][[i]], tolerance = 1e-15, label = paste(m, i))
      } else {
        expect_identical(.correct_pvalues(p[[i]], m), python[[m]][[i]], label = paste(m, i))
      }
    }
  }
  expect_error(.correct_pvalues(0.05, "BH"),
               "Unknown correction 'BH'. Choose from ('holm', 'bonferroni', 'fdr_bh', 'none')", fixed = TRUE)
})

test_that("median split: ties are shared out with the shared generator, as in Python", {
  # Python: ssdiff.backends.group.median_split on the same values.
  y <- c(1, 2, 2, 2, 2, 3, 2, 5, 2, 0.5)
  expect_equal(.median_split(y, 2137), c("low", "high", "low", "low", "high", "high", "high", "high", "low", "low"))
  expect_equal(.median_split(y, 7), c("low", "low", "low", "high", "high", "high", "high", "high", "low", "low"))
  y2 <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4)
  expect_equal(.median_split(y2, 2137),
               c("low", "low", "low", "low", "high", "high", "low", "high", "high", "low", "high", "high",
                 "high", "high", "high", "low", "low", "low", "high", "low"))
  expect_error(.median_split(c(5, 5, 5)), "Cannot median-split: all y values are identical (5.0)", fixed = TRUE)
})

test_that("small groups are dropped with a warning; fewer than 2 groups raise", with_empty_cache({
  g <- groups_text()
  labels <- rep(c("a", "b", "c", "d"), c(110, 115, 15, 0))[seq_len(240)]
  labels[c(1, 2)] <- c(NA, "")
  expect_warning(res <- fit_groups(g$s, labels, n_perm = 20),
                 "Group 'c' dropped: only 15 docs (minimum 20)", fixed = TRUE)
  expect_equal(res$G, 2L)
  expect_equal(res$n_kept, 223L)
  expect_equal(res$texts[[1]]$keep, !is.na(labels) & labels %in% c("a", "b"))
  expect_equal(nrow(res$texts[[1]]$x), 223L)
  expect_equal(res$texts[[1]]$x, g$s$x[res$texts[[1]]$keep, ])
  expect_equal(res$groups, ifelse(labels[res$texts[[1]]$keep] == "a", "g1", "g2"))

  few <- rep(c("a", "b"), c(225, 15))
  expect_error(suppressWarnings(fit_groups(g$s, few, n_perm = 20)), "Need at least 2 groups after filtering, got 1",
               fixed = TRUE)
}))

test_that("groups are sorted by code point and pairs follow combn order", with_empty_cache({
  g <- groups_text()
  labels <- rep(c("zeta", "alpha", "Beta"), 80)
  res <- fit_groups(g$s, labels, n_perm = 20)
  expect_equal(res$group_labels, c(g1 = "Beta", g2 = "alpha", g3 = "zeta"))
  expect_equal(res$pairs$contrast, c("g1_g2", "g1_g3", "g2_g3"))
  expect_equal(names(res$units), res$pairs$contrast)
  expect_equal(res$units$g1_g3[c("g1", "g2")], list(g1 = "g1", g2 = "g3"))
  # Numbers as labels sort as Python's str(): "10" before "2".
  nums <- fit_groups(g$s, rep(c(2, 10), 120), n_perm = 20)
  expect_equal(nums$group_labels, c(g1 = "10", g2 = "2"))
}))

test_that("p-values are deterministic for a seed and ignore set.seed()", with_empty_cache({
  g <- groups_text()
  set.seed(1)
  a <- fit_groups(g$s, "group", data = g$data, n_perm = 200, random_state = 11)
  set.seed(2)
  b <- fit_groups(g$s, "group", data = g$data, n_perm = 200, random_state = 11)
  expect_identical(a$pairs, b$pairs)
  expect_identical(a$test, b$test)
  # Every p-value is a count over n_perm + 1.
  expect_equal(a$pairs$p_raw * 201, round(a$pairs$p_raw * 201))
}))

test_that("the result's fields, tables and units", with_empty_cache({
  g <- groups_text()
  res <- fit_groups(g$s, "group", data = g$data, n_perm = 50, correction = "fdr_bh", random_state = 3)
  expect_s3_class(res, c("ssd_groups", "ssd_result"), exact = TRUE)
  expect_named(res, c("G", "n_kept", "group_labels", "groups", "stats", "test", "pairs", "units", "texts",
                      "fit_args"))
  expect_named(res$stats, c("G", "n_kept", "n_perm", "correction", "random_state", "pvalue"))
  expect_named(res$test, c("name", "pvalue", "omnibus_T", "omnibus_p", "G", "n_kept", "n_perm", "correction",
                           "random_state"))
  expect_named(res$pairs, c("contrast", "g1", "g2", "T", "p_raw", "p_corrected", "cohens_d", "n_g1", "n_g2",
                            "contrast_norm"))
  expect_equal(res$pairs$p_corrected, p.adjust(res$pairs$p_raw, "BH"))
  expect_equal(res$test$omnibus_T, mean(res$pairs$T))
  expect_equal(res$fit_args, list(median_split = FALSE, n_perm = 50, correction = "fdr_bh", random_state = 3))
  expect_named(res$texts[[1]], c("name", "docs", "pre_docs", "tokenizer", "pack", "lexicon", "window", "sif_a",
                                 "use_full_doc", "x", "keep", "n_no_vector"))
  expect_null(res$texts[[1]]$name)
  expect_equal(res$texts[[1]]$n_no_vector, 0L)

  u <- res$units$g1_g2
  expect_named(u, c("g1", "g2", "beta", "gradient", "beta_norm", "alignment_scores"))
  x <- g$s$x
  high <- g$data$group == "high"
  low <- g$data$group == "low"
  expect_equal(u$beta, colMeans(x[high, ]) - colMeans(x[low, ]))
  expect_equal(u$gradient, u$beta / sqrt(sum(u$beta^2)))
  expect_equal(u$alignment_scores, drop(x[high | low, ] %*% u$gradient))
  # Cohen's d on the projections onto the unit contrast of the centroids, pooled SD with n - 1.
  cen <- function(m) colMeans(m) / sqrt(sum(colMeans(m)^2))
  proj <- drop(x %*% .unit_vector(cen(x[high, ]) - cen(x[low, ])))
  pooled <- sqrt((69 * var(proj[high]) + 79 * var(proj[low])) / 148)
  expect_equal(res$pairs$cohens_d[1], (mean(proj[high]) - mean(proj[low])) / pooled)
  expect_equal(res$pairs$T[1], 1 - sum(cen(x[high, ]) * cen(x[low, ])))
}))

test_that("fit_groups() checks its arguments", with_empty_cache({
  g <- groups_text()
  expect_error(fit_groups(g$s, "group", data = g$data, correction = "fdr"), "Unknown correction 'fdr'", fixed = TRUE)
  expect_error(fit_groups(g$s, "group", data = g$data, n_perm = 1.5), "n_perm must be a whole number")
  expect_error(fit_groups(g$s, "group", data = g$data, n_perm = Inf), "n_perm must be a whole number >= 0, got Inf",
               fixed = TRUE)
  expect_error(fit_groups(g$s, "group", data = g$data, random_state = -1), "seed must be a whole number")
  expect_error(fit_groups(g$s, "group", data = g$data, median_split = TRUE),
               "fit_groups(median_split = TRUE) requires a numeric outcome", fixed = TRUE)
}))
