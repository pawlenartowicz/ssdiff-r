# contrast() (ssdiff/results/contrast_result.py): the statistics on
# hand-built gradients and its accessors. Python's parity case is in test-parity.R.

test_that("contrast(): cos_g1_g2 and diff_norm on hand-built gradients", with_empty_cache({
  emb <- load_wordpack(test_pack())
  D <- emb$pack$dims
  g1 <- c(3, 0, rep(0, D - 2))                       # unit (1, 0, ...)
  g2 <- c(1, 1, rep(0, D - 2))                       # unit (1, 1, ...) / sqrt(2)
  cc <- contrast(g1, g2, emb, labels = c("a", "b"))
  expect_s3_class(cc, c("ssd_contrast", "ssd_result"))
  expect_equal(cc$stats$cos_g1_g2, 1 / sqrt(2))
  expect_equal(cc$stats$diff_norm, sqrt((1 - 1 / sqrt(2))^2 + 1 / 2))
  expect_equal(sqrt(sum(cc$beta^2)), 1)
  expect_equal(cc$beta[1:2], c(1 - 1 / sqrt(2), -1 / sqrt(2)) / cc$stats$diff_norm)
  expect_identical(cc$labels, c("a", "b"))
  # The scale of a gradient does not enter.
  expect_equal(contrast(10 * g1, 0.1 * g2, emb)$beta, cc$beta)
}))

test_that("contrast(): words and clusters work; documents raise", with_empty_cache({
  emb <- load_wordpack(test_pack())
  set.seed(1)
  cc <- contrast(stats::rnorm(emb$pack$dims), stats::rnorm(emb$pack$dims), emb)
  expect_identical(nrow(words(cc, n = 5)), 10L)
  expect_true(nrow(clusters(cc, side = "pos")) > 0L)
  expect_error(docs(cc), "A contrast has no documents.", fixed = TRUE)
  expect_error(doc_text(cc, 1), "A contrast has no documents.", fixed = TRUE)
  expect_error(snippets(cc), "A contrast has no documents, so it has no snippets.", fixed = TRUE)
  expect_error(cluster_snippets(cc), "A contrast has no documents, so it has no snippets.", fixed = TRUE)
  expect_error(words(cc, unit = 1), "A contrast has one unit", fixed = TRUE)
  expect_null(cc$texts[[1]]$name)
  expect_identical(cc$texts[[1]]$pack, .pack_ref(emb))
}))

test_that("contrast(): the checks", with_empty_cache({
  emb <- load_wordpack(test_pack())
  D <- emb$pack$dims
  g <- c(1, rep(0, D - 1))
  expect_error(contrast(g, g[-1], emb), sprintf("g1 and g2 must have the same length, got %d and %d.", D, D - 1),
               fixed = TRUE)
  expect_error(contrast(g, 2 * g, emb), "g1 and g2 are the same direction (or zero): the contrast is empty.",
               fixed = TRUE)
  expect_error(contrast(g, -g, emb, labels = "a"), "labels must hold two names, one per gradient.", fixed = TRUE)
  expect_error(contrast(g[1:3], g[1:3], emb), "the embeddings have", fixed = TRUE)
  expect_error(contrast(g, -g, list()), "emb must be embeddings from load_wordpack()", fixed = TRUE)
}))
