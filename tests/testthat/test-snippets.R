# Snippets (ssdiff/utils/snippets.py): windows, the fallbacks without a
# lexicon, and snippets, cluster snippets, clusters and docs against the
# page fixtures, which Python wrote from fits on the golden texts and the test
# pack. Python's ids there are 0-based.

# A snippet table from a page fixture, ids made 1-based.
py_snippets <- function(page, table, i) {
  py <- page_table(page, table, i)
  for (col in c("snippet_id", "doc_id", "start_token_idx", "end_token_idx", "start_sent_idx", "end_sent_idx",
                "cluster_id")) {
    py[[col]] <- as.integer(py[[col]]) + 1L
  }
  py
}

# Rows of one side and cluster whose cosines tie (to 1e-6) come in either
# order: R and Python break exact ties differently, and the small test pack
# makes such ties common. Within each tie the rows go in doc and token order.
tie_order <- function(tab, cos) {
  key <- paste(tab$side, tab$cluster_id)
  run <- cumsum(c(TRUE, abs(diff(cos)) > 1e-6 | key[-1L] != key[-length(key)]))
  `rownames<-`(tab[order(run, tab$doc_id, tab$start_token_idx), , drop = FALSE], NULL)
}

expect_snippets_match <- function(r_tab, py, label) {
  cols <- c("side", "doc_id", "seed", "start_token_idx", "end_token_idx", "start_sent_idx",
            "end_sent_idx", "text_window", "text_surface", "cluster_id")
  r_tab <- as.data.frame(r_tab)
  expect_identical(r_tab$snippet_id, py$snippet_id, label = label)
  expect_equal(r_tab$cosine, py$cosine, tolerance = 1e-6, label = label)
  expect_identical(tie_order(r_tab, py$cosine)[cols], tie_order(py, py$cosine)[cols], label = label)
}

# A cluster's top words in alphabetical order: words that tie on cos_centroid
# (the two words of a two-word cluster always do) come in either order.
sorted_words <- function(w) vapply(strsplit(w, ", ", fixed = TRUE), function(x) paste(sort(x), collapse = ", "), "")

# The rows of `sides` only, numbered from 1.
side_rows <- function(tab, sides) `rownames<-`(tab[tab$side %in% sides, , drop = FALSE], NULL)

# `sides`: the sides whose cluster tables are compared.
expect_page_unit <- function(res, page, i, unit, sides = c("pos", "neg")) {
  key <- page$units[[i]]$key
  expect_snippets_match(snippets(res, top_per_side = 10, unit = unit), py_snippets(page, "snippets", i), key)
  py <- page_table(page, "docs", i)
  d <- docs(res, unit = unit)
  expect_identical(d$doc_id, as.integer(py$doc_id) + 1L, label = key)
  for (col in setdiff(names(d), c("doc_id", "group"))) expect_equal(d[[col]], py[[col]], tolerance = 1e-6, label = col)
  if ("group" %in% names(d)) expect_identical(d$group, py$group)
  expect_snippets_match(side_rows(cluster_snippets(res, top_per_cluster = 5, unit = unit), sides),
                        side_rows(py_snippets(page, "cluster_snippets", i), sides), key)
  py <- side_rows(page_table(page, "clusters", i), sides)
  cl <- side_rows(clusters(res, unit = unit), sides)
  expect_identical(cl$cluster_id, as.integer(py$cluster_id) + 1L, label = key)
  expect_identical(data.frame(side = cl$side, size = cl$size, top_words = sorted_words(cl$top_words),
                              stringsAsFactors = FALSE),
                   data.frame(side = py$side, size = as.integer(py$size), top_words = sorted_words(py$top_words),
                              stringsAsFactors = FALSE), label = key)
  expect_equal(cl$coherence, py$coherence, tolerance = 1e-6)
}

test_that("groups: snippets, cluster snippets, clusters and docs match Python's page", with_empty_cache({
  g <- groups_fit()
  page <- page_fixture("groups")
  for (i in seq_along(page$units)) {
    u <- g$units[[page$units[[i]]$key]]
    expect_page_unit(g, page, i, c(u$g1, u$g2))
  }
}))

test_that("pca_ols: snippets, cluster snippets, clusters and docs match Python's page", with_empty_cache({
  expect_page_unit(pca_ols_fit(), page_fixture("pca_ols"), 1L, NULL)
}))

test_that("pls: snippets, cluster snippets, clusters and docs match Python's page", with_empty_cache({
  skip_if_not_installed("plskit")
  r <- results_text()
  p <- fit_pls(r$s, "y", data = r$data, k = 2, test = FALSE)
  expect_page_unit(p, page_fixture("pls"), 1L, NULL)
}))

test_that("multipls: snippets and docs of each dimension match Python's page", with_empty_cache({
  skip_if_not_installed("plskit")
  m <- multipls_fit()
  page <- page_fixture("multipls")
  for (i in seq_along(page$units)) {
    expect_page_unit(m, page, i, page$units[[i]]$key)
  }
}))

test_that("a window joins the previous or the next sentence when it reaches into it", {
  pre <- list(list(sents_surface = c("A b c.", "D e f.", "G h."), token_to_sent = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L)))
  emb_rows <- diag(8)
  tk <- list(tokens = c("a", "b", "c", "d", "e", "f", "g", "h"), doc = rep(1L, 8), pos = 1:8, w = rep(1, 8),
             wv = emb_rows, row = 1:8, sent = pre[[1]]$token_to_sent, n_tokens = 8L, n_sents = 3L, first = 0L,
             surf = pre[[1]]$sents_surface, essay = "A b c. D e f. G h.")
  occ <- .seed_occurrences(tk, c("b", "d", "e", "h"), window = 1L)
  expect_identical(occ$meta$seed, c("b", "d", "e", "h"))
  expect_identical(occ$meta$text_window, c("A b c.", "A b c. D e f.", "D e f.", "G h."))
  expect_identical(occ$meta$start_sent_idx, c(1L, 1L, 2L, 3L))
  expect_identical(occ$meta$end_sent_idx, c(1L, 2L, 2L, 3L))
  expect_identical(occ$meta$start_token_idx, c(1L, 3L, 4L, 7L))
  expect_identical(occ$meta$end_token_idx, c(3L, 5L, 6L, 9L))     # not clipped at the end, as in Python
  # The window's own word is left out: "b" is the mean of "a" and "c".
  expect_equal(occ$x[1, ], c(1, 0, 1, 0, 0, 0, 0, 0) / sqrt(2))
  # A window that only reaches the next sentence joins it.
  occ <- .seed_occurrences(tk, "c", window = 1L)
  expect_identical(occ$meta$text_window, "A b c. D e f.")
})

test_that("without a lexicon: one snippet per sentence, one cluster snippet per document", with_empty_cache({
  r <- results_text()
  s <- ssd(r$d$text, r$emb)
  g <- fit_groups(s, "group", data = r$data, n_perm = 10)
  sn <- snippets(g, top_per_side = 5, unit = c("g1", "g2"))
  expect_true(all(sn$seed == "<SENT>"))
  expect_identical(sn$start_sent_idx, sn$end_sent_idx)
  first <- sn[1L, ]
  expect_identical(first$text_window, s$pre_docs[[first$doc_id]]$sents_surface[first$start_sent_idx])
  cs <- cluster_snippets(g, top_per_cluster = 2, unit = c("g1", "g2"))
  expect_true(all(cs$seed == "<DOC>"))
  expect_identical(cs$text_window, cs$text_surface)
  expect_identical(cs$start_token_idx, rep(1L, nrow(cs)))
  expect_identical(cs$end_token_idx, lengths(s$docs)[cs$doc_id])
}))

test_that("min_cosine and top_per_side limit the rows; printing clips the windows", with_empty_cache({
  g <- groups_fit()
  all <- snippets(g, top_per_side = NULL, unit = c("g1", "g2"))
  high <- snippets(g, top_per_side = NULL, min_cosine = 0.3, unit = c("g1", "g2"))
  expect_true(all(high$cosine >= 0.3))
  expect_identical(nrow(high), sum(all$cosine >= 0.3))
  expect_identical(high$snippet_id, seq_len(nrow(high)))
  expect_identical(snippets(g, top_per_side = 2, unit = c("g1", "g2"))$cosine,
                   c(all$cosine[all$side == "pos"][1:2], all$cosine[all$side == "neg"][1:2]))
  out <- capture.output(print(snippets(g, top_per_side = 1, unit = c("g1", "g2"))))
  expect_true(any(grepl("…", out)))
  expect_true(nchar(all$text_surface[1L]) > 40L)
}))
