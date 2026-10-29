# suggest_lexicon() and evaluate_lexicon() (ssdiff/utils/lexicon.py) against Python's
# results on a small text; the parity text's cases are in test-parity.R.

small_ref <- function() jsonlite::read_json(test_path("fixtures", "lexicon", "python_small.json"))

# The reference file holds Python's results for a 12-document token-list
# text on the fixture pack (SSD(docs, embeddings=Embeddings.load(fixture))):
# quantile bins with ties, 3 and 2 groups (Yates' correction), seeds with and
# without an outcome.
test_that("ranking, group statistics and seed rows match Python's on a small text", with_empty_cache({
  ref <- small_ref()
  s <- ssd(lapply(ref$docs, unlist), load_wordpack(test_pack()), use_full_doc = FALSE)
  y <- unlist(ref$y)
  g3 <- unlist(ref$g3)
  g2 <- unlist(ref$g2)
  expect_lexicon_equal(suggest_lexicon(s, y, min_docs = 2, top_k = 50, n_bins = 3), ref$continuous)
  expect_lexicon_equal(suggest_lexicon(s, g3, min_docs = 2, var_type = "categorical"), ref$categorical3)
  expect_lexicon_equal(suggest_lexicon(s, seeds = c("dog", "cat"), top_k = 3, min_docs = 1), ref$seeds)
  expect_lexicon_equal(suggest_lexicon(s, y, seeds = c("cat", "dog"), top_k = 2, min_docs = 1, n_bins = 2),
                       ref$seeds_y)
  expect_lexicon_equal(evaluate_lexicon(s, g2, c("dog", "cat", "swim", "dog"), var_type = "categorical"),
                       ref$evaluate2)
  expect_lexicon_equal(evaluate_lexicon(s, y, c("dog", "duck"), n_bins = 3), ref$evaluate_cont)
}))

test_that("y as values and as a column name give the same rows; missing y rows are dropped", with_empty_cache({
  ref <- small_ref()
  s <- ssd(lapply(ref$docs, unlist), load_wordpack(test_pack()), use_full_doc = FALSE)
  df <- data.frame(y = unlist(ref$y), g = unlist(ref$g3))
  expect_identical(suggest_lexicon(s, "y", data = df, min_docs = 2), suggest_lexicon(s, df$y, min_docs = 2))
  expect_identical(evaluate_lexicon(s, "g", "dog", data = df, var_type = "categorical"),
                   evaluate_lexicon(s, df$g, "dog", var_type = "categorical"))
  df$g[c(1, 5)] <- c(NA, "")
  ev <- evaluate_lexicon(s, "g", "dog", data = df, var_type = "categorical")
  expect_equal(ev$stats$n_docs, 10L)
  expect_equal(ev$suggestions$freq, 4L)       # "dog" is in documents 1, 2, 4, 7 and 11
}))

test_that("y = NULL: seeds or evaluate_lexicon only, outcome columns NA", with_empty_cache({
  ref <- small_ref()
  s <- ssd(lapply(ref$docs, unlist), load_wordpack(test_pack()), use_full_doc = FALSE)
  expect_error(suggest_lexicon(s), "Ranking needs an outcome. To list neighbours of words without one, pass seeds.",
               fixed = TRUE)
  ev <- evaluate_lexicon(s, NULL, c("duck", "dog", "swim", "cat"))
  expect_identical(ev$suggestions$token, c("cat", "dog", "duck", "swim"))   # by freq, then token
  expect_true(all(is.na(ev$suggestions[c("cov_bal", "corr", "pvalue", "direction", "rank")])))
  expect_equal(ev$stats$n_docs, 12L)
  expect_true(all(is.na(ev$summary[c("q1", "q4", "corr_any")])))
  expect_equal(ev$summary$docs_any, 11L)
}))

test_that("seeds are checked", with_empty_cache({
  s <- ssd(list(c("dog", "cat")), load_wordpack(test_pack()), use_full_doc = FALSE)
  expect_error(suggest_lexicon(s, seeds = c("dog", "zzz", "qq")),
               "These seeds are not in the embeddings: \"zzz\", \"qq\"", fixed = TRUE)
  expect_error(suggest_lexicon(s, seeds = character(0)), "seeds is empty", fixed = TRUE)
  expect_error(suggest_lexicon(s, 1:1, top_k = -1), "top_k must be >= 0, got -1", fixed = TRUE)
  expect_error(suggest_lexicon(s, 1, var_type = "ordinal"), "var_type must be 'continuous' or 'categorical'",
               fixed = TRUE)
}))

test_that("quantile bins follow numpy's percentile, ties included", {
  # np.percentile(1..8, [0, 25, 50, 75, 100]) = 1, 2.75, 4.5, 6.25, 8.
  expect_equal(.quantile_bins(1:8, 4), c(0, 0, 1, 1, 2, 2, 3, 3))
  # Edges 1, 1, 1, 1.75, 3 collapse to 1, 1.75, 3.
  expect_equal(.quantile_bins(c(1, 1, 1, 1, 2, 3), 4), c(0, 0, 0, 0, 1, 1))
  # One edge left: the median split.
  expect_equal(.quantile_bins(c(5, 5, 5), 4), c(0, 0, 0))
})

test_that("print() shows the stats and the first 20 rows", with_empty_cache({
  emb <- load_wordpack(test_pack())
  d <- golden_texts()
  s <- ssd(d$text, emb, use_full_doc = FALSE)
  lex <- suggest_lexicon(s, d$y)
  out <- capture.output(print(lex))
  expect_match(out[1], "25 tokens, 240 documents, continuous", fixed = TRUE)
  expect_length(grep("^ +[a-z]+ +[0-9]+ ", out), 20L)
  expect_identical(lexicon_tokens(lex), lex$suggestions$token)
}))

test_that("cov_bal averages as np.mean() does: numpy's pairwise double sum, not R's mean()", {
  # numpy: np.mean(np.array([1, 1, 3, 2]) / 7) and np.mean(np.arange(1, 21) / 7).
  v <- c(1, 1, 3, 2) / 7
  expect_identical(.np_mean(v), 0.24999999999999997)
  expect_identical(mean(v), 0.25)
  expect_identical(.np_mean((1:20) / 7), 1.5000000000000002)
  expect_identical(.np_mean(numeric(0)), 0)
})

test_that("the summary's medians are doubles; numeric group labels read as Python's", with_empty_cache({
  ref <- small_ref()
  s <- ssd(lapply(ref$docs, unlist), load_wordpack(test_pack()), use_full_doc = FALSE)
  sm <- evaluate_lexicon(s, NULL, c("dog", "cat"))$summary
  expect_type(sm$hits_median, "double")
  expect_type(sm$types_median, "double")
  ev <- evaluate_lexicon(s, rep(c(0.1 + 0.2, 2), 6), "dog", var_type = "categorical")
  expect_identical(names(ev$summary$group_cov[[1]]), c("0.30000000000000004", "2"))
}))
