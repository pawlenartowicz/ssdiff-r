# The accessors (spec 04, "Accessors" and Tests): columns against the page
# fixtures, `unit`, term units, determinism, docs and the pack after readRDS().

expect_columns <- function(result, page, unit = NULL, i = 1L) {
  expect_identical(names(words(result, unit = unit)), page_columns(page, "words", i))
  expect_identical(names(clusters(result, unit = unit)), setdiff(page_columns(page, "clusters", i), "top_snippet"))
  expect_identical(names(cluster_words(result, unit = unit)), page_columns(page, "cluster_words", i))
  if (inherits(result, "ssd_contrast")) return(invisible())
  expect_identical(names(snippets(result, unit = unit)), page_columns(page, "snippets", i))
  expect_identical(names(cluster_snippets(result, top_per_cluster = 2, unit = unit)),
                   page_columns(page, "cluster_snippets", i))
  expect_identical(names(docs(result, unit = unit)), setdiff(page_columns(page, "docs", i), "text"))
}

test_that("the tables have Python's columns: groups, pca_ols, several terms, contrast", with_empty_cache({
  expect_columns(groups_fit(), page_fixture("groups"), unit = c("g1", "g2"))
  expect_columns(pca_ols_fit(), page_fixture("pca_ols"))
  expect_columns(interaction_fit(), page_fixture("regression"), unit = "ssd()")
  r <- pca_ols_fit()
  cc <- contrast(r$gradient, groups_fit()$units[[1]]$gradient, load_wordpack(test_pack()))
  expect_columns(cc, page_fixture("contrast"))
}))

test_that("the tables have Python's columns: pls and multipls", with_empty_cache({
  skip_if_not_installed("plskit")
  expect_columns(pls_fit(), page_fixture("pls"))
  m <- multipls_fit()
  page <- page_fixture("multipls")
  expect_columns(m, page, unit = "dim-1")
  expect_identical(names(outcomes(m, "dim-1")), page_columns(page, "outcomes"))
  expect_identical(names(docs(pls_fit(), unit = 2)), c("doc_id", "alignment_score"))
}))

test_that("unit: a pair in either order, unknown pairs, and NULL stacking every pair", with_empty_cache({
  g <- groups_fit()
  expect_identical(words(g, unit = c("g3", "g1")), words(g, unit = c("g1", "g3")))
  expect_identical(words(g, unit = c("mid", "high")), words(g, unit = c("g1", "g3")))
  expect_identical(docs(g, unit = "g2_g3"), docs(g, unit = c("g2", "g3")))
  expect_error(words(g, unit = c("g1", "g4")), "Unknown pair", fixed = TRUE)
  expect_error(words(g, unit = c("g1", "g1")), "Unknown pair", fixed = TRUE)
  # A label that is another group's internal name is ambiguous.
  g_clash <- g
  g_clash$group_labels[] <- c("g3", "x", "y")
  expect_error(words(g_clash, unit = c("g3", "x")), "'g3' is a group's label and the internal name", fixed = TRUE)
  all <- words(g)
  expect_identical(unique(all$contrast), c("g1_g2", "g1_g3", "g2_g3"))
  expect_identical(all[all$contrast == "g1_g3", ], `rownames<-`(words(g, unit = c("g1", "g3")), which(all$contrast == "g1_g3")))
  d <- docs(g)
  expect_identical(names(d)[1L], "contrast")
  expect_identical(as.vector(table(d$contrast)), as.vector(g$pairs$n_g1 + g$pairs$n_g2))
}))

test_that("unit: a term key and its position give the same rows; NULL stacks terms", with_empty_cache({
  r <- interaction_fit()
  keys <- c("ssd()", "grouplow:ssd()", "groupmid:ssd()")
  expect_identical(names(r$units), keys)
  expect_identical(words(r, unit = "grouplow:ssd()"), words(r, unit = 2))
  expect_identical(docs(r, unit = "groupmid:ssd()"), docs(r, unit = 3))
  all <- words(r, n = 3)
  expect_identical(names(all)[1L], "term")
  expect_identical(unique(all$term), keys)
  expect_identical(all[all$term == "ssd()", -1L], words(r, n = 3, unit = 1))
  expect_error(words(r, unit = "sex:ssd()"), "Unknown unit 'sex:ssd()'", fixed = TRUE)
  expect_error(words(r, unit = 4), "Unknown unit", fixed = TRUE)
  expect_error(words(pca_ols_fit(), unit = 1), "This result has one unit", fixed = TRUE)
}))

test_that("term units share y_hat and residual; each unit reads its own text", with_empty_cache({
  r <- interaction_fit()
  d <- docs(r)
  per_term <- split(d[c("doc_id", "y_true", "y_hat", "residual")], d$term)
  for (p in per_term) expect_identical(`rownames<-`(p, NULL), `rownames<-`(per_term[[1]], NULL))
  expect_equal(per_term[[1]]$residual, unname(residuals(r)))
  expect_false(identical(d$alignment_score[d$term == "ssd()"], d$alignment_score[d$term == "grouplow:ssd()"]))

  t <- results_text()
  rev_tokens <- lapply(rev(t$s$docs), identity)
  diary <- ssd(rev_tokens, t$emb, lexicon = t$d$lexicon)
  two <- fit_regression(list(essay = t$s, diary = diary), y ~ ssd(essay) + ssd(diary), data = t$data, fixed_k = 4)
  expect_identical(vapply(two$units, `[[`, 0L, "text"), c(`ssd(essay)` = 1L, `ssd(diary)` = 2L))
  expect_identical(doc_text(two, 5, unit = "ssd(diary)"), paste(rev_tokens[[5]], collapse = " "))
  expect_identical(doc_text(two, 5, unit = "ssd(essay)"), t$d$text[5])
  expect_error(doc_text(two, 5), "several texts", fixed = TRUE)
  expect_error(snippets(two, unit = "ssd(diary)"), "No snippets: this text came in as token lists", fixed = TRUE)
  expect_s3_class(snippets(two, top_per_side = 2, unit = "ssd(essay)"), "data.frame")
}))

test_that("a text from token lists: snippets raise, doc_text joins the tokens", with_empty_cache({
  t <- results_text()
  s <- ssd(t$s$docs, t$emb, lexicon = t$d$lexicon)
  g <- fit_groups(s, "group", data = t$data, n_perm = 10)
  expect_error(snippets(g, unit = c("g1", "g2")), "No snippets: this text came in as token lists", fixed = TRUE)
  expect_error(cluster_snippets(g, unit = c("g1", "g2")), "No snippets: this text came in as token lists", fixed = TRUE)
  expect_identical(doc_text(g, c(3, 7), unit = c("g1", "g2")),
                   c(paste(s$docs[[3]], collapse = " "), paste(s$docs[[7]], collapse = " ")))
  mid <- which(t$data$group == "mid")[1L]              # g3: not in the pair g1_g2
  expect_error(doc_text(g, mid, unit = c("g1", "g2")), sprintf("No doc with doc_id=%d", mid), fixed = TRUE)
}))

test_that("clusters are the same on every call, and cluster_snippets() uses their ids", with_empty_cache({
  g <- groups_fit()
  a <- clusters(g, topn = 60, unit = c("g1", "g2"))
  expect_identical(clusters(g, topn = 60, unit = c("g1", "g2")), a)
  expect_identical(clusters(g, topn = 60, random_state = 7, unit = c("g1", "g2")),
                   clusters(g, topn = 60, random_state = 7, unit = c("g1", "g2")))
  cs <- cluster_snippets(g, top_per_cluster = 3, topn = 60, unit = c("g1", "g2"))
  for (side in c("pos", "neg")) {
    expect_setequal(cs$cluster_id[cs$side == side], a$cluster_id[a$side == side])
    expect_identical(cs$snippet_id[cs$side == side], seq_len(sum(cs$side == side)))
  }
  cw <- cluster_words(g, topn = 60, unit = c("g1", "g2"))
  expect_identical(as.vector(table(factor(paste(cw$side, cw$cluster_id), unique(paste(a$side, a$cluster_id))))), a$size)
  expect_identical(a$top_words[1L], paste(utils::head(cw$word, 5L), collapse = ", "))
}))

test_that("words: n per side, the sides, cosines signed by side", with_empty_cache({
  r <- pca_ols_fit()
  w <- words(r, n = NULL)
  expect_identical(as.vector(table(w$side)[c("pos", "neg")]), c(23L, 23L))
  expect_true(all(w$cos_beta[w$side == "pos"] >= 0) && all(w$cos_beta[w$side == "neg"] <= 0))
  expect_identical(words(r, side = "neg", n = 5), `rownames<-`(w[w$side == "neg", ][1:5, ], NULL))
  expect_identical(w$rank[w$side == "pos"], 1:23)
  expect_error(words(r, side = "up"), "side must be 'pos', 'neg' or 'both'", fixed = TRUE)
}))

test_that("docs: pos, neg and misdiagnosed rows", with_empty_cache({
  r <- pca_ols_fit()
  all <- docs(r)
  expect_identical(all$doc_id, 1:240)
  over <- docs(r, "misdiagnosed", n = 10, direction = "over")
  expect_true(all(over$residual < 0))
  expect_identical(over$residual, sort(all$residual)[1:10])
  under <- docs(r, "misdiagnosed", n = 10, direction = "under")
  expect_true(all(under$residual > 0))
  both <- docs(r, "misdiagnosed", n = 10)
  expect_identical(abs(both$residual), sort(abs(all$residual), decreasing = TRUE)[1:10])
  expect_identical(docs(r, "pos", n = 3)$y_hat, sort(all$y_hat, decreasing = TRUE)[1:3])
  expect_identical(docs(r, "neg", n = 3)$y_hat, sort(all$y_hat)[1:3])
  g <- groups_fit()
  expect_error(docs(g, "misdiagnosed", unit = c("g1", "g2")), "needs residuals", fixed = TRUE)
  expect_identical(docs(g, "pos", n = 2, unit = c("g1", "g2"))$alignment_score,
                   sort(g$units$g1_g2$alignment_scores, decreasing = TRUE)[1:2])
  expect_identical(doc_text(r, 4), golden_texts()$text[4])
}))

test_that("sweep_table() and plot_sweep()", with_empty_cache({
  r <- pca_ols_fit()
  expect_identical(sweep_table(r), r$texts[[1]]$sweep)
  expect_identical(sweep_table(r, 1), r$texts[[1]]$sweep)
  expect_error(sweep_table(r, 2), "Unknown text 2; the fit has one text, at position 1.", fixed = TRUE)
  expect_identical(nrow(sweep_table(interaction_fit())), 0L)
  expect_identical(names(sweep_table(interaction_fit())), names(r$texts[[1]]$sweep))
  expect_error(plot_sweep(interaction_fit()), "fixed_k", fixed = TRUE)
  skip_if_not(capabilities("png"))
  path <- tempfile(fileext = ".png")
  expect_null(plot_sweep(r, path = path, dpi = 50))
  expect_identical(readBin(path, "raw", 4L), as.raw(c(0x89, 0x50, 0x4e, 0x47)))
  # One K: no detrended score or change to draw, the axes still are.
  t <- results_text()
  one <- fit_regression(t$s, "y", data = t$data, k_min = 4, k_max = 4)
  expect_identical(sweep_table(one)$k, 4L)
  unlink(path)
  expect_null(plot_sweep(one, path = path, dpi = 50))
  expect_identical(readBin(path, "raw", 4L), as.raw(c(0x89, 0x50, 0x4e, 0x47)))
}))

test_that("a saved result finds its pack in the cache in a fresh session", with_empty_cache({
  g <- groups_fit()
  expected <- words(g, unit = c("g1", "g2"))
  cache <- Sys.getenv("WORDPACK_CACHE")
  dir.create(cache, recursive = TRUE)
  file.copy(test_pack(), file.path(cache, "test-1.wordpack"))
  rds <- tempfile(fileext = ".rds")
  out <- tempfile(fileext = ".rds")
  saveRDS(g, rds)
  log <- in_fresh_session(sprintf("saveRDS(ssdiff::words(readRDS(%s), unit = c('g1', 'g2')), %s)",
                                  deparse(rds), deparse(out)))
  expect_true(file.exists(out), info = paste(log, collapse = "\n"))
  expect_identical(readRDS(out), expected)

  unlink(file.path(cache, "test-1.wordpack"))
  rm(list = ls(.registry, all.names = TRUE), envir = .registry)
  expect_error(words(readRDS(rds)), "This result was fitted on wordpack test-1. Load it with load_wordpack(\"test-1\")",
               fixed = TRUE)
}))

test_that("a saved mixed-model result works in a fresh session", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- regression_text()
  res <- suppressWarnings(fit_regression(r$s, y ~ cov + (1 | g), data = r$data, fixed_k = 3))
  rds <- tempfile(fileext = ".rds")
  saveRDS(res, rds)
  # Each generic first in its own session, before anything loads fastglmm.
  calls <- list(residuals = "residuals(res)", logLik = "logLik(res)",
                summary = "{options(useFancyQuotes = FALSE); capture.output(summary(res))}")
  for (name in names(calls)) {
    out <- tempfile(fileext = ".rds")
    log <- in_fresh_session(sprintf("res <- readRDS(%s); saveRDS(%s, %s)", deparse(rds), calls[[name]], deparse(out)))
    expect_true(file.exists(out), info = paste(c(name, log), collapse = "\n"))
  }
  expect_identical(readRDS(out), capture.output(summary(res)))
}))

test_that("a mixed-model result without fastglmm says how to install it", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- regression_text()
  res <- suppressWarnings(fit_regression(r$s, y ~ cov + (1 | g), data = r$data, fixed_k = 3))
  fake <- need_pkg
  environment(fake) <- list2env(list(requireNamespace = function(...) FALSE), parent = environment(need_pkg))
  local_mocked_bindings(need_pkg = fake)
  expect_error(coef(res), "this fit_regression() result needs the 'fastglmm' package", fixed = TRUE)
}))
