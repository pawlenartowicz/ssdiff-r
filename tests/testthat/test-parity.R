# R against Python's parity cases (spec 01, "Cross-language tests"), one
# test per case file: the committed golden cases on the test pack, or, when
# SSDIFF_PARITY_DIR is set (tests/parity/run.py), the reference cases on the
# real pack in that folder. Identical: tokens, n_kept, chosen K, word lists,
# cluster membership, lexicon tokens and counts, term keys and test labels.
# Relative 1e-6 on R's vectors: document vectors, betas, r2, effect sizes,
# cosines. Group permutation p-values identical, on R's vectors and on
# Python's. On Python's vectors put into s$x: plskit p-values and statistics
# 1e-9, fit_regression p-values 1e-8.

# The text-level checks: the first tokens, n_kept, the kept rows and every document vector.
expect_parity_texts <- function(case, texts) {
  tokens <- per_text(case, "tokens")
  n_kept <- per_text(case, "n_kept")
  for (name in names(texts)) {
    s <- texts[[name]]
    ref <- lapply(tokens[[name]], function(doc) as.character(unlist(doc)))
    expect_identical(s$docs[seq_along(ref)], ref, label = paste(name, "tokens"))
    expect_identical(s$n_kept, as.integer(n_kept[[name]]), label = paste(name, "n_kept"))
    v <- parity_vectors(case$texts[[name]]$vectors)
    expect_identical(which(s$keep), as.integer(v$doc_id) + 1L, label = paste(name, "kept rows"))
    expect_equal(s$x, v$x, tolerance = 1e-6, label = paste(name, "document vectors"))
  }
}

# Python's units keyed as R keys them: the text name of a two-text case is ssd(name).
py_units <- function(case) {
  u <- case$units
  if (length(case$texts) > 1L) names(u) <- paste0("ssd(", names(u), ")")
  u
}

# Betas of every unit (relative 1e-6) and their words and clusters (on Python's betas).
expect_parity_units <- function(case, res) {
  emb <- load_wordpack(case_pack())
  units <- py_units(case)
  keys <- setdiff(names(units), "main")
  if (length(keys)) expect_identical(names(res$units), keys)
  for (key in names(units)) {
    beta <- if (key == "main") res$beta else res$units[[key]]$beta
    expect_equal(beta, unlist(units[[key]]$beta), tolerance = 1e-6, label = paste(key, "beta"))
    expect_unit_words(units[[key]], emb, key, case$cluster_args)
  }
}

# pls_info of a gradients fit; order made 1-based.
expect_pls_info <- function(res, case) expect_py_table(res$pls_info, case$pls_info, 1e-6, ids = "order", label = "pls_info")

# The sweep table of text i against Python's: K and n_clusters identical, the
# rest relative 1e-6. The fixtures stop the sweep below the rank of the
# document vectors (k_max = min(120, D - 2): the ABTT takes one direction out of
# every word vector), where every K is a real fit in both languages.
expect_parity_sweep <- function(res, i, ref) {
  expect_py_table(sweep_table(res, i), ref, 1e-6, label = paste("text", i, "sweep"))
}

check_regression <- function(case) {
  run <- run_parity_case(case)
  res <- run$res
  expect_parity_texts(case, run$texts)
  fit_info <- per_text(case, "fit_info")
  sweep <- per_text(case, "sweep")
  for (i in seq_along(run$texts)) {
    name <- names(run$texts)[i]
    # Python keeps the formula string it was given; R deparses the formula it
    # was given, which names the text (ssd(text)) where Python names it after its column.
    expect_py_table(.text_fit_info(res, i), fit_info[[name]], 1e-6, skip = "formula", label = paste(name, "fit_info"))
    expect_identical(res$texts[[i]]$k, as.integer(py_table(fit_info[[name]])$n_components))
    expect_parity_sweep(res, i, sweep[[name]])
    expect_identical(res$texts[[i]]$n_no_vector, as.integer(per_text(case, "n_no_vector")[[name]]))
  }
  expect_identical(which(res$texts[[1L]]$keep), as.integer(unlist(case$kept_doc_id)) + 1L)
  # r2_adj and the p-value are compared on Python's vectors: near 0, r2_adj's
  # relative change is larger than r2's (tokens_raw: 9e-6).
  expect_py_table(res$stats, case$stats, 1e-6, skip = c("pvalue", "r2_adj"), label = "stats")
  test <- py_table(case$test)
  expect_identical(res$test$label, if (is.null(test$label)) NA_character_ else test$label)
  expect_identical(res$test$name, test$name)
  expect_identical(res$terms$term, py_table(case$terms)$label)
  expect_parity_units(case, res)
  if (!is.null(case$pls_info)) expect_pls_info(res, case)

  res <- run_parity_case(case, python = TRUE)$res
  expect_identical(vapply(res$texts, `[[`, 0L, "k"), vapply(run$res$texts, `[[`, 0L, "k"))
  st <- py_table(case$stats)
  expect_equal(c(res$stats$r2, res$stats$r2_adj, res$stats$pvalue), c(st$r2, st$r2_adj, st$pvalue), tolerance = 1e-8,
               label = "r2, r2_adj, pvalue")
  expect_py_table(res$test, case$test, 1e-8, label = "test")
  terms <- res$terms
  names(terms)[names(terms) == "term"] <- "label"
  expect_py_table(terms, case$terms, 1e-8, label = "terms")
  mt <- case$model_test
  expect_equal(unlist(res$model_test), py_num(mt[c("fvalue", "f_pvalue", "llr", "llr_pvalue")]), tolerance = 1e-8,
               label = "model_test")
  expect_py_table(res$covariates, case$covariates, 1e-8, label = "covariates")
}

check_pls <- function(case) {
  run <- run_parity_case(case)
  res <- run$res
  expect_parity_texts(case, run$texts)
  expect_py_table(res$fit_info, case$fit_info, 1e-6, label = "fit_info")
  expect_py_table(res$stats, case$stats, 1e-6, skip = "pvalue", label = "stats")
  expect_parity_units(case, res)
  if (!is.null(case$pls_info)) expect_pls_info(res, case)

  res <- run_parity_case(case, python = TRUE)$res
  expect_identical(res$fit_info$n_components, as.integer(py_table(case$fit_info)$n_components))
  expect_py_table(res$test, case$test, 1e-9, label = "test")
  expect_equal(res$stats$pvalue, py_table(case$stats)$pvalue, tolerance = 1e-9)
}

check_multipls <- function(case) {
  run <- run_parity_case(case)
  res <- run$res
  expect_parity_texts(case, run$texts)
  expect_py_table(res$stats, case$stats, 1e-6, skip = "pvalue", label = "stats")
  expect_pls_info(res, case)
  expect_parity_units(case, res)

  res <- run_parity_case(case, python = TRUE)$res
  expect_py_table(res$test, case$test, 1e-9, label = "test")
  expect_equal(res$stats$pvalue, py_table(case$stats)$pvalue, tolerance = 1e-9)
}

check_groups <- function(case) {
  run <- run_parity_case(case)
  res <- run$res
  expect_parity_texts(case, run$texts)
  expect_identical(names(res$units), unlist(case$keys))
  expect_py_table(res$pairs, case$pairs, 1e-6, label = "pairs")
  # Permutation p-values are counts: identical on R's vectors too.
  expect_py_table(res$stats, case$stats, 0, label = "stats")
  expect_identical(res$test$omnibus_p, py_table(case$test)$omnibus_p)
  expect_identical(res$pairs[c("p_raw", "p_corrected")], py_table(case$pairs)[c("p_raw", "p_corrected")])
  expect_parity_units(case, res)

  # On Python's vectors the statistics are the same doubles up to summation
  # order, and the permutation p-values (counts) are identical.
  res <- run_parity_case(case, python = TRUE)$res
  expect_py_table(res$stats, case$stats, 0, label = "stats")
  expect_py_table(res$test, case$test, 1e-9, label = "test")
  expect_identical(res$test$omnibus_p, py_table(case$test)$omnibus_p)
  expect_py_table(res$pairs, case$pairs, 1e-9, label = "pairs")
  expect_identical(res$pairs[c("p_raw", "p_corrected")], py_table(case$pairs)[c("p_raw", "p_corrected")])
  for (key in names(res$units)) {
    expect_equal(res$units[[key]]$beta, unlist(case$units[[key]]$beta), tolerance = 1e-9, label = key)
  }
}

check_contrast <- function(case) {
  emb <- load_wordpack(case_pack())
  cc <- contrast(unlist(case$g1$gradient), unlist(case$g2$gradient), emb, labels = unlist(case$kwargs$labels))
  expect_py_table(cc$stats, case$stats, 1e-9, label = "stats")
  expect_identical(cc$labels, unlist(case$kwargs$labels))
  expect_parity_units(case, cc)
}

check_lexicon <- function(case) {
  s <- ssd(case_texts()$text, load_wordpack(case_pack()), use_full_doc = FALSE)
  args <- unlist(case$args[-1L])
  got <- do.call(case$call, c(list(s, case_texts()[[case$args[[1L]]]]), if (length(args)) list(args), case$kwargs))
  expect_lexicon_equal(got, case)
}

checks <- list(fit_regression = check_regression, ssdiff.fit_regression = check_regression, fit_pls = check_pls,
               fit_multipls = check_multipls, fit_groups = check_groups, ssdiff.contrast = check_contrast,
               suggest_lexicon = check_lexicon, evaluate_lexicon = check_lexicon)

test_that("every parity case has a check", {
  calls <- vapply(parity_case_names(), function(name) parity_case(name)$call, "")
  expect_true(all(calls %in% names(checks)), info = paste(setdiff(calls, names(checks)), collapse = ", "))
})

# The Suggests packages a case's R call needs (spec 01, "Which functions need
# which Suggests package"): plskit for the PLS fits and for ssd_gradients > 1,
# fastglmm for random effects.
parity_needs <- function(case) {
  plskit <- case$call %in% c("fit_pls", "fit_multipls") || isTRUE(case$kwargs$ssd_gradients > 1)
  fastglmm <- any(grepl("|", unlist(case$args), fixed = TRUE))
  c("plskit", "fastglmm")[c(plskit, fastglmm)]
}

for (name in parity_case_names()) {
  test_that(paste("R matches Python:", name), with_empty_cache({
    case <- parity_case(name)
    for (pkg in parity_needs(case)) skip_if_not_installed(pkg)
    checks[[case$call]](case)
  }))
}
