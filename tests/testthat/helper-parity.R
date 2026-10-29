# Python's parity cases (spec 01, "Cross-language tests"): reading a case,
# rebuilding its texts, running its R call, and comparing tables and units.

# Where the cases are: SSDIFF_PARITY_DIR, the reference cases on the real pack
# that tests/parity/run.py writes, else the committed golden cases on the test pack.
case_dir <- function() {
  dir <- Sys.getenv("SSDIFF_PARITY_DIR")
  if (nzchar(dir)) dir else test_path("fixtures", "golden")
}

# The pack the cases were fitted on: the file Python read (pack.json of the
# reference cases), else the test pack.
case_pack <- function() {
  if (!nzchar(Sys.getenv("SSDIFF_PARITY_DIR"))) return(test_pack())
  jsonlite::read_json(file.path(case_dir(), "pack.json"))$path
}

case_texts <- function() jsonlite::read_json(file.path(case_dir(), "texts.json"), simplifyVector = TRUE)

# Every case file; texts.json, pack.json, versions.json and the vectors_*.json files are not cases.
parity_case_names <- function() {
  f <- setdiff(list.files(case_dir(), pattern = "\\.json$"), c("texts.json", "pack.json", "versions.json"))
  sort(sub("\\.json$", "", f[!startsWith(f, "vectors_")]), method = "radix")
}

.parity_vectors <- new.env()

# Python's document vectors of a text configuration (read once per run).
parity_vectors <- function(file) {
  if (is.null(.parity_vectors[[file]])) {
    .parity_vectors[[file]] <- jsonlite::read_json(file.path(case_dir(), file), simplifyVector = TRUE)
  }
  .parity_vectors[[file]]
}

# A text's vectors replaced by Python's (0-based doc_id); the rows must be the text's own.
with_python_vectors <- function(s, file) {
  v <- parity_vectors(file)
  testthat::expect_identical(as.integer(v$doc_id) + 1L, which(s$keep))
  s$x <- v$x
  s
}

# A text-level field of the case as a list keyed by text name. A case with one
# text holds the field itself, a case with several holds it per name.
per_text <- function(case, field) {
  names <- names(case$texts)
  if (length(names) == 1L) return(stats::setNames(list(case[[field]]), names))
  case[[field]][names]
}

# Text `name` of the case, built as Python built it: the strings of its column
# in texts.json or the case's token lists, with its lexicon, window, sif_a and
# use_full_doc.
parity_text <- function(case, name, emb) {
  cfg <- case$texts[[name]]
  docs <- if (isTRUE(cfg$token_lists)) {
    lapply(per_text(case, "token_lists")[[name]], function(doc) as.character(unlist(doc)))
  } else {
    case_texts()[[cfg$column]]
  }
  ssd(docs, emb, lexicon = if (is.null(cfg$lexicon)) NULL else unlist(cfg$lexicon), window = cfg$window,
      sif_a = cfg$sif_a, use_full_doc = cfg$use_full_doc)
}

parity_data <- function() {
  d <- case_texts()
  data.frame(y = d$y, y2 = d$y2, cov = d$cov, group = d$group, stringsAsFactors = FALSE)
}

# The case's R call on its rebuilt texts; `python = TRUE` puts Python's
# document vectors into each text first. Returns the result and the texts.
run_parity_case <- function(case, python = FALSE) {
  emb <- load_wordpack(case_pack())
  texts <- lapply(stats::setNames(nm = names(case$texts)), function(name) {
    s <- parity_text(case, name, emb)
    if (python) with_python_vectors(s, case$texts[[name]]$vectors) else s
  })
  data <- parity_data()
  args <- unlist(case$args)
  call <- sub("^ssdiff\\.", "", case$call)
  res <- if (length(texts) > 1L) {
    do.call(fit_regression, c(list(texts, stats::as.formula(args[[1L]]), data = data), case$kwargs))
  } else if (call == "fit_regression") {
    # A text from strings is named after its column, as Python names it; token lists have no name.
    term <- if (isTRUE(case$texts[[1L]]$token_lists)) "ssd()" else "ssd(text)"
    rhs <- if (length(args)) sub("^~\\s*", "", args[[1L]]) else term
    if (!grepl("ssd(", rhs, fixed = TRUE)) rhs <- paste(rhs, "+", term)
    do.call(fit_regression, c(list(texts[[1L]], stats::as.formula(paste("y ~", rhs)), data = data), case$kwargs))
  } else {
    y <- unlist(case$ssd$y)
    if (length(args)) y <- stats::as.formula(paste(y, args[[1L]]))
    do.call(call, c(list(texts[[1L]], y, data = data), case$kwargs))
  }
  list(res = res, texts = texts)
}

# Python's {"columns", "rows"} tables as data frames (null becomes NA).
py_num <- function(v) vapply(v, function(x) if (is.null(x) || identical(x, NA)) NA_real_ else as.numeric(x), 0)

# R's table against Python's, column by column for Python's columns: JSON
# integers, strings and booleans identical, JSON floats to `tol` (relative).
# `ids` names columns of 0-based positions; `skip` names columns left out.
expect_py_table <- function(got, ref, tol, ids = character(0), skip = character(0), label = "") {
  py <- py_table(ref)
  testthat::expect_identical(nrow(got), nrow(py), label = paste(label, "rows"))
  for (col in setdiff(names(py), skip)) {
    lab <- paste(label, col)
    v <- py[[col]]
    if (is.list(v)) {
      g <- got[[col]]
      for (i in seq_along(v)) {
        ref_i <- unlist(v[[i]])
        if (col %in% ids) ref_i <- as.integer(ref_i) + 1L
        testthat::expect_equal(as.numeric(g[[i]]), as.numeric(ref_i), tolerance = tol, label = lab)
      }
    } else if (col %in% ids) {
      testthat::expect_identical(as.integer(got[[col]]), as.integer(v) + 1L, label = lab)
    } else if (is.double(v)) {
      testthat::expect_equal(as.numeric(got[[col]]), v, tolerance = tol, label = lab)
    } else if (is.integer(v)) {
      testthat::expect_identical(as.integer(got[[col]]), v, label = lab)
    } else if (is.character(v)) {
      testthat::expect_identical(as.character(got[[col]]), v, label = lab)
    } else if (is.logical(v) && !all(is.na(v))) {
      testthat::expect_identical(as.logical(got[[col]]), v, label = lab)
    } else {
      testthat::expect_true(all(is.na(got[[col]])), label = lab)  # null in every row
    }
  }
}

# Words, cluster words and clusters of a unit against Python's, on Python's
# beta: word lists, ranks, cluster ids, sides and membership identical, cosines
# relative 1e-6. `cluster_args` are the clustering arguments Python used.
expect_unit_words <- function(u, emb, key, cluster_args = list()) {
  res <- gradient_result(unlist(u$beta), emb)
  for (side in c("pos", "neg")) {
    py <- py_table(u$words[[side]])
    w <- words(res, side = side, n = nrow(py))
    testthat::expect_identical(w$word, py$word, label = paste(key, side, "words"))
    testthat::expect_identical(w$rank, as.integer(py$rank), label = paste(key, side, "ranks"))
    testthat::expect_equal(w$cos_beta, py$cos_beta, tolerance = 1e-6, label = paste(key, side, "cos_beta"))
  }
  cw <- do.call(cluster_words, c(list(res), cluster_args))
  expect_py_table(cw, u$cluster_words, 1e-6, ids = "cluster_id", skip = "contrast", label = paste(key, "cluster_words"))
  if (!is.null(u$clusters)) {
    cl <- do.call(clusters, c(list(res), cluster_args))
    expect_py_table(cl, u$clusters, 1e-6, ids = "cluster_id", skip = "contrast", label = paste(key, "clusters"))
  }
}

# Compares an ssd_lexicon with a Python LexiconResult written as tables:
# columns and order, tokens, seeds and counts exactly, the rest to `tol`.
expect_lexicon_equal <- function(got, ref, tol = 1e-9) {
  stats <- py_table(ref$stats)
  expect_identical(names(got$stats), names(stats))
  expect_identical(got$stats$var_type, stats$var_type)
  expect_equal(c(got$stats$n_docs, got$stats$n_tokens), c(stats$n_docs, stats$n_tokens))

  sug <- py_table(ref$suggestions)
  expect_identical(names(got$suggestions), names(sug))
  for (col in intersect(c("seed", "token", "direction"), names(sug))) {
    expect_identical(got$suggestions[[col]], as.character(sug[[col]]), label = col)
  }
  expect_identical(got$suggestions$freq, as.integer(sug$freq))
  for (col in intersect(c("cov_all", "cov_bal", "corr", "pvalue", "rank"), names(sug))) {
    expect_equal(got$suggestions[[col]], as.numeric(sug[[col]]), tolerance = tol, label = col)
  }
  # Python's cosines are float32.
  if ("cos" %in% names(sug)) expect_equal(got$suggestions$cos, as.numeric(sug$cos), tolerance = 1e-6)

  if (is.null(ref$summary)) return(expect_null(got$summary))
  sm <- py_table(ref$summary)
  group_cov <- sm$group_cov[[1]]
  sm$group_cov <- NULL
  if (!is.null(group_cov) && !identical(group_cov, NA)) {
    expect_equal(got$summary$group_cov[[1]], unlist(group_cov), tolerance = tol)
    got$summary$group_cov <- NULL
  }
  expect_identical(names(got$summary), names(sm))
  expect_equal(unlist(got$summary), unlist(lapply(sm, as.numeric)), tolerance = tol)
}
