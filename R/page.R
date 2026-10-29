# Port of ssdiff/results/page.py. Change the two together.

# The version of the page's JSON; ssdiff/results/page.py has the same number.
REPORT_FORMAT <- 2L

# The cluster arguments every table of the page uses: the accessors' defaults.
.PAGE_CLUSTER_ARGS <- list(topn = 100, k = NULL, k_min = 2, k_max = 10, random_state = 2137, min_cluster_size = 2)

.PAGE_KINDS <- c(ssd_pls = "pls", ssd_pca_ols = "pca_ols", ssd_groups = "groups", ssd_multipls = "multipls",
                 ssd_contrast = "contrast")

#' The results page
#'
#' `save_page()` writes a result as one self-contained HTML file: the
#' results page of Python ssdiff, which opens from disk in a browser, with no
#' network. It shows the words, clusters, snippets and documents of every
#' unit of the result, its tests and its settings. `page_data()` gives the
#' data the page reads, as a nested list.
#'
#' Every table is computed when the page is written, with the accessors'
#' default cluster arguments ([clusters()]). A text that came in as token
#' lists has no snippet tables, and its documents' text is the tokens joined
#' by spaces.
#'
#' Of the documents, the page keeps the text of the `doc_texts` highest and
#' lowest by alignment score, and for [fit_pls()] and [fit_regression()] also
#' the `doc_texts` most over- and under-predicted; the other documents are
#' listed without text.
#'
#' In `page_data()`, a JSON object is a named list and a JSON array an
#' unnamed list or an [I()] vector. Every table is
#' `list(columns = I(names(df)), rows = df)`, with `df` a data frame. Ids
#' are as in the page: `cluster_id` and a unit's `text` count from 0, as
#' the page counts them; `doc_id` is the document's row in the text, from 1.
#'
#' @param result A result of [fit_pls()], [fit_multipls()],
#'   [fit_regression()], [fit_groups()] or [contrast()].
#' @param path The file to write.
#' @param title The page's title; `NULL` for `"SSD results"`.
#' @param dataset A named list or named character vector describing the data,
#'   shown on the page; `NULL` for none.
#' @param snippets_per_side The number of snippets per side of each unit.
#' @param snippets_per_cluster The number of snippets per cluster.
#' @param doc_texts The number of documents at each end of each ranking whose
#'   text goes into the page; `Inf` for every document, `0` for none.
#' @return `page_data()`: a named list. `save_page()`: `path`, invisibly.
#' @examples
#' \dontrun{
#' r <- fit_regression(s, "score", data = df)
#' save_page(r, "results.html", title = "Essays", dataset = c(File = "essays.csv"))
#' str(page_data(r), max.level = 1)
#' }
#' @export
page_data <- function(result, title = NULL, dataset = NULL, snippets_per_side = 200, snippets_per_cluster = 50,
                      doc_texts = 50) {
  hit <- intersect(class(result), names(.PAGE_KINDS))
  if (!inherits(result, "ssd_result") || !length(hit)) {
    stop(sprintf("no results page for %s", class(result)[1L]), call. = FALSE)
  }
  kind <- .PAGE_KINDS[[hit[1L]]]
  if (!is.null(title) && (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("title must be a string or NULL, got ", deparse1(title), call. = FALSE)
  }
  if (!is.null(dataset) && (!(is.list(dataset) || is.atomic(dataset)) ||
                            (length(dataset) && (is.null(names(dataset)) || any(!nzchar(names(dataset))))))) {
    stop("dataset must be a named list or a named character vector, or NULL", call. = FALSE)
  }
  .check_n(snippets_per_side, "snippets_per_side")
  .check_n(snippets_per_cluster, "snippets_per_cluster", null_ok = FALSE)
  .check_n(doc_texts, "doc_texts", null_ok = FALSE)
  a <- list(sps = snippets_per_side, spc = snippets_per_cluster, doc_texts = doc_texts)

  built <- if (kind == "contrast") {
    .page_contrast(result, a)
  } else if (kind == "groups") {
    .page_groups(result, a)
  } else if (kind == "multipls") {
    .page_multipls(result, a)
  } else if (is.null(result$pls_info) && !is.null(result$units)) {
    .page_terms(result, a)
  } else {
    .page_continuous(result, kind, a)
  }
  data <- list(
    report_format = REPORT_FORMAT,
    kind = kind,
    title = if (is.null(title) || !nzchar(title)) "SSD results" else title,
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    software = list(package = "ssdiff", language = "r", version = as.character(utils::packageVersion("ssdiff"))),
    texts = lapply(result$texts, .page_text_entry, kind = kind),
    settings = list(fit = .page_fit_args(result$fit_args), clusters = .PAGE_CLUSTER_ARGS),
    tables = built$tables,
    units = built$units
  )
  if (!is.null(dataset)) data$dataset <- if (length(dataset)) as.list(dataset) else stats::setNames(list(), character(0))
  data
}

#' @rdname page_data
#' @export
save_page <- function(result, path = "ssd_results.html", title = NULL, dataset = NULL, snippets_per_side = 200,
                      snippets_per_cluster = 50, doc_texts = 50) {
  json <- page_json(page_data(result, title = title, dataset = dataset, snippets_per_side = snippets_per_side,
                              snippets_per_cluster = snippets_per_cluster, doc_texts = doc_texts))
  tpl_path <- system.file("page_template.html", package = "ssdiff")
  tpl <- rawToChar(readBin(tpl_path, "raw", file.size(tpl_path)))
  Encoding(tpl) <- "UTF-8"
  # fixed = TRUE: a backslash in the JSON is not read as a backreference.
  html <- sub("__SSD_DATA__", json, tpl, fixed = TRUE)
  writeBin(charToRaw(enc2utf8(html)), path)
  invisible(path)
}

# The page's JSON text. Every "<" becomes <, so no "</script>" or
# "<!--<script" in a string can end the data element early.
page_json <- function(data) {
  json <- jsonlite::toJSON(data, auto_unbox = TRUE, null = "null", na = "null", digits = NA, dataframe = "values")
  gsub("<", "\\u003c", enc2utf8(as.character(json)), fixed = TRUE)
}

# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

.page_table <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  class(df) <- "data.frame"
  rownames(df) <- NULL
  list(columns = I(names(df)), rows = df)
}

# pls_info's order and signs cells stay arrays when they hold one value.
.page_pls_info <- function(df) {
  for (col in intersect(c("order", "signs"), names(df))) df[[col]] <- I(lapply(df[[col]], function(v) I(unclass(v))))
  .page_table(df)
}

# Python's test view of a fit with one text term has no `label` column.
.page_test <- function(df, labelled) {
  .page_table(if (labelled) df else df[setdiff(names(df), "label")])
}

# Words, clusters and, when the text has sentences, snippets of one unit.
# `cluster_id` counts from 0; the tables other than words carry no contrast,
# as Python writes them.
.page_unit_tables <- function(result, unit, text, a) {
  no_contrast <- function(df) {
    df$contrast <- rep(NA_character_, nrow(df))
    df
  }
  from0 <- function(df) {
    df$cluster_id <- df$cluster_id - 1L
    df
  }
  cl <- from0(no_contrast(clusters(result, unit = unit)))
  cl$top_snippet <- rep("", nrow(cl))
  t <- list(words = .page_table(words(result, n = NULL, unit = unit)), clusters = .page_table(cl),
            cluster_words = .page_table(from0(no_contrast(cluster_words(result, unit = unit)))))
  if (!is.null(result$texts[[text]]$pre_docs)) {
    t$snippets <- .page_table(no_contrast(snippets(result, top_per_side = a$sps, unit = unit)))
    t$cluster_snippets <- .page_table(from0(no_contrast(cluster_snippets(result, top_per_cluster = a$spc, unit = unit))))
  }
  t
}

# A unit's docs table with `text`: the text of the documents at both ends of
# each ranking (alignment score, and residual where there is one).
.page_docs <- function(result, unit, text, doc_texts) {
  d <- docs(result, unit = unit)
  rankings <- list(d$alignment_score)
  if ("residual" %in% names(d)) rankings <- c(rankings, list(d$residual))
  keep <- if (is.infinite(doc_texts)) {
    d$doc_id
  } else if (doc_texts <= 0) {
    integer(0)
  } else {
    unlist(lapply(rankings, function(v) {
      o <- order(v)
      d$doc_id[c(utils::head(o, doc_texts), utils::tail(o, doc_texts))]
    }))
  }
  t <- result$texts[[text]]
  d$text <- rep(NA_character_, nrow(d))
  has <- d$doc_id %in% keep
  ids <- d$doc_id[has]
  d$text[has] <- if (is.null(t$pre_docs)) {
    vapply(t$docs[ids], paste, "", collapse = " ")
  } else {
    vapply(t$pre_docs[ids], `[[`, "", "raw")
  }
  .page_table(d)
}

# ---------------------------------------------------------------------------
# Result kinds
# ---------------------------------------------------------------------------

.page_unit <- function(key, label, pos_label, neg_label, tables, text) {
  list(key = key, label = label, pos_label = pos_label, neg_label = neg_label, tables = tables, text = text - 1L)
}

# ssd_pls, and ssd_pca_ols with one text term: the fit's own unit, and with
# ssd_gradients > 1 one more unit per rotated gradient.
.page_continuous <- function(result, kind, a) {
  tables <- list(stats = .page_table(result$stats), fit_info = .page_table(result$fit_info),
                 test = .page_test(result$test, labelled = FALSE))
  if (!is.null(result$covariates) && nrow(result$covariates)) tables$covariates <- .page_table(result$covariates)
  dims <- !is.null(result$pls_info)
  main <- c(.page_unit_tables(result, NULL, 1L, a), list(docs = .page_docs(result, NULL, 1L, a$doc_texts)))
  units <- list(.page_unit("main", if (dims) "\u03b2" else "", "+\u03b2", "\u2212\u03b2", main, 1L))
  if (dims) {
    tables$pls_info <- .page_pls_info(result$pls_info)
    keys <- names(result$units)
    for (i in seq_along(keys)) {
      dim <- c(list(stats = .page_table(result$units[[i]]$stats)), .page_unit_tables(result, keys[[i]], 1L, a),
               list(docs = .page_docs(result, keys[[i]], 1L, a$doc_texts)))
      units[[i + 1L]] <- .page_unit(keys[[i]], sprintf("Dim %d", i), "+", "\u2212", dim, 1L)
    }
  }
  list(tables = tables, units = units)
}

# ssd_pca_ols with several text terms: one unit per term, each with its own
# stats and fit_info; the top-level fit_info is the first unit's.
.page_terms <- function(result, a) {
  tables <- list(stats = .page_table(result$stats), fit_info = .page_table(result$units[[1L]]$fit_info),
                 test = .page_test(result$test, labelled = TRUE))
  if (!is.null(result$covariates) && nrow(result$covariates)) tables$covariates <- .page_table(result$covariates)
  units <- lapply(names(result$units), function(key) {
    u <- result$units[[key]]
    tabs <- c(list(stats = .page_table(u$stats), fit_info = .page_table(u$fit_info)),
              .page_unit_tables(result, key, u$text, a), list(docs = .page_docs(result, key, u$text, a$doc_texts)))
    list(key = key, label = key, pos_label = "+\u03b2", neg_label = "\u2212\u03b2", text = u$text - 1L,
         test = if (is.na(u$test)) NULL else u$test, tables = tabs)
  })
  list(tables = tables, units = units)
}

.page_groups <- function(result, a) {
  tables <- list(stats = .page_table(result$stats), test = .page_table(result$test),
                 pairs = .page_table(result$pairs))
  labels <- result$group_labels
  lab <- function(g) as.character(if (g %in% names(labels)) labels[[g]] else g)
  units <- lapply(names(result$units), function(key) {
    u <- result$units[[key]]
    tabs <- c(.page_unit_tables(result, key, 1L, a), list(docs = .page_docs(result, key, 1L, a$doc_texts)))
    .page_unit(key, paste(lab(u$g1), "vs", lab(u$g2)), lab(u$g1), lab(u$g2), tabs, 1L)
  })
  list(tables = tables, units = units)
}

.page_multipls <- function(result, a) {
  tables <- list(stats = .page_table(result$stats), pls_info = .page_pls_info(result$pls_info),
                 test = .page_table(result$test))
  keys <- names(result$units)
  units <- lapply(seq_along(keys), function(i) {
    key <- keys[[i]]
    tabs <- c(list(stats = .page_table(result$units[[i]]$stats), outcomes = .page_table(outcomes(result, key))),
              .page_unit_tables(result, key, 1L, a), list(docs = .page_docs(result, key, 1L, a$doc_texts)))
    .page_unit(key, sprintf("Dim %d", i), "+", "\u2212", tabs, 1L)
  })
  list(tables = tables, units = units)
}

.page_contrast <- function(result, a) {
  lab <- result$labels
  units <- list(.page_unit("contrast", paste(lab[1L], "\u2212", lab[2L]), lab[1L], lab[2L],
                           .page_unit_tables(result, NULL, 1L, a), 1L))
  list(tables = list(stats = .page_table(result$stats)), units = units)
}

# ---------------------------------------------------------------------------
# Texts and settings
# ---------------------------------------------------------------------------

.page_text_entry <- function(t, kind) {
  ref <- t$pack
  emb <- .embeddings_from_ref(ref)
  embeddings <- list(name = ref$name, version = ref$version, language = emb$language, n = emb$pack$n,
                     dims = emb$pack$dims, abtt = ref$abtt)
  if (kind == "contrast") {
    return(list(name = NULL, embeddings = embeddings, settings = NULL, k = NULL, sweep = NULL, n_no_vector = NULL,
                snippets = NULL))
  }
  pca <- kind == "pca_ols"
  list(
    name = t$name,
    embeddings = embeddings,
    settings = list(lexicon = if (length(t$lexicon)) I(sort(t$lexicon, method = "radix")),
                    window = t$window, sif_a = t$sif_a, use_full_doc = t$use_full_doc),
    k = if (pca) t$k,
    sweep = if (pca && !is.null(t$sweep) && nrow(t$sweep)) .page_table(t$sweep),
    n_no_vector = t$n_no_vector,
    snippets = !is.null(t$pre_docs)
  )
}

# fit_args as the page shows them; each group of outcome_groups stays an array.
.page_fit_args <- function(fit_args) {
  if (length(fit_args$outcome_groups)) fit_args$outcome_groups <- lapply(fit_args$outcome_groups, I)
  fit_args
}
