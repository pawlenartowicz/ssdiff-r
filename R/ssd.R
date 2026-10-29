# Port of the SSD class of ssdiff/ssd.py (construction, with_lexicon, repr) and
# its outcome helpers. Change them together.

#' Build a text for SSD analysis
#'
#' `ssd()` tokenizes a text once and builds its document vectors. The result
#' knows nothing about an outcome: every fit ([fit_pls()], [fit_multipls()],
#' [fit_regression()], [fit_groups()]) takes it first and the outcome second.
#'
#' A document's vector is built from its words' embeddings, each weighted by
#' SIF, `sif_a / (sif_a + count / total)`, with word counts over all documents
#' of the text. With a lexicon (`use_full_doc = FALSE`), each occurrence of a
#' lexicon word gives the weighted mean of the words within `window` on each
#' side, the occurrence itself left out, and the document's vector is the mean
#' of those. With `use_full_doc = TRUE` it is the weighted mean of all the
#' document's words. A document with no vector is not kept; the others are
#' scaled to unit length.
#'
#' Without a lexicon and with `use_full_doc = FALSE`, the object holds tokens
#' only. [suggest_lexicon()] and [evaluate_lexicon()] work on it, which makes it
#' the cheap way to tokenize before choosing a lexicon; [with_lexicon()] then
#' adds the vectors. The fits need vectors and raise on it.
#'
#' @section Two uses of the name `ssd`:
#' `ssd()` is also the text term in a [fit_regression()] formula, as in
#' `score ~ age + sex * ssd()` or `score ~ ssd(essay)`. There it is never
#' called: `fit_regression()` reads the term from the formula and replaces it
#' with the text's principal components before any model function sees the
#' formula. Passing such a formula to [lm()] or [glm()] directly calls this
#' constructor, which fails on its missing `texts`.
#'
#' @param texts A character vector (or a factor), one text per document (`NA`
#'   is an empty document), tokenized with the pack's tokenizer and
#'   stopwords. Or a list of character vectors, one token list per document,
#'   used as given.
#' @param embeddings Embeddings from [load_wordpack()]. `NULL` loads
#'   `load_wordpack()`, the `"default"` tag.
#' @param lexicon Character vector of lexicon words, or `NULL`.
#' @param window Context window: the number of words on each side of a lexicon
#'   word, a whole number of at least 1.
#' @param sif_a SIF smoothing parameter, greater than 0. `NULL` weighs every
#'   word 1 (a plain mean).
#' @param use_full_doc `TRUE` builds each document's vector from all its words,
#'   `FALSE` from the windows around the lexicon's words. `NULL` means `TRUE`
#'   without a lexicon and `FALSE` with one.
#' @return An object of class `ssd`, a list with `docs` (the tokens of each
#'   document), `pre_docs` (sentence structure for snippets, `NULL` for token
#'   lists), `tokenizer`, `x` (one unit row per document with a vector, `NULL`
#'   for tokens only), `keep` (whether each document has a vector), `n_docs`,
#'   `n_kept`, `lexicon`, `window`, `sif_a`, `use_full_doc`, and internally the
#'   SIF counts and a reference to the pack (`pack`).
#' @seealso [with_lexicon()], [suggest_lexicon()]
#' @examples
#' \dontrun{
#' emb <- load_wordpack()
#' s <- ssd(df$essay, embeddings = emb)                    # full documents
#' s <- ssd(df$essay, embeddings = emb, lexicon = c("work", "job"))
#' s0 <- ssd(df$essay, embeddings = emb, use_full_doc = FALSE)   # tokens only
#' s1 <- with_lexicon(s0, c("work", "job"))
#' }
#' @export
ssd <- function(texts, embeddings = NULL, lexicon = NULL, window = 3, sif_a = 1e-3,
                use_full_doc = NULL) {
  .check_window_sif(window, sif_a)
  if (!is.null(use_full_doc)) .check_flag(use_full_doc, "use_full_doc")
  .check_lexicon(lexicon)
  if (!length(texts)) stop("texts is empty; pass one text per document.", call. = FALSE)
  if (is.factor(texts)) texts <- as.character(texts)
  .texts_form(texts)
  if (is.null(embeddings)) embeddings <- load_wordpack()
  .check_embeddings(embeddings, "embeddings")
  tok <- .tokenize_texts(texts, embeddings)
  if (is.null(use_full_doc)) use_full_doc <- is.null(lexicon)
  .new_ssd(tok, embeddings, lexicon, window, sif_a, use_full_doc)
}

#' Give a text another lexicon
#'
#' A new [ssd()] object on the same tokens, with new document vectors. The
#' defaults are `ssd()`'s, not `s`'s settings. `ssd(texts, emb, lexicon)`
#' equals `with_lexicon(ssd(texts, emb), lexicon)`, so a text can be tokenized
#' once and fitted with several lexicons.
#'
#' @param s A text from [ssd()].
#' @inheritParams ssd
#' @param use_full_doc `TRUE` builds each document's vector from all its words,
#'   `FALSE` from the windows around the lexicon's words.
#' @return An object of class `ssd`.
#' @examples
#' \dontrun{
#' s0 <- ssd(df$essay, use_full_doc = FALSE)
#' lex <- suggest_lexicon(s0, "score", data = df)
#' s <- with_lexicon(s0, lexicon_tokens(lex)[1:10])
#' }
#' @export
with_lexicon <- function(s, lexicon, window = 3, sif_a = 1e-3, use_full_doc = FALSE) {
  .check_ssd(s)
  .check_window_sif(window, sif_a)
  .check_flag(use_full_doc, "use_full_doc")
  .check_lexicon(lexicon)
  .new_ssd(s[c("docs", "pre_docs", "tokenizer")], .embeddings_from_ref(s$pack), lexicon,
           window, sif_a, use_full_doc)
}

#' Print a text
#'
#' Shows the number of documents and of documents with a vector, the vector
#' dimension, the lexicon size, the tokenizer and the fit functions.
#'
#' @param x A text from [ssd()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.ssd <- function(x, ...) {
  vectors <- if (is.null(x$x)) "no vectors (tokens only)" else
    sprintf("%s with a vector, D = %d", format(x$n_kept, big.mark = ","), ncol(x$x))
  cat(sprintf("<ssd> %s documents, %s, |L| = %d\n", format(x$n_docs, big.mark = ","), vectors,
              length(x$lexicon)))
  cat("  tokenizer: ", if (is.null(x$tokenizer)) "none (token lists)" else x$tokenizer, "\n", sep = "")
  cat("  fits: fit_pls()  fit_multipls()  fit_regression()  fit_groups()\n")
  invisible(x)
}

# The ssd object from tokens (the docs, pre_docs and tokenizer of
# .tokenize_texts) and settings. Vectors are built when there is a lexicon or
# use_full_doc is TRUE.
.new_ssd <- function(tok, emb, lexicon, window, sif_a, use_full_doc) {
  lexicon <- if (is.null(lexicon)) character(0) else sort(unique(lexicon[!is.na(lexicon)]), method = "radix")
  window <- as.integer(window)
  sif_counts <- x <- keep <- NULL
  if (use_full_doc || length(lexicon)) {
    if (!is.null(sif_a)) sif_counts <- .sif_counts(tok$docs)
    v <- .doc_vectors(tok$docs, emb, lexicon, window, sif_a, use_full_doc, sif_counts)
    x <- v$x
    keep <- v$keep
  }
  structure(list(
    docs = tok$docs,
    pre_docs = tok$pre_docs,
    tokenizer = tok$tokenizer,
    x = x,
    keep = keep,
    n_docs = length(tok$docs),
    n_kept = if (is.null(keep)) NA_integer_ else sum(keep),
    lexicon = lexicon,
    window = window,
    sif_a = sif_a,
    use_full_doc = use_full_doc,
    sif_counts = sif_counts,
    pack = .pack_ref(emb)
  ), class = "ssd")
}

.check_ssd <- function(s, arg = "s") {
  if (!inherits(s, "ssd")) stop(arg, " must be a text from ssd()", call. = FALSE)
}

.check_lexicon <- function(lexicon) {
  if (!is.null(lexicon) && !is.character(lexicon)) {
    stop("lexicon must be a character vector of words", call. = FALSE)
  }
}

.check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(arg, " must be TRUE or FALSE", call. = FALSE)
}

.check_window_sif <- function(window, sif_a) {
  if (!is.numeric(window) || length(window) != 1L || is.na(window) || window < 1 || window != round(window)) {
    stop("window must be >= 1, got ", format(window), call. = FALSE)
  }
  if (!is.null(sif_a) && (!is.numeric(sif_a) || length(sif_a) != 1L || is.na(sif_a) || sif_a <= 0)) {
    stop("sif_a must be > 0, got ", format(sif_a), call. = FALSE)
  }
}

# Stops unless the text has document vectors; every fit calls it.
.require_vectors <- function(s) {
  .check_ssd(s)
  if (is.null(s$x)) {
    stop("This text has no lexicon. Use with_lexicon(s, lexicon) or use_full_doc = TRUE.", call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# The outcome argument
# ---------------------------------------------------------------------------

# What each fit takes: its name in messages, whether a formula is allowed, and
# the form of the values ("numeric": one number per document, "matrix": several
# numeric columns, "labels": group labels as strings, "any": labels as given).
.OUTCOME_KINDS <- list(
  pls = list(method = "fit_pls", formula = TRUE, type = "numeric"),
  regression = list(method = "fit_regression", formula = TRUE, type = "numeric"),
  multipls = list(method = "fit_multipls", formula = FALSE, type = "matrix"),
  groups = list(method = "fit_groups", formula = FALSE, type = "labels"),
  groups_median = list(method = "fit_groups", formula = FALSE, type = "numeric"),
  continuous = list(method = NULL, formula = FALSE, type = "numeric"),
  categorical = list(method = NULL, formula = FALSE, type = "any")
)

# The outcome of a fit over all n documents (Python's _outcome_values and the
# outcome part of SSD._rows). `y` is values, a column name of `data` (several
# for multipls), or a formula whose left side names the column. Returns
# `values` (numeric vector, numeric matrix with column names, or labels),
# `valid` (finite numbers, labels not NA or ""; for a matrix every column),
# `outcome_name`, `outcome_names` (matrix columns) and `formula` (or NULL).
# `check_data = FALSE` skips the row-count and reserved-name checks of `data`,
# as Python's lexicon helpers do.
.outcome_values <- function(y, data, n, kind, method = NULL, check_data = TRUE) {
  k <- .OUTCOME_KINDS[[kind]]
  if (is.null(method)) method <- k$method
  formula <- NULL
  if (inherits(y, "formula")) {
    if (!k$formula) stop(sprintf("%s() takes no formula", method), call. = FALSE)
    if (length(y) != 3L) {
      stop(sprintf("The formula %s needs the outcome column on its left side, e.g. score ~ age.",
                   deparse1(y)), call. = FALSE)
    }
    formula <- y
    y <- deparse1(y[[2L]])
  }
  if (check_data) .check_data(data, n, method)

  name <- "y"
  names <- NULL
  is_name <- is.character(y) && is.null(dim(y)) && (length(y) == 1L || kind == "multipls")
  if (is_name) {
    shown <- paste0("\"", y, "\"", collapse = ", ")
    if (is.null(data)) {
      example <- if (length(y) == 1L) paste0("df$", y) else paste0("df[c(", shown, ")]")
      stop(sprintf("%s() got the column name %s but no data. Pass data = df, or pass the values, e.g. %s.",
                   method, shown, example), call. = FALSE)
    }
    missing <- setdiff(y, names(data))
    if (length(missing)) {
      stop(sprintf("%s(): data has no column %s.", method, paste0("\"", missing, "\"", collapse = ", ")),
           call. = FALSE)
    }
    if (kind == "multipls") {
      names <- y
      y <- data[y]
    } else {
      name <- y
      y <- data[[y]]
    }
  }

  two_d <- is.data.frame(y) || is.matrix(y)
  len <- if (two_d) nrow(y) else length(y)
  if (len != n) {
    stop(sprintf("%s(): the outcome has %d values but the text has %d documents.", method, len, n),
         call. = FALSE)
  }

  if (k$type == "matrix") {
    if (!two_d || ncol(y) < 2L) {
      stop("fit_multipls() needs several outcomes (at least two columns). For one outcome use fit_pls().",
           call. = FALSE)
    }
    cols <- if (is.data.frame(y)) as.list(y) else list(y)
    if (!all(vapply(cols, function(v) is.numeric(v) || is.logical(v), TRUE))) {
      stop("2-D outcome (several outcomes) must be numeric.", call. = FALSE)
    }
    if (is.null(names)) names <- colnames(y)
    if (is.null(names)) names <- paste0("y", seq_len(ncol(y)))
    values <- matrix(as.double(as.matrix(y)), n, ncol(y), dimnames = list(NULL, names))
    valid <- rowSums(!is.finite(values)) == 0
  } else {
    if (two_d) {
      if (ncol(y) == 1L) {
        stop(sprintf("%s() takes one outcome as a vector, but it is a %d x 1 matrix. Pass the column itself, e.g. df$col.",
                     method, n), call. = FALSE)
      }
      stop(sprintf("%s() takes one outcome, but it has %d columns. Use fit_multipls() for several outcomes.",
                   method, ncol(y)), call. = FALSE)
    }
    if (is.factor(y)) y <- as.character(y)
    if (is.logical(y)) y <- as.double(y)
    if (!is.numeric(y) && !is.character(y)) {
      stop(sprintf("%s(): the outcome must be numbers or labels, got %s.", method, class(y)[1]),
           call. = FALSE)
    }
    if (k$type == "numeric" && !is.numeric(y)) {
      stop(switch(kind,
        groups_median = "fit_groups(median_split = TRUE) requires a numeric outcome, but it holds categorical labels.",
        continuous = sprintf("%s() requires a numeric outcome with var_type = \"continuous\", but it holds labels; use var_type = \"categorical\".", method),
        sprintf("%s() requires a numeric outcome, but it holds categorical labels; use fit_groups() instead.", method)
      ), call. = FALSE)
    }
    valid <- if (is.numeric(y)) is.finite(y) else !is.na(y) & y != ""
    values <- if (is.numeric(y)) as.double(y) else y
    if (k$type == "labels" && is.numeric(values)) values <- .number_labels(values)
    values <- unname(values)
  }
  list(values = values, valid = valid, outcome_name = name, outcome_names = names, formula = formula)
}

# Labels of numeric groups, as Python writes them: a whole number as its
# integer ("3", not "3.0"; 1e20 in full), any other number as str() writes it.
.number_labels <- function(v) {
  whole <- is.finite(v) & v == round(v)
  out <- .py_repr(v)
  out[whole] <- sprintf("%.0f", v[whole] + 0)          # + 0 turns -0 into 0
  out
}

# Python's str() of doubles: the fewest digits that read back to the same
# double, in fixed notation for exponents -4 to 15 ("5.0", "0.30000000000000004")
# and in scientific notation otherwise ("1e-05", "1e+20").
.py_repr <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("nan")
    if (is.infinite(v)) return(if (v > 0) "inf" else "-inf")
    if (v == 0) return(if (1 / v < 0) "-0.0" else "0.0")
    for (d in 1:17) {
      s <- sprintf("%.*e", d - 1L, v)
      if (as.numeric(s) == v) break
    }
    m <- regmatches(s, regexec("^(-?)([0-9])\\.?([0-9]*)e([-+][0-9]+)$", s))[[1L]]
    digits <- paste0(m[3L], m[4L])
    e <- as.integer(m[5L])
    if (e < -4L || e >= 16L) {
      mant <- if (nchar(digits) > 1L) paste0(substr(digits, 1L, 1L), ".", substring(digits, 2L)) else digits
      return(sprintf("%s%se%s%02d", m[2L], mant, if (e < 0L) "-" else "+", abs(e)))
    }
    if (e < 0L) return(paste0(m[2L], "0.", strrep("0", -e - 1L), digits))
    digits <- formatC(digits, width = e + 1L, flag = "-")
    digits <- gsub(" ", "0", digits, fixed = TRUE)
    whole <- substr(digits, 1L, e + 1L)
    frac <- substring(digits, e + 2L)
    paste0(m[2L], whole, ".", if (nzchar(frac)) frac else "0")
  }, "", USE.NAMES = FALSE)
}

# Stops unless `data` (if given) is a data frame with n rows and no column
# starting with the reserved prefix ".ssd_".
.check_data <- function(data, n, method) {
  if (is.null(data)) return(invisible())
  if (!is.data.frame(data)) stop("data must be a data frame", call. = FALSE)
  if (nrow(data) != n) {
    stop(sprintf("nrow(data)=%d != s$n_docs=%d. data must have one row per document.", nrow(data), n),
         call. = FALSE)
  }
  clash <- sort(names(data)[startsWith(names(data), ".ssd_")], method = "radix")
  if (length(clash)) {
    stop(sprintf("data columns %s start with '.ssd_', which %s() reserves for its own columns. Rename them.",
                 paste0("\"", clash, "\"", collapse = ", "), method), call. = FALSE)
  }
  invisible()
}

# The rows of one fit (Python's SSD._rows): documents with a valid outcome and
# a vector. Adds to .outcome_values()'s list: `keep` (one per document), `x`
# and `y` (kept rows, in text order), `data` (kept rows, renumbered; NULL
# without data), `n_raw`, `n_kept`, `n_dropped` and `n_no_vector` (documents
# with a valid outcome and no vector).
.outcome_rows <- function(s, y, data, kind) {
  .require_vectors(s)
  out <- .outcome_values(y, data, s$n_docs, kind)
  keep <- out$valid & s$keep
  if (!is.null(data)) {
    data <- data[keep, , drop = FALSE]
    rownames(data) <- NULL
  }
  c(out, list(
    keep = keep,
    x = .x_rows(s, keep),
    y = if (is.matrix(out$values)) out$values[keep, , drop = FALSE] else out$values[keep],
    data = data,
    n_raw = s$n_docs,
    n_kept = sum(keep),
    n_dropped = s$n_docs - sum(keep),
    n_no_vector = sum(out$valid & !s$keep)
  ))
}

# Vectors of the documents where `keep` is TRUE; each must have a vector.
.x_rows <- function(s, keep) s$x[keep[s$keep], , drop = FALSE]

# Stops when a column the formula's right side uses (outside ssd() terms), or
# a fixed-effect expression such as log(x), is missing in the kept rows `data`
# (Python's _check_nobs message). The outcome's own gaps are already out of
# `data`. Python checks the fit's row count afterwards; checking first also
# covers lm()/glm(), which would drop the rows.
.check_nobs <- function(formula, data) {
  vars <- intersect(.vars_outside_ssd(formula[[3L]]), names(data))
  if (!length(vars)) return(invisible())
  fixed <- .split_extras(formula[[3L]])$fixed
  exprs <- if (is.null(fixed)) list() else
    as.list(attr(stats::terms(stats::as.formula(call("~", fixed), env = baseenv())), "variables"))[-1L]
  exprs <- Filter(function(e) is.call(e) && !.is_ssd_call(e), exprs)
  vals <- suppressWarnings(lapply(exprs, eval, data, environment(formula)))
  bad <- sum(!do.call(stats::complete.cases, c(list(data[vars]), vals)))
  if (bad) {
    stop(sprintf(paste0("The formula's columns have missing values in %d of the %d kept documents. ",
                        "Fill or drop those rows before fitting."), bad, nrow(data)), call. = FALSE)
  }
  invisible()
}

# Names a formula or expression uses, leaving out the arguments of ssd() terms.
.vars_outside_ssd <- function(e) {
  if (missing(e)) return(character(0))
  if (is.name(e)) return(as.character(e))
  if (!is.call(e) || identical(e[[1L]], as.name("ssd"))) return(character(0))
  unique(unlist(lapply(as.list(e)[-1L], .vars_outside_ssd)))
}
