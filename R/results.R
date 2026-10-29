# Port of the views of ssdiff/results/* (single_result.py, continuous_result.py,
# group_result.py, multi_pls_result.py): what each accessor's table holds,
# without the view machinery. Change them together.

# One entry of a result's `texts` (04, "Every result has texts"): the text's
# tokens, settings and pack, and the vectors `x` of the rows `keep` the fit used.
.result_text <- function(s, keep, x, n_no_vector, name = NULL) {
  list(name = name, docs = s$docs, pre_docs = s$pre_docs, tokenizer = s$tokenizer, pack = s$pack,
       lexicon = s$lexicon, window = s$window, sif_a = s$sif_a, use_full_doc = s$use_full_doc,
       x = x, keep = keep, n_no_vector = n_no_vector)
}

# ---------------------------------------------------------------------------
# Units
# ---------------------------------------------------------------------------

# One unit (Python's leaf): its key, gradient, text (position in
# result$texts), the `contrast` its words carry and the one its other tables
# carry, and its kind ("main", "term", "dim", "pair", "contrast").
.view <- function(key, gradient, text, kind, word_contrast = NA_character_, contrast = NA_character_) {
  list(key = key, gradient = gradient, text = text, kind = kind, word_contrast = word_contrast, contrast = contrast)
}

# The units an accessor works on: `views`, and `stack`, the column that tells
# stacked units apart ("contrast", "dim", "term"; NULL for one unit).
.views <- function(result, unit) {
  if (inherits(result, "ssd_contrast")) {
    if (!is.null(unit)) stop("A contrast has one unit; leave unit = NULL.", call. = FALSE)
    return(list(views = list(.view("contrast", result$gradient, 1L, "contrast")), stack = NULL))
  }
  if (inherits(result, "ssd_groups")) return(.views_pairs(result, unit))
  if (inherits(result, "ssd_multipls")) return(.views_dims(result, unit, stack = TRUE))
  if (inherits(result, "ssd_pca_ols") && is.null(result$pls_info) && !is.null(result$units)) {
    return(.views_terms(result, unit))
  }
  if (is.null(unit)) return(list(views = list(.view("main", result$gradient, 1L, "main")), stack = NULL))
  if (is.null(result$units)) stop("This result has one unit; leave unit = NULL.", call. = FALSE)
  .views_dims(result, unit, stack = FALSE)
}

# The key of `unit` among `keys`: a key, or a position.
.unit_key <- function(unit, keys) {
  if (is.numeric(unit) && length(unit) == 1L && !is.na(unit) && unit == round(unit) && unit >= 1 &&
      unit <= length(keys)) {
    return(keys[[unit]])
  }
  if (is.character(unit) && length(unit) == 1L && unit %in% keys) return(unit)
  stop(sprintf("Unknown unit %s; the units are %s, or their positions 1 to %d.", .shown(unit), .quoted(keys),
               length(keys)), call. = FALSE)
}

.views_dims <- function(result, unit, stack) {
  keys <- names(result$units)
  one <- function(key) {
    .view(key, result$units[[key]]$gradient, 1L, "dim", word_contrast = key)
  }
  if (is.null(unit)) return(list(views = lapply(keys, one), stack = "dim"))
  list(views = list(one(.unit_key(unit, keys))), stack = NULL)
}

.views_terms <- function(result, unit) {
  keys <- names(result$units)
  one <- function(key) .view(key, result$units[[key]]$gradient, result$units[[key]]$text, "term")
  if (is.null(unit)) return(list(views = lapply(keys, one), stack = "term"))
  list(views = list(one(.unit_key(unit, keys))), stack = NULL)
}

# A pair as c("g1", "g2") in either order, with canonical names or the
# user's labels, or as its key "g1_g2".
.views_pairs <- function(result, unit) {
  keys <- names(result$units)
  one <- function(key) .view(key, result$units[[key]]$gradient, 1L, "pair", word_contrast = key, contrast = key)
  if (is.null(unit)) return(list(views = lapply(keys, one), stack = if (length(keys) > 1L) "contrast"))
  key <- NULL
  if (is.character(unit) && length(unit) == 1L && unit %in% keys) key <- unit
  if (is.character(unit) && length(unit) == 2L && !anyNA(unit)) {
    gl <- result$group_labels
    clash <- unit[unit %in% names(gl) & unit %in% gl & unit != gl[match(unit, names(gl))]]
    if (length(clash)) {
      stop(sprintf("'%s' is a group's label and the internal name of another group; pass the pair by its key, one of %s.",
                   clash[1L], .quoted(keys)), call. = FALSE)
    }
    canon <- ifelse(unit %in% names(gl), unit, names(gl)[match(unit, gl)])
    if (!anyNA(canon)) {
      canon <- canon[order(as.integer(substring(canon, 2L)))]
      key <- paste(canon, collapse = "_")
    }
  }
  if (is.null(key) || !key %in% keys) {
    shown <- vapply(keys, function(k) {
      u <- result$units[[k]]
      sprintf("c(\"%s\", \"%s\")", u$g1, u$g2)
    }, "")
    stop(sprintf("Unknown pair %s; the pairs are %s.", deparse1(unit), paste(shown, collapse = ", ")), call. = FALSE)
  }
  list(views = list(one(key)), stack = NULL)
}

# Runs `table_of(view)` on each unit and stacks the tables. Stacked units get
# their key in the `stack` column, put first, unless the table has it already.
.stack <- function(v, table_of) {
  tabs <- lapply(v$views, table_of)
  if (is.null(v$stack)) return(tabs[[1L]])
  if (!v$stack %in% names(tabs[[1L]])) {
    tabs <- Map(function(t, view) {
      cls <- class(t)
      t <- cbind(stats::setNames(data.frame(rep(view$key, nrow(t)), stringsAsFactors = FALSE), v$stack), t)
      class(t) <- cls
      t
    }, tabs, v$views)
  }
  out <- do.call(rbind, tabs)
  rownames(out) <- NULL
  class(out) <- class(tabs[[1L]])
  out
}

# Embeddings of a result's texts, each pack read once per accessor call.
.pack_cache <- function(result) {
  cache <- new.env(parent = emptyenv())
  function(text) {
    ref <- result$texts[[text]]$pack
    key <- paste(ref$name, ref$version, ref$abtt)
    if (is.null(cache[[key]])) cache[[key]] <- .embeddings_from_ref(ref)
    cache[[key]]
  }
}

.check_result <- function(result) {
  if (!inherits(result, "ssd_result")) {
    stop("result must be a result of fit_pls(), fit_multipls(), fit_regression(), fit_groups() or contrast()",
         call. = FALSE)
  }
}

.check_side <- function(side) {
  if (!is.character(side) || length(side) != 1L || !side %in% c("pos", "neg", "both")) {
    stop(sprintf("side must be 'pos', 'neg' or 'both', got %s", .shown(side)), call. = FALSE)
  }
  if (side == "both") c("pos", "neg") else side
}

.check_n <- function(n, arg, null_ok = TRUE) {
  if (is.null(n) && null_ok) return(invisible())
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 || (is.finite(n) && n != round(n))) {
    stop(sprintf("%s must be a whole number >= 0%s, got %s", arg, if (null_ok) " or NULL" else "", deparse1(n)),
         call. = FALSE)
  }
}

.check_min_cosine <- function(min_cosine) {
  if (!is.null(min_cosine) && (!is.numeric(min_cosine) || length(min_cosine) != 1L || is.na(min_cosine))) {
    stop("min_cosine must be a number or NULL, got ", deparse1(min_cosine), call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Words and clusters
# ---------------------------------------------------------------------------

#' Words, clusters and snippets of a result
#'
#' Read a fitted gradient in the embeddings and in the texts. Each accessor
#' computes its table when called; the tables are data frames with the
#' columns of Python ssdiff's views, in the same order.
#'
#' - `words()`: the nearest words of each pole, by cosine with the gradient:
#'   `side`, `rank`, `word`, `cos_beta`, `contrast`. Of the 100 nearest
#'   words per side (`n = NULL`), the first `n`.
#' - `clusters()`: k-means clusters of the `topn` nearest words of each pole,
#'   with the k in `k_min`..`k_max` of the best silhouette unless `k` is
#'   given: `cluster_id`, `side`, `size`, `coherence` (the mean cosine of the
#'   words with their centroid), `centroid_cos_beta`, `contrast`, `top_words`
#'   (the five words nearest the centroid). The positive side's clusters come
#'   in descending, the negative side's in ascending order of
#'   `centroid_cos_beta`.
#' - `cluster_words()`: every word of those clusters: `cluster_id`, `side`,
#'   `word`, `cos_centroid`, `cos_beta`, `contrast`.
#' - `snippets()`: text windows scored by cosine with each pole. Each
#'   occurrence of a lexicon word gives one window, the SIF-weighted words
#'   within the text's `window` on each side; without a lexicon each
#'   sentence is one window. The `top_per_side` best per side, at least
#'   `min_cosine`.
#' - `cluster_snippets()`: windows scored against each cluster's centroid,
#'   the `top_per_cluster` best per cluster; `cluster_id` matches
#'   `clusters()` called with the same arguments. Without a lexicon each
#'   document is one window.
#'
#' The snippet tables have `snippet_id`, `side`, `doc_id`, `cosine`, `seed`
#' (the lexicon word), `start_token_idx`, `end_token_idx`, `start_sent_idx`,
#' `end_sent_idx`, `text_window`, `text_surface` (the whole document),
#' `cluster_id` and `contrast`. Printing clips the two text columns to 40
#' characters; the data frame keeps the full text.
#'
#' Neighbour lists leave out the words the pack's neighbour filter matches,
#' and search the 50,000 rows of the pack held in memory. `clusters()`,
#' `cluster_words()` and `cluster_snippets()` with the same arguments give
#' the same clusters: k-means draws from ssdiff's own generator, seeded by
#' `random_state`.
#'
#' A text that came in as token lists has no sentences, so its snippets
#' raise. A [contrast()] has words and clusters only.
#'
#' Ids count from 1: `doc_id` is the document's row in the text, and
#' `cluster_id` the k-means cluster, from 1.
#'
#' @section Units:
#' `unit` picks one unit of a result with several:
#' - [fit_groups()]: a pair, `c("g1", "g2")` in either order (the canonical
#'   names, or the group labels);
#' - [fit_multipls()], and [fit_pls()] or [fit_regression()] with
#'   `ssd_gradients > 1`: a dimension, `"dim-2"` or `2`;
#' - [fit_regression()] with several text terms: a term key such as
#'   `"sexM:ssd(essay)"`, or its position in `names(res$units)`.
#'
#' `unit = NULL` is the fit's own gradient for `fit_pls()` and
#' `fit_regression()` with one text term, and the only pair of a two-group
#' fit. Otherwise it stacks every unit: groups by `contrast`, dimensions
#' with a `dim` column, text terms with a `term` column.
#'
#' The embeddings come from the session (see [load_wordpack()]); after
#' [readRDS()] in a new session the pack is loaded from the cache by its full
#' name, or the accessor stops with a message naming it.
#'
#' @param result A result of [fit_pls()], [fit_multipls()],
#'   [fit_regression()], [fit_groups()] or [contrast()].
#' @param side `"pos"`, `"neg"` or `"both"`.
#' @param n The number of words per side; `NULL` for all 100.
#' @param topn The number of nearest words per side that are clustered.
#' @param k The number of clusters, or `NULL` to choose it by silhouette.
#' @param k_min,k_max The range of k tried when `k` is `NULL`.
#' @param min_cluster_size Clusters with fewer words are left out.
#' @param random_state Seed of k-means, a whole number from 0 to 2^53 - 1.
#' @param top_per_side The number of snippets per side; `NULL` for all.
#' @param top_per_cluster The number of snippets per cluster.
#' @param min_cosine The lowest cosine a snippet may have, or `NULL`.
#' @param unit The unit, see "Units".
#' @return A data frame.
#' @seealso [docs()], [doc_text()], [outcomes()]
#' @examples
#' \dontrun{
#' r <- fit_pls(s, "score", data = df)
#' words(r, n = 10)
#' clusters(r, side = "pos")
#' snippets(r, top_per_side = 5)
#' g <- fit_groups(s, "condition", data = df)
#' words(g, unit = c("g1", "g3"))
#' }
#' @name accessors
NULL

#' @rdname accessors
#' @export
words <- function(result, side = "both", n = 20, unit = NULL) {
  .check_result(result)
  sides <- .check_side(side)
  .check_n(n, "n")
  emb_of <- .pack_cache(result)
  v <- .views(result, unit)
  .stack(v, function(view) .words_table(emb_of(view$text), view$gradient, view$word_contrast, sides, n))
}

# The 100 nearest words of each side, cosine >= 0, the first n (NULL: all).
.words_table <- function(emb, gradient, contrast, sides, n) {
  tabs <- lapply(sides, function(side) {
    nb <- .filtered_neighbors(emb, if (side == "pos") gradient else -gradient, topn = 100L)
    nb <- nb[nb$cos >= 0, , drop = FALSE]
    if (!is.null(n)) nb <- utils::head(nb, n)
    k <- nrow(nb)
    data.frame(side = rep(side, k), rank = seq_len(k), word = nb$word,
               cos_beta = if (side == "pos") nb$cos else -nb$cos, contrast = rep(contrast, k),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, tabs)
  rownames(out) <- NULL
  out
}

.cluster_args <- function(topn, k, k_min, k_max, min_cluster_size, random_state) {
  .check_n(topn, "topn", null_ok = FALSE)
  if (!is.null(k)) .check_count(k, "k")
  .check_count(k_min, "k_min")
  .check_count(k_max, "k_max")
  .check_n(min_cluster_size, "min_cluster_size", null_ok = FALSE)
  .rng(random_state)                                   # checks the seed
  list(topn = topn, k = k, k_min = k_min, k_max = k_max, min_cluster_size = min_cluster_size,
       random_state = random_state)
}

# .cluster_top_neighbors() of one side with the accessors' arguments.
.side_clusters <- function(emb, gradient, side, a) {
  .cluster_top_neighbors(emb, gradient, topn = a$topn, k = a$k, k_min = a$k_min, k_max = a$k_max,
                         restrict_vocab = .RAM_TOP_N, random_state = a$random_state,
                         min_cluster_size = a$min_cluster_size, side = side)
}

.clusters_table <- function(cl, side, contrast) {
  k <- length(cl)
  data.frame(cluster_id = vapply(cl, function(c) as.integer(c$id), 0L), side = rep(side, k),
             size = vapply(cl, function(c) as.integer(c$size), 0L), coherence = vapply(cl, `[[`, 0, "coherence"),
             centroid_cos_beta = vapply(cl, `[[`, 0, "centroid_cos_beta"), contrast = rep(contrast, k),
             top_words = vapply(cl, function(c) paste(utils::head(c$words$word, 5L), collapse = ", "), ""),
             stringsAsFactors = FALSE)
}

.cluster_words_table <- function(cl, side, contrast) {
  tabs <- lapply(cl, function(c) {
    w <- c$words
    data.frame(cluster_id = rep(as.integer(c$id), nrow(w)), side = rep(side, nrow(w)), word = w$word,
               cos_centroid = w$cos_centroid, cos_beta = w$cos_beta, contrast = rep(contrast, nrow(w)),
               stringsAsFactors = FALSE)
  })
  if (!length(tabs)) {
    return(data.frame(cluster_id = integer(0), side = character(0), word = character(0), cos_centroid = numeric(0),
                      cos_beta = numeric(0), contrast = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, tabs)
}

# One clusters-type table per unit and side, stacked.
.cluster_accessor <- function(result, side, a, unit, table_of) {
  .check_result(result)
  sides <- .check_side(side)
  emb_of <- .pack_cache(result)
  v <- .views(result, unit)
  .stack(v, function(view) {
    emb <- emb_of(view$text)
    out <- do.call(rbind, lapply(sides, function(sd) table_of(.side_clusters(emb, view$gradient, sd, a), sd, view$contrast)))
    rownames(out) <- NULL
    out
  })
}

#' @rdname accessors
#' @export
clusters <- function(result, side = "both", topn = 100, k = NULL, k_min = 2, k_max = 10, min_cluster_size = 2,
                     random_state = 2137, unit = NULL) {
  a <- .cluster_args(topn, k, k_min, k_max, min_cluster_size, random_state)
  .cluster_accessor(result, side, a, unit, .clusters_table)
}

#' @rdname accessors
#' @export
cluster_words <- function(result, side = "both", topn = 100, k = NULL, k_min = 2, k_max = 10,
                          min_cluster_size = 2, random_state = 2137, unit = NULL) {
  a <- .cluster_args(topn, k, k_min, k_max, min_cluster_size, random_state)
  .cluster_accessor(result, side, a, unit, .cluster_words_table)
}

# ---------------------------------------------------------------------------
# Snippets
# ---------------------------------------------------------------------------

# The text a unit's snippets come from; stops for a contrast and token lists.
.snippet_text <- function(result, view) {
  if (view$kind == "contrast") stop("A contrast has no documents, so it has no snippets.", call. = FALSE)
  t <- result$texts[[view$text]]
  if (is.null(t$pre_docs)) stop("No snippets: this text came in as token lists", call. = FALSE)
  t
}

# Occurrences of a result's texts, each built once per accessor call.
.occurrence_cache <- function(result, emb_of, mode) {
  cache <- new.env(parent = emptyenv())
  function(view) {
    key <- as.character(view$text)
    if (is.null(cache[[key]])) {
      t <- .snippet_text(result, view)
      tk <- .snippet_tokens(t, emb_of(view$text))
      occ <- if (length(t$lexicon)) .seed_occurrences(tk, t$lexicon, t$window) else mode(tk)
      cache[[key]] <- c(occ, list(essay = tk$essay))
    }
    cache[[key]]
  }
}

#' @rdname accessors
#' @export
snippets <- function(result, side = "both", top_per_side = 30, min_cosine = NULL, unit = NULL) {
  .check_result(result)
  sides <- .check_side(side)
  .check_n(top_per_side, "top_per_side")
  .check_min_cosine(min_cosine)
  emb_of <- .pack_cache(result)
  occ_of <- .occurrence_cache(result, emb_of, .sentence_occurrences)
  v <- .views(result, unit)
  .stack(v, function(view) {
    occ <- occ_of(view)
    .snippets_along(occ, occ$essay, view$gradient, sides, top_per_side, min_cosine, view$contrast)
  })
}

#' @rdname accessors
#' @export
cluster_snippets <- function(result, side = "both", top_per_cluster = 100, min_cosine = NULL, topn = 100, k = NULL,
                             k_min = 2, k_max = 10, min_cluster_size = 2, random_state = 2137, unit = NULL) {
  .check_result(result)
  sides <- .check_side(side)
  .check_n(top_per_cluster, "top_per_cluster", null_ok = FALSE)
  .check_min_cosine(min_cosine)
  a <- .cluster_args(topn, k, k_min, k_max, min_cluster_size, random_state)
  emb_of <- .pack_cache(result)
  occ_of <- .occurrence_cache(result, emb_of, .document_occurrences)
  v <- .views(result, unit)
  .stack(v, function(view) {
    occ <- occ_of(view)
    emb <- emb_of(view$text)
    out <- do.call(rbind, lapply(sides, function(sd) {
      .snippets_by_centroids(occ, occ$essay, emb, .side_clusters(emb, view$gradient, sd, a), sd, top_per_cluster,
                             min_cosine, view$contrast)
    }))
    rownames(out) <- NULL
    structure(out, class = c("ssd_snippets", "data.frame"))
  })
}

# ---------------------------------------------------------------------------
# Documents
# ---------------------------------------------------------------------------

# Cosines of the rows of x with a unit gradient.
.alignment <- function(x, gradient) drop(.l2_normalize_rows(x) %*% gradient)

# A unit's documents: `doc_id`, `y_true`, `y_hat`, `residual`,
# `alignment_score` for a fit's own gradient and its text terms; `doc_id`,
# `group`, `alignment_score` for a pair; `doc_id`, `alignment_score` for a
# dimension.
.docs_table <- function(result, view) {
  if (view$kind == "contrast") stop("A contrast has no documents.", call. = FALSE)
  t <- result$texts[[view$text]]
  doc_id <- which(t$keep)
  if (view$kind == "pair") {
    u <- result$units[[view$key]]
    in_pair <- result$groups %in% c(u$g1, u$g2)
    return(data.frame(doc_id = doc_id[in_pair], group = unname(result$group_labels[result$groups[in_pair]]),
                      alignment_score = u$alignment_scores, stringsAsFactors = FALSE))
  }
  if (view$kind == "dim") {
    return(data.frame(doc_id = doc_id, alignment_score = .alignment(t$x, view$gradient)))
  }
  if (inherits(result, "ssd_pls")) {
    residual <- result$residuals
    align <- result$alignment_scores
  } else {
    residual <- unname(.reg_resid(result$model))
    align <- if (view$kind == "term") result$units[[view$key]]$alignment_scores else result$alignment_scores
  }
  data.frame(doc_id = doc_id, y_true = result$y, y_hat = result$y_hat, residual = unname(residual),
             alignment_score = align)
}

#' Documents of a result
#'
#' `docs()` gives one row per document a unit was fitted on. For
#' [fit_pls()] and [fit_regression()]: `doc_id`, `y_true`, `y_hat`,
#' `residual` and `alignment_score` (the cosine of the document's vector with
#' the gradient). `y_hat` is the text's part of the prediction on the
#' outcome's scale, or the model's fitted mean (on the response scale) for a
#' model with several text terms, an interaction or a family other than the
#' Gaussian; `residual` is `y_true` minus the model's fitted value, or for
#' `fit_pls()` with a formula the PLS residual with the covariates taken
#' out. For a pair of [fit_groups()]: `doc_id`, `group` (the label) and
#' `alignment_score`, the pair's documents only. For a dimension
#' ([fit_multipls()], `ssd_gradients > 1`): `doc_id` and `alignment_score`.
#'
#' `which = "pos"` and `"neg"` give the `n` documents with the highest and
#' the lowest `y_hat` (the alignment score where there is no `y_hat`).
#' `which = "misdiagnosed"` gives the `n` largest residuals: `direction =
#' "both"` by absolute value, `"over"` the most negative (the model predicts
#' too high), `"under"` the most positive.
#'
#' `doc_text()` gives the text of documents of the unit: the raw string, or
#' for a text from token lists the tokens joined by spaces.
#'
#' `doc_id` is the document's row in the text, from 1.
#'
#' @inheritParams accessors
#' @param which `"all"`, `"pos"`, `"neg"` or `"misdiagnosed"`.
#' @param n The number of documents for `which` other than `"all"`.
#' @param direction For `"misdiagnosed"`: `"both"`, `"over"` or `"under"`.
#' @param doc_id Document ids (rows of the text) among the unit's documents.
#' @return `docs()`: a data frame. `doc_text()`: a character vector.
#' @seealso [words()]
#' @examples
#' \dontrun{
#' r <- fit_pls(s, "score", data = df)
#' docs(r, "pos", n = 5)
#' docs(r, "misdiagnosed", direction = "over")
#' doc_text(r, 12)
#' }
#' @export
docs <- function(result, which = "all", n = 20, direction = "both", unit = NULL) {
  .check_result(result)
  if (!is.character(which) || length(which) != 1L || !which %in% c("all", "pos", "neg", "misdiagnosed")) {
    stop(sprintf("which must be 'all', 'pos', 'neg' or 'misdiagnosed', got %s", .shown(which)), call. = FALSE)
  }
  if (!is.character(direction) || length(direction) != 1L || !direction %in% c("both", "over", "under")) {
    stop(sprintf("direction must be one of 'both'/'over'/'under', got %s", .shown(direction)), call. = FALSE)
  }
  .check_n(n, "n", null_ok = FALSE)
  v <- .views(result, unit)
  .stack(v, function(view) {
    d <- .docs_table(result, view)
    if (which == "all") return(d)
    score <- if ("y_hat" %in% names(d)) d$y_hat else d$alignment_score
    if (which == "misdiagnosed") {
      if (!"residual" %in% names(d)) {
        stop("docs(which = \"misdiagnosed\") needs residuals, which only fit_pls() and fit_regression() have.",
             call. = FALSE)
      }
      score <- switch(direction, both = abs(d$residual), over = -d$residual, under = d$residual)
    }
    o <- order(if (which == "neg") score else -score, method = "radix")
    out <- d[utils::head(o, n), , drop = FALSE]
    rownames(out) <- NULL
    out
  })
}

#' @rdname docs
#' @export
doc_text <- function(result, doc_id, unit = NULL) {
  .check_result(result)
  v <- .views(result, unit)
  view <- v$views[[1L]]
  if (view$kind == "contrast") stop("A contrast has no documents.", call. = FALSE)
  texts <- unique(vapply(v$views, `[[`, 0L, "text"))
  if (length(texts) > 1L) {
    stop("doc_text(): the units come from several texts; pass unit = to pick one.", call. = FALSE)
  }
  ids <- unique(unlist(lapply(v$views, function(view) .docs_table(result, view)$doc_id)))
  if (!is.numeric(doc_id) || anyNA(doc_id) || any(!doc_id %in% ids)) {
    bad <- if (is.numeric(doc_id)) doc_id[is.na(doc_id) | !doc_id %in% ids][1L] else deparse1(doc_id)
    stop(sprintf("No doc with doc_id=%s", format(bad)), call. = FALSE)
  }
  t <- result$texts[[texts]]
  if (is.null(t$pre_docs)) {
    return(vapply(t$docs[doc_id], paste, "", collapse = " "))
  }
  vapply(t$pre_docs[doc_id], `[[`, "", "raw")
}

# ---------------------------------------------------------------------------
# Outcomes and the sweep
# ---------------------------------------------------------------------------

#' Outcome loadings of a multi-outcome fit
#'
#' The loadings of the outcomes on each dimension of [fit_multipls()]
#' (`V_rot`): `outcome`, `group` (the outcome's group of `outcome_groups`, or
#' `NA`) and `loading`. With `outcome_groups`, one row per group follows,
#' with the plain mean of its outcomes' loadings (outcomes scored in opposite
#' directions cancel).
#'
#' @param result A result of [fit_multipls()].
#' @param unit A dimension, `"dim-2"` or `2`; `NULL` stacks every dimension
#'   with a `dim` column.
#' @return A data frame.
#' @examples
#' \dontrun{
#' r <- fit_multipls(s, c("valence", "arousal"), data = df)
#' outcomes(r, "dim-1")
#' }
#' @export
outcomes <- function(result, unit = NULL) {
  if (!inherits(result, "ssd_multipls")) stop("outcomes() needs a result of fit_multipls().", call. = FALSE)
  v <- .views(result, unit)
  groups <- result$outcome_groups
  group_of <- stats::setNames(rep(names(groups), lengths(groups)), unlist(groups, use.names = FALSE))
  .stack(v, function(view) {
    loading <- result$V_rot[, match(view$key, names(result$units))]
    names(loading) <- result$outcome_names
    own <- if (length(groups)) unname(group_of[result$outcome_names]) else rep(NA_character_, length(loading))
    data.frame(outcome = c(result$outcome_names, names(groups)), group = c(own, names(groups)),
               loading = c(unname(loading), vapply(groups, function(m) mean(loading[m]), 0, USE.NAMES = FALSE)),
               stringsAsFactors = FALSE)
  })
}

# The position of a text of an ssd_pca_ols: a name or a position; NULL with one text.
.text_index <- function(result, text) {
  if (!inherits(result, "ssd_pca_ols")) stop("This needs a result of fit_regression().", call. = FALSE)
  n <- length(result$texts)
  names <- vapply(result$texts, function(t) if (is.null(t$name)) NA_character_ else t$name, "")
  if (is.null(text)) {
    if (n == 1L) return(1L)
    stop(sprintf("The fit has %d texts; pass text = one of %s, or a position.", n, .quoted(names)), call. = FALSE)
  }
  if (is.numeric(text) && length(text) == 1L && !is.na(text) && text == round(text) && text >= 1 && text <= n) {
    return(as.integer(text))
  }
  if (is.character(text) && length(text) == 1L && text %in% names) return(match(text, names))
  named <- names[!is.na(names)]
  where <- if (length(named)) {
    sprintf("the texts are %s, or their positions 1 to %d", .quoted(named), n)
  } else if (n == 1L) {
    "the fit has one text, at position 1"
  } else {
    sprintf("the texts are at positions 1 to %d", n)
  }
  stop(sprintf("Unknown text %s; %s.", .shown(text), where), call. = FALSE)
}

#' The K sweep of a regression
#'
#' `sweep_table()` gives the per-K table of [fit_regression()]'s sweep for
#' one text: `k`, `var_explained`, the interpretability columns, the
#' gradient's change from the previous K (`beta_delta_1_minus_cos`), their
#' z-scores and window means, and `joint_score`, which picked K. It has no
#' rows when K came from `fixed_k`.
#'
#' `plot_sweep()` draws the chart of the sweep, in base graphics: the
#' detrended interpretability (z) against K on the left axis, the change of
#' the gradient (1 - cosine, a running median over 7 values of K) on the
#' right axis, and the chosen K as a vertical line. With `path` it writes a
#' PNG file at `dpi`, 6.4 by 4.8 inches; otherwise it draws on the current
#' device.
#'
#' @param result A result of [fit_regression()].
#' @param text A text's name or position; `NULL` when the fit has one text.
#' @param path A file path for a PNG, or `NULL`.
#' @param dpi The resolution of the PNG.
#' @return `sweep_table()`: a data frame. `plot_sweep()`: `NULL`, invisibly.
#' @examples
#' \dontrun{
#' r <- fit_regression(s, "score", data = df)
#' sweep_table(r)
#' plot_sweep(r)
#' plot_sweep(r, path = "sweep.png")
#' }
#' @export
sweep_table <- function(result, text = NULL) {
  i <- .text_index(result, text)
  tab <- result$texts[[i]]$sweep
  if (is.null(tab)) {
    tab <- as.data.frame(stats::setNames(lapply(.SWEEP_COLUMNS, function(c) numeric(0)), .SWEEP_COLUMNS))
    for (col in c("k", "n_clusters", "total_size")) tab[[col]] <- integer(0)
  }
  tab
}

#' @rdname sweep_table
#' @export
plot_sweep <- function(result, text = NULL, path = NULL, dpi = 300) {
  i <- .text_index(result, text)
  tab <- result$texts[[i]]$sweep
  if (is.null(tab)) {
    stop("No sweep data: fit_regression() was called with fixed_k. Refit with fixed_k = NULL to run the sweep.",
         call. = FALSE)
  }
  best <- .text_fit_info(result, i)$best_k
  if (!is.null(path)) {
    grDevices::png(path, width = 6.4 * dpi, height = 4.8 * dpi, res = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  blue <- "#1f77b4"
  orange <- "#ff7f0e"
  op <- graphics::par(mar = c(4.5, 4.5, 1, 4.5))
  on.exit(graphics::par(op), add = TRUE, after = FALSE)  # before dev.off(), on the PNG device
  # A sweep of one or two K has no detrended score and no smoothed change: empty axes.
  lim <- function(v) if (any(is.finite(v))) range(v, finite = TRUE) else c(-1, 1)
  graphics::plot(tab$k, tab$interp_resid_z, type = "o", pch = 16, col = blue, xlab = "K", ylab = "",
                 ylim = lim(tab$interp_resid_z))
  graphics::mtext("Detrended interpretability (z)", side = 2, line = 3, col = blue)
  graphics::abline(h = 0, col = "grey60")
  graphics::abline(v = best, col = "red", lwd = 2)
  graphics::par(new = TRUE)
  smooth <- .rolling_median(tab$beta_delta_1_minus_cos, 7L)
  graphics::plot(tab$k, smooth, type = "l", lwd = 2, col = orange, axes = FALSE, xlab = "", ylab = "",
                 ylim = lim(smooth))
  graphics::axis(4, col.axis = orange)
  graphics::mtext("Beta change (smoothed 1 - cosine)", side = 4, line = 3, col = orange)
  invisible(NULL)
}

# The fit_info of the units of text i (the result's own with one text).
.text_fit_info <- function(result, i) {
  for (u in result$units) if (identical(u$text, i) && !is.null(u$fit_info)) return(u$fit_info)
  result$fit_info
}

# The median of the finite values within window %/% 2 positions of each position.
.rolling_median <- function(x, window = 7L) {
  half <- window %/% 2L
  n <- length(x)
  vapply(seq_len(n), function(i) {
    w <- x[max(1L, i - half):min(n, i + half)]
    w <- w[is.finite(w)]
    if (length(w)) stats::median(w) else NA_real_
  }, 0)
}
