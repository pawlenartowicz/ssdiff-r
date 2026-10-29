# Port of ssdiff/utils/snippets.py. Change the two together. Single-threaded;
# positions (tokens, sentences, documents, snippets, clusters) count from 1.

# The tokens of a text from strings, with what the snippets need: per token its
# document `doc`, position `pos` in the document, SIF weight `w`, sentence
# `sent` and the row `row` of its unit vector in `wv` (NA when the pack has no
# vector for it); per document its sentences (`surf`, flat, from `first` + 1)
# and the whole text (`essay`). SIF counts are over all documents of the text.
.snippet_tokens <- function(t, emb) {
  docs <- t$docs
  pre <- t$pre_docs
  tokens <- unlist(docs, use.names = FALSE)
  wv <- .word_vectors(emb, tokens)
  w <- .sif_weights(tokens, t$sif_a, .sif_counts(docs))
  sents <- lapply(pre, `[[`, "sents_surface")
  n_sents <- lengths(sents)
  list(tokens = tokens, doc = rep.int(seq_along(docs), lengths(docs)), pos = sequence(lengths(docs)),
       w = w, wv = wv, row = match(tokens, rownames(wv)),
       sent = unlist(lapply(pre, `[[`, "token_to_sent"), use.names = FALSE),
       n_tokens = lengths(docs), n_sents = n_sents, first = cumsum(c(0L, n_sents))[seq_along(docs)],
       surf = unlist(sents, use.names = FALSE),
       essay = vapply(sents, paste, "", collapse = " "))
}

# Unit vectors of groups of tokens: row g is the SIF-weighted sum of the
# vectors of the tokens `idx` with `group == g`, scaled to unit length. Groups
# whose sum is below 1e-12 are dropped; `keep` says which stayed.
.occurrence_vectors <- function(tk, idx, group, n) {
  has <- !is.na(tk$row[idx])
  sums <- .weighted_row_sums(group[has], tk$row[idx][has], tk$w[idx][has], n, tk$wv)
  norms <- sqrt(rowSums(sums^2))
  keep <- norms > 1e-12
  list(x = sums[keep, , drop = FALSE] / norms[keep], keep = keep)
}

# One occurrence per lexicon hit: the context words within `window` of it in
# its document, the hit left out (Python's _collect_occurrences_for_doc).
# Returns the occurrences' unit vectors `x` and a data frame `meta` of doc_id,
# seed, start/end token and sentence, and the text window.
.seed_occurrences <- function(tk, lexicon, window) {
  hit <- which(tk$tokens %in% lexicon)
  offsets <- c(-rev(seq_len(window)), seq_len(window))
  h <- rep(seq_along(hit), each = length(offsets))
  j <- rep(hit, each = length(offsets)) + rep(offsets, length(hit))
  ok <- j >= 1L & j <= length(tk$tokens)
  ok[ok] <- tk$doc[j[ok]] == tk$doc[hit[h[ok]]]
  occ <- .occurrence_vectors(tk, j[ok], h[ok], length(hit))
  hit <- hit[occ$keep]

  d <- tk$doc[hit]
  i <- tk$pos[hit]
  lo <- pmax(1L, i - window)
  hi <- pmin(tk$n_tokens[d], i + window)
  s <- tk$sent[hit]
  s_lo <- tk$sent[hit - (i - lo)]
  s_hi <- tk$sent[hit + (hi - i)]
  # _make_snippet_anchor: the hit's sentence, joined with the previous one
  # when the window reaches back into it, else with the next one when there is one.
  back <- s_lo < s
  fwd <- !back & s_hi != s & s < tk$n_sents[d]
  start <- ifelse(back, s - 1L, s)
  end <- ifelse(fwd, s + 1L, s)
  base <- tk$first[d]
  text <- tk$surf[base + s]
  text[back] <- trimws(paste(tk$surf[(base + s - 1L)[back]], text[back]), whitespace = "[\\h\\v]")
  text[fwd] <- trimws(paste(text[fwd], tk$surf[(base + s + 1L)[fwd]]), whitespace = "[\\h\\v]")
  list(x = occ$x, meta = data.frame(doc_id = d, seed = tk$tokens[hit], start_token_idx = lo, end_token_idx = i + window,
                                    start_sent_idx = start, end_sent_idx = end, text_window = text,
                                    stringsAsFactors = FALSE))
}

# Without a lexicon, snippets() takes each sentence as one occurrence
# (Python's _collect_sentence_occurrences_for_doc).
.sentence_occurrences <- function(tk) {
  key <- (tk$doc - 1) * (max(tk$sent, 0L) + 1) + tk$sent
  keys <- sort(unique(key))
  g <- match(key, keys)
  occ <- .occurrence_vectors(tk, seq_along(tk$tokens), g, length(keys))
  first_tok <- match(seq_along(keys), g)
  d <- tk$doc[first_tok][occ$keep]
  s <- tk$sent[first_tok][occ$keep]
  pos <- split(tk$pos, g)[occ$keep]
  text <- ifelse(s <= tk$n_sents[d], tk$surf[tk$first[d] + s], "")
  list(x = occ$x, meta = data.frame(doc_id = d, seed = rep("<SENT>", length(d)),
                                    start_token_idx = vapply(pos, min, 0L), end_token_idx = vapply(pos, max, 0L),
                                    start_sent_idx = s, end_sent_idx = s, text_window = text,
                                    stringsAsFactors = FALSE))
}

# Without a lexicon, cluster_snippets() takes each document as one occurrence
# (Python's _collect_doc_occurrences_for_doc).
.document_occurrences <- function(tk) {
  n <- length(tk$n_tokens)
  occ <- .occurrence_vectors(tk, seq_along(tk$tokens), tk$doc, n)
  d <- which(occ$keep)
  list(x = occ$x, meta = data.frame(doc_id = d, seed = rep("<DOC>", length(d)), start_token_idx = rep(1L, length(d)),
                                    end_token_idx = tk$n_tokens[d], start_sent_idx = rep(1L, length(d)),
                                    end_sent_idx = pmax(1L, tk$n_sents[d]), text_window = tk$essay[d],
                                    stringsAsFactors = FALSE))
}

# A snippet table: Python's Snippet columns without text_lemmas.
.snippet_table <- function(meta, rows, side, cosine, cluster_id, contrast, essay) {
  m <- meta[rows, , drop = FALSE]
  n <- length(rows)
  structure(data.frame(snippet_id = seq_len(n), side = rep(side, n), doc_id = m$doc_id, cosine = cosine,
                       seed = m$seed, start_token_idx = m$start_token_idx, end_token_idx = m$end_token_idx,
                       start_sent_idx = m$start_sent_idx, end_sent_idx = m$end_sent_idx,
                       text_window = m$text_window, text_surface = essay[m$doc_id],
                       cluster_id = rep(as.integer(cluster_id), length.out = n), contrast = rep(contrast, n),
                       stringsAsFactors = FALSE),
            class = c("ssd_snippets", "data.frame"))
}

# snippets_along_beta(): the occurrences scored by cosine with +gradient
# ("pos") and -gradient ("neg"), at least `min_cosine`, the `top_per_side`
# highest per side (NULL: all). snippet_id runs over the sides in order.
.snippets_along <- function(occ, essay, gradient, sides, top_per_side, min_cosine, contrast) {
  cos_pos <- drop(occ$x %*% .unit_vector(gradient))
  out <- lapply(sides, function(side) {
    cos <- if (side == "pos") cos_pos else -cos_pos
    idx <- if (is.null(min_cosine)) seq_along(cos) else which(cos >= min_cosine)
    idx <- idx[order(-cos[idx], method = "radix")]
    if (!is.null(top_per_side)) idx <- utils::head(idx, top_per_side)
    .snippet_table(occ$meta, idx, side, cos[idx], NA_integer_, contrast, essay)
  })
  .renumber(do.call(rbind, out))
}

# Stacked snippet tables numbered again from 1.
.renumber <- function(out) {
  out$snippet_id <- seq_len(nrow(out))
  rownames(out) <- NULL
  structure(out, class = c("ssd_snippets", "data.frame"))
}

# cluster_snippets_by_centroids() for one side: each cluster's centroid (the
# unit mean of its words' vectors) scores every occurrence; per cluster the
# `top_per_cluster` best, the clusters in the order of Python's labels
# "<side>_cluster_<rank>" as strings; then `min_cosine`. `clusters` are
# .cluster_top_neighbors()'s, in rank order.
.snippets_by_centroids <- function(occ, essay, emb, clusters, side, top_per_cluster, min_cosine, contrast) {
  empty <- .snippet_table(occ$meta, integer(0), side, numeric(0), NA_integer_, contrast, essay)
  targets <- lapply(clusters, function(cl) {
    v <- .word_vectors(emb, cl$words$word)
    if (nrow(v)) .unit_vector(colMeans(v)) else numeric(0)
  })
  ok <- which(vapply(targets, function(v) length(v) && any(v != 0), TRUE))
  if (!length(ok) || !nrow(occ$x)) return(empty)
  scores <- occ$x %*% do.call(cbind, targets[ok])
  labels <- paste0(side, "_cluster_", ok)
  parts <- lapply(order(labels, method = "radix"), function(k) {
    cos <- scores[, k]
    idx <- utils::head(order(-cos, method = "radix"), top_per_cluster)
    .snippet_table(occ$meta, idx, side, cos[idx], clusters[[ok[k]]]$id, contrast, essay)
  })
  out <- do.call(rbind, parts)
  if (!is.null(min_cosine)) out <- out[out$cosine >= min_cosine, , drop = FALSE]
  .renumber(out)
}

# Printing clips the text columns to 40 characters; the data frame keeps them whole.
#' @rdname result-print
#' @export
print.ssd_snippets <- function(x, ...) {
  clip <- function(v) ifelse(!is.na(v) & nchar(v) > 40L, paste0(substr(v, 1L, 39L), "\u2026"), v)
  shown <- as.data.frame(x)
  for (col in c("text_window", "text_surface")) if (col %in% names(shown)) shown[[col]] <- clip(shown[[col]])
  print(shown, ...)
  invisible(x)
}
