# Port of ssdiff/utils/vectors.py. Change the two together.

# SIF word counts over all documents of a text (compute_global_sif): `counts`,
# a named integer vector, and `total`, the number of tokens.
.sif_counts <- function(docs) {
  tokens <- unlist(docs, use.names = FALSE)
  keys <- unique(tokens)
  list(counts = stats::setNames(tabulate(match(tokens, keys), length(keys)), keys),
       total = length(tokens))
}

# SIF weights a / (a + p(word)) of `tokens` under `sif_counts`; 1 for every
# token when `sif_a` is NULL (`sif_counts` is then not evaluated).
.sif_weights <- function(tokens, sif_a, sif_counts) {
  if (is.null(sif_a)) return(rep(1, length(tokens)))
  count <- unname(sif_counts$counts[tokens])
  count[is.na(count)] <- 0
  sif_a / (sif_a + count / max(sif_counts$total, 1))
}

# Document vectors (build_and_normalize_doc_vectors): `x` (one unit row per
# document with a vector) and `keep` (one logical per document). Seed mode
# averages, per lexicon hit, the SIF-weighted context words within `window`
# (the hit excluded), then averages the hits; full mode takes the SIF-weighted
# mean of all words. `sif_a = NULL` weighs every word 1. `sif_counts` are the
# training counts (.sif_counts of the fitted docs).
.doc_vectors <- function(docs, emb, lexicon, window, sif_a, use_full_doc, sif_counts) {
  n <- length(docs)
  tokens <- unlist(docs, use.names = FALSE)
  doc <- rep.int(seq_len(n), lengths(docs))
  wv <- .word_vectors(emb, tokens)
  row <- match(tokens, rownames(wv))
  w <- .sif_weights(tokens, sif_a, sif_counts)

  if (use_full_doc) {
    has <- which(!is.na(row))
    w_sum <- .sum_by(w[has], doc[has], n)
    keep <- w_sum > 0
    sums <- .weighted_row_sums(doc[has], row[has], w[has] / w_sum[doc[has]], n, wv)
  } else {
    hit <- which(tokens %in% lexicon)
    # One entry per (hit, context word with a vector).
    pairs <- lapply(c(-rev(seq_len(window)), seq_len(window)), function(o) {
      j <- hit + o
      ok <- j >= 1L & j <= length(tokens)
      ok[ok] <- doc[j[ok]] == doc[hit[ok]] & !is.na(row[j[ok]])
      cbind(h = seq_along(hit)[ok], j = j[ok])
    })
    pairs <- do.call(rbind, pairs)
    h <- pairs[, "h"]
    j <- pairs[, "j"]
    occ_w <- .sum_by(w[j], h, length(hit))
    occ_doc <- doc[hit]
    n_occ <- tabulate(occ_doc[occ_w > 0], n)
    keep <- n_occ > 0
    sums <- .weighted_row_sums(occ_doc[h], row[j], w[j] / occ_w[h] / n_occ[occ_doc[h]], n, wv)
  }
  x <- sums[keep, , drop = FALSE]
  norms <- pmax(sqrt(rowSums(x^2)), 1e-12)
  list(x = x / norms, keep = keep)
}

# Sums of `v` per group 1..n.
.sum_by <- function(v, g, n) {
  out <- numeric(n)
  if (length(v)) {
    s <- rowsum(v, g)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

# n x D matrix whose row d is the sum of weight * wv[row, ] over the entries
# of document d. Entries are summed in blocks, so memory stays bounded.
.weighted_row_sums <- function(doc, row, weight, n, wv, block = 100000L) {
  out <- matrix(0, n, ncol(wv))
  starts <- if (length(doc)) seq.int(1L, length(doc), by = block) else integer(0)
  for (start in starts) {
    i <- start:min(start + block - 1L, length(doc))
    part <- rowsum(wv[row[i], , drop = FALSE] * weight[i], doc[i])
    d <- as.integer(rownames(part))
    out[d, ] <- out[d, , drop = FALSE] + part
  }
  out
}
