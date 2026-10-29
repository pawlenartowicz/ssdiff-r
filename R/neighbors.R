# Port of ssdiff/utils/neighbors.py and Embeddings.similar_by_vector in
# ssdiff/embeddings.py. Change them together.

# TRUE for words hidden from neighbour lists: words the pack's
# neighbor_filter matches at their start (Python's re.match), or without a
# filter, words with a digit or starting with an uppercase letter.
.skip_neighbor <- function(emb, words) {
  nf <- emb$neighbor_filter
  if (!is.null(nf)) return(grepl(paste0("(*UCP)^(?:", nf, ")"), words, perl = TRUE))
  grepl("(*UCP)\\d|^\\p{Lu}", words, perl = TRUE)
}

# Keys of 1-based rows of the pack.
.row_keys <- function(emb, rows) vapply(rows - 1L, function(r) .wp_key(emb$pack$vocab, r), "")

# The `topn` in-memory rows most similar to `vec`, by cosine: a data frame
# with `row` (1-based) and `cos`, most similar first. Only the first
# `restrict_vocab` rows are searched, and never more than the rows in RAM.
.similar_by_vector <- function(emb, vec, topn = 10L, restrict_vocab = NULL) {
  held <- nrow(emb$vectors)
  if (is.null(restrict_vocab)) restrict_vocab <- held
  if (restrict_vocab > held && emb$pack$n > held) {
    warning(sprintf(paste0("restrict_vocab=%d clamped to %d: only the %d most frequent wordpack rows ",
                           "are held in RAM and searchable."), restrict_vocab, held, held), call. = FALSE)
  }
  restrict_vocab <- min(restrict_vocab, held)
  norm <- sqrt(sum(vec^2))
  if (norm < 1e-12 || restrict_vocab < 1L) return(data.frame(row = integer(0), cos = numeric(0)))
  sims <- drop(emb$vectors[seq_len(restrict_vocab), , drop = FALSE] %*% (vec / norm))
  top <- order(-sims, method = "radix")[seq_len(min(topn, length(sims)))]
  data.frame(row = top, cos = sims[top])
}

# The `topn` nearest words to `vec` that .skip_neighbor() keeps, from the
# `cand` nearest of the first `restrict` rows: a data frame with `word`, `cos`
# and `row` (1-based), most similar first.
.filtered_neighbors <- function(emb, vec, topn = 20L, cand = 2000L, restrict = 10000L) {
  nb <- .similar_by_vector(emb, vec, topn = cand, restrict_vocab = restrict)
  words <- .row_keys(emb, nb$row)
  ok <- which(!.skip_neighbor(emb, words))
  ok <- ok[seq_len(min(topn, length(ok)))]
  data.frame(word = words[ok], cos = nb$cos[ok], row = nb$row[ok], stringsAsFactors = FALSE)
}

# Clusters of the `topn` filtered neighbours of +beta (side "pos") or -beta
# ("neg"): k-means with a fixed `k`, or the k in k_min..k_max with the best
# silhouette. Clusters smaller than min_cluster_size are dropped. Returns a
# list of clusters, each with `id` (the k-means label, from 1), `size`,
# `centroid_cos_beta`, `coherence` and `words` (a data frame of `word`,
# `cos_centroid`, `cos_beta`, by cos_centroid descending). Clusters are sorted
# by centroid_cos_beta, descending for "pos" and ascending for "neg".
.cluster_top_neighbors <- function(emb, beta, topn = 100L, k = NULL, k_min = 2L, k_max = 10L,
                                   restrict_vocab = .RAM_TOP_N, random_state = 2137,
                                   min_cluster_size = 2L, side = "pos") {
  bu <- .unit_vector(beta)
  nb <- .filtered_neighbors(emb, if (side == "pos") bu else -bu, topn = topn, restrict = restrict_vocab)
  if (nrow(nb) < max(2L, k_min)) stop("Not enough neighbors to cluster.", call. = FALSE)
  w <- emb$vectors[nb$row, , drop = FALSE]
  km <- if (!is.null(k)) {
    .kmeans(w, min(as.integer(k), nrow(w)), random_state = random_state)
  } else {
    .kmeans_auto_k(w, k_min = k_min, k_max = min(k_max, nrow(w)), random_state = random_state)
  }
  clusters <- list()
  for (cid in sort(unique(km$labels))) {
    idx <- which(km$labels == cid)
    if (length(idx) < min_cluster_size) next
    wc <- w[idx, , drop = FALSE]
    centroid <- .unit_vector(colMeans(wc))
    cos_centroid <- drop(wc %*% centroid)
    o <- order(-cos_centroid, method = "radix")
    clusters[[length(clusters) + 1L]] <- list(
      id = cid,
      size = length(idx),
      centroid_cos_beta = sum(centroid * bu),
      coherence = mean(cos_centroid),
      words = data.frame(word = nb$word[idx][o], cos_centroid = cos_centroid[o],
                         cos_beta = drop(wc %*% bu)[o], stringsAsFactors = FALSE)
    )
  }
  key <- vapply(clusters, `[[`, 0, "centroid_cos_beta")
  clusters[order(if (side == "pos") -key else key, method = "radix")]
}
