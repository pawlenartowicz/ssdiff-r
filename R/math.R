# Port of ssdiff/utils/math.py. Change the two together.
# Python's chi2_sf, f_sf and t_sf are not ported: R has pchisq, pf and pt.

# v / ||v||, or zeros when ||v|| < eps.
.unit_vector <- function(v, eps = 1e-12) {
  n <- sqrt(sum(v^2))
  if (n < eps) v * 0 else v / n
}

# Column z-scores with the population SD (ddof = 0); a column with SD <= 1e-12
# is scaled by 1. Returns `x`, `mean` and `scale`.
.standardize <- function(x) {
  x <- as.matrix(x)
  mean <- colMeans(x)
  xc <- sweep(x, 2L, mean)
  scale <- sqrt(colMeans(xc^2))
  scale[!(scale > 1e-12)] <- 1
  list(x = sweep(xc, 2L, scale, "/"), mean = mean, scale = scale)
}

# PCA by SVD of the centred matrix: `z` (scores, n x n_components),
# `components` (n_components x D) and `evr` (explained variance ratio).
# A component's sign may differ from numpy's.
.pca_fit_transform <- function(x, n_components) {
  x <- as.matrix(x)
  xc <- sweep(x, 2L, colMeans(x))
  s <- svd(xc, nu = 0L, nv = n_components)
  v <- s$v[, seq_len(n_components), drop = FALSE]
  var <- s$d^2 / (nrow(x) - 1)
  total <- sum(var)
  evr <- if (total > 0) var[seq_len(n_components)] / total else numeric(n_components)
  list(z = xc %*% v, components = t(v), evr = evr)
}

# Rows scaled to unit length; a row with norm below 1e-12 is divided by 1e-12.
.l2_normalize_rows <- function(x) x / pmax(sqrt(rowSums(x^2)), 1e-12)

# Squared distances between the rows of x and of centers, by Python's
# formula (|x|^2 + |c|^2 - 2 x.c, floored at 0), so the k-means steps round
# as Python's do.
.sq_dists <- function(x, centers, x_sq = rowSums(x^2)) {
  d <- outer(x_sq, rowSums(centers^2), "+") - 2 * tcrossprod(x, centers)
  d[d < 0] <- 0
  d
}

# Squared distance of every row of x to the vector v, summed per row.
.sq_dist_to <- function(x, v) rowSums((x - rep(v, each = nrow(x)))^2)

# k-means++ start with several candidates per step (05: "R call sites").
# Returns the row numbers of the k starting centres.
.kmeans_plus_plus <- function(x, k, r) {
  n <- nrow(x)
  n_local_trials <- 2L + as.integer(floor(log(k)))
  centers <- .rng_integers(r, n) + 1
  min_dists <- rep(Inf, n)
  for (i in seq_len(k - 1L)) {
    min_dists <- pmin(min_dists, .sq_dist_to(x, x[centers[length(centers)], ]))
    if (sum(min_dists) == 0) {
      idx <- .rng_integers(r, n) + 1
    } else {
      # Weights are not divided by their total: that would change the cumulative sums.
      candidates <- .rng_choice(r, min_dists, n_local_trials) + 1
      idx <- candidates[1L]
      best <- Inf
      for (cand in candidates) {
        potential <- sum(pmin(min_dists, .sq_dist_to(x, x[cand, ])))
        if (potential < best) {
          best <- potential
          idx <- cand
        }
      }
    }
    centers <- c(centers, idx)
  }
  centers
}

# Means of the rows of x per label 1..k (sequential sums, as numpy's mean over
# axis 0); a label with no rows gets a zero row. Returns `centers` and `counts`.
.label_means <- function(x, labels, k) {
  counts <- tabulate(labels, k)
  centers <- matrix(0, k, ncol(x))
  s <- rowsum(x, labels)
  present <- as.integer(rownames(s))
  centers[present, ] <- s / counts[present]
  list(centers = centers, counts = counts)
}

# Lloyd's k-means from a k-means++ start, the best of n_init runs by inertia.
# An empty cluster is re-seeded with the point farthest from its centre.
# Returns `labels` (1..k), `centers` (k x D) and `inertia`.
.kmeans <- function(x, k, random_state = 2137, max_iter = 300L, n_init = 1L) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (k > n) stop(sprintf("Cannot request k=%d clusters from n=%d samples.", k, n), call. = FALSE)
  r <- .rng(random_state)
  x_sq <- rowSums(x^2)
  best <- list(labels = integer(n), centers = matrix(0, k, ncol(x)), inertia = Inf)
  for (init in seq_len(n_init)) {
    centers <- x[.kmeans_plus_plus(x, k, r), , drop = FALSE]
    prev <- rep(-1L, n)
    for (it in seq_len(max_iter)) {
      d <- .sq_dists(x, centers, x_sq)
      labels <- max.col(-d, ties.method = "first")
      if (identical(labels, prev)) break
      prev <- labels
      upd <- .label_means(x, labels, k)
      empty <- which(upd$counts == 0L)
      if (length(empty)) {
        min_d <- d[cbind(seq_len(n), labels)]
        for (j in empty) {
          far <- which.max(min_d)
          upd$centers[j, ] <- x[far, ]
          min_d[far] <- -1
        }
      }
      centers <- upd$centers
    }
    d <- .sq_dists(x, centers, x_sq)
    labels <- max.col(-d, ties.method = "first")
    inertia <- sum(d[cbind(seq_len(n), labels)])
    if (inertia < best$inertia) best <- list(labels = labels, centers = centers, inertia = inertia)
  }
  best
}

# Euclidean distances between all rows of x (dot-product form, floored at 0).
.pairwise_euclidean <- function(x) {
  x_sq <- rowSums(x^2)
  d <- outer(x_sq, x_sq, "+") - 2 * tcrossprod(x)
  d[d < 0] <- 0
  sqrt(d)
}

# Mean silhouette coefficient from a distance matrix; 0 with fewer than 2 labels.
.silhouette <- function(dists, labels) {
  n <- length(labels)
  labs <- sort(unique(labels))
  if (length(labs) < 2L) return(0)
  a <- numeric(n)
  b <- rep(Inf, n)
  for (lab in labs) {
    inside <- labels == lab
    size <- sum(inside)
    if (size > 1L) a[inside] <- rowSums(dists[inside, inside, drop = FALSE]) / (size - 1)
    if (any(!inside)) b[!inside] <- pmin(b[!inside], rowMeans(dists[!inside, inside, drop = FALSE]))
  }
  denom <- pmax(a, b)
  ok <- is.finite(b) & denom > 0
  sil <- numeric(n)
  sil[ok] <- (b[ok] - a[ok]) / denom[ok]
  mean(sil)
}

# k-means for each k in k_min..k_max (capped at n - 1), keeping the k with the
# highest mean silhouette. Returns `.kmeans()`'s list plus `k`.
.kmeans_auto_k <- function(x, k_min = 2L, k_max = 5L, random_state = 2137, n_init = 1L) {
  x <- as.matrix(x)
  upper <- min(k_max, max(k_min, nrow(x) - 1L))
  dists <- .pairwise_euclidean(x)
  best_k <- max(2L, k_min)
  best_s <- -1
  best <- NULL
  # As Python's range(): with k_min > k_max no k is tried and k_min is fitted below.
  ks <- if (max(2L, k_min) <= max(2L, upper)) seq.int(max(2L, k_min), max(2L, upper)) else integer(0)
  for (k in ks) {
    km <- .kmeans(x, k, random_state = random_state, n_init = n_init)
    if (length(unique(km$labels)) <= 1L || max(tabulate(km$labels)) <= 1L) next
    s <- .silhouette(dists, km$labels)
    if (s > best_s) {
      best_s <- s
      best_k <- k
      best <- km
    }
  }
  if (is.null(best)) best <- .kmeans(x, best_k, random_state = random_state)
  c(best, list(k = as.integer(best_k)))
}
