# Port of ssdiff/backends/group.py and SSD.fit_groups in ssdiff/ssd.py (with
# the fields of GroupResult in ssdiff/results/group_result.py). Change them together.

.MIN_GROUP_SIZE <- 20L

# Python's correction names and the p.adjust() method each maps to.
.CORRECTIONS <- c(holm = "holm", bonferroni = "bonferroni", fdr_bh = "BH", none = "none")

#' Compare groups of documents
#'
#' Tests whether groups of documents differ in meaning: each group's
#' centroid is the unit-length mean of its document vectors, and two groups
#' differ by `T = 1 - cos` of their centroids. One permutation test gives the
#' omnibus p-value (for the mean `T` over all pairs) and a p-value for every
#' pair, each `(count(null >= observed) + 1) / (n_perm + 1)`. The pair
#' p-values are corrected for multiple comparisons.
#'
#' Groups with fewer than 20 documents are dropped, each with a warning. The
#' groups are named `g1`, `g2`, ... in sorted order of their labels
#' (`group_labels` holds the labels), and every pair `g1_g2`, ... is a unit
#' with its own gradient, the difference of the two groups' mean vectors.
#'
#' With `median_split = TRUE`, `y` is numeric and is split at its median into
#' `"low"` and `"high"`. Documents at the median are shared out at random
#' (with `random_state`) so that the two groups are as equal in size as
#' possible. No group is dropped then.
#'
#' The permutations come from ssdiff's own random generator, the one the
#' Python package uses, so the same data and `random_state` give the same
#' p-values in R and Python. R's `set.seed()` has no effect on them.
#'
#' @param s A text from [ssd()], with document vectors.
#' @param y The group labels: a factor, character or numeric vector with one
#'   value per document, or the name of a column of `data`. With
#'   `median_split = TRUE`, a numeric outcome. Documents with `NA` or `""`
#'   drop out.
#' @param data A data frame with one row per document, or `NULL`.
#' @param median_split `TRUE` splits a numeric `y` at its median.
#' @param n_perm Number of permutations.
#' @param correction Correction of the pair p-values: `"holm"`,
#'   `"bonferroni"`, `"fdr_bh"` (Benjamini-Hochberg) or `"none"`, computed
#'   with [stats::p.adjust()].
#' @param random_state Seed of the permutations and of the median split's
#'   tie-break, a whole number from 0 to 2^53 - 1.
#' @param verbose `TRUE` shows a progress bar and the omnibus p-value.
#' @return An object of class `c("ssd_groups", "ssd_result")`, a list with:
#'   \describe{
#'   \item{`G`, `n_kept`}{the number of groups and of documents used;}
#'   \item{`group_labels`}{a named character vector, `g1`, `g2`, ... to the
#'     labels;}
#'   \item{`groups`}{the group (`g1`, ...) of each document used;}
#'   \item{`stats`}{a one-row data frame: `G`, `n_kept`, `n_perm`,
#'     `correction`, `random_state`, `pvalue` (the omnibus p-value);}
#'   \item{`test`}{a one-row data frame: `name`, `pvalue`, `omnibus_T`,
#'     `omnibus_p`, `G`, `n_kept`, `n_perm`, `correction`, `random_state`;}
#'   \item{`pairs`}{a data frame with one row per pair: `contrast`, `g1`,
#'     `g2`, `T`, `p_raw`, `p_corrected`, `cohens_d` (on the documents'
#'     projections onto the pair's contrast), `n_g1`, `n_g2`, `contrast_norm`;}
#'   \item{`units`}{one entry per pair, named by `contrast`, with `g1`, `g2`,
#'     `beta` (mean of `g1` minus mean of `g2`), `gradient` (`beta` at unit
#'     length), `beta_norm` and `alignment_scores` (the cosine of each of the
#'     pair's documents with `gradient`);}
#'   \item{`texts`}{the text the fit used, with its tokens, settings, pack,
#'     and the vectors (`x`) and rows (`keep`) of the documents used;}
#'   \item{`fit_args`}{`median_split`, `n_perm`, `correction`, `random_state`.}
#'   }
#' @examples
#' \dontrun{
#' s <- ssd(df$essay, lexicon = c("work", "job"))
#' g <- fit_groups(s, "condition", data = df)
#' g$pairs
#' h <- fit_groups(s, df$score, median_split = TRUE)
#' }
#' @export
fit_groups <- function(s, y, data = NULL, median_split = FALSE, n_perm = 5000, correction = "holm",
                       random_state = 2137, verbose = FALSE) {
  fit_args <- mget(setdiff(names(formals()), c("s", "y", "data", "verbose")))
  .check_flag(median_split, "median_split")
  .check_flag(verbose, "verbose")
  if (!is.numeric(n_perm) || length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 0 || n_perm != round(n_perm)) {
    stop("n_perm must be a whole number >= 0, got ", format(n_perm), call. = FALSE)
  }
  .check_correction(correction)
  .rng(random_state)                                   # checks the seed

  rows <- .outcome_rows(s, y, data, if (median_split) "groups_median" else "groups")
  keep <- rows$keep
  x <- rows$x
  if (median_split) {
    groups <- .median_split(rows$y, random_state)
  } else {
    groups <- rows$y
    sizes <- table(groups)
    small <- sort(names(sizes)[sizes < .MIN_GROUP_SIZE], method = "radix")
    for (g in small) {
      warning(sprintf("Group '%s' dropped: only %d docs (minimum %d)", g, sizes[[g]], .MIN_GROUP_SIZE),
              call. = FALSE)
    }
    if (length(small)) {
      used <- !groups %in% small
      keep[keep] <- used
      x <- x[used, , drop = FALSE]
      groups <- groups[used]
    }
  }
  labels <- sort(unique(groups), method = "radix")
  G <- length(labels)
  if (G < 2L) stop(sprintf("Need at least 2 groups after filtering, got %d", G), call. = FALSE)

  canon <- paste0("g", seq_len(G))
  group_idx <- match(groups, labels)
  res <- .permutation_test(x, group_idx, G, n_perm, correction, random_state, verbose)

  means <- .label_means(x, group_idx, G)$centers
  units <- lapply(seq_len(nrow(res$pairs)), function(k) {
    a <- res$pairs$a[k]
    b <- res$pairs$b[k]
    beta <- means[a, ] - means[b, ]
    gradient <- .unit_vector(beta)
    in_pair <- group_idx == a | group_idx == b
    list(g1 = canon[a], g2 = canon[b], beta = beta, gradient = gradient, beta_norm = sqrt(sum(beta^2)),
         alignment_scores = .alignment(x[in_pair, , drop = FALSE], gradient))
  })
  contrast <- paste0(canon[res$pairs$a], "_", canon[res$pairs$b])
  names(units) <- contrast

  n_kept <- nrow(x)
  pairs <- data.frame(contrast = contrast, g1 = canon[res$pairs$a], g2 = canon[res$pairs$b],
                      res$pairs[c("T", "p_raw", "p_corrected", "cohens_d", "n_g1", "n_g2", "contrast_norm")],
                      stringsAsFactors = FALSE)
  structure(list(
    G = G,
    n_kept = n_kept,
    group_labels = stats::setNames(labels, canon),
    groups = canon[group_idx],
    stats = data.frame(G = G, n_kept = n_kept, n_perm = as.integer(n_perm), correction = correction,
                       random_state = random_state, pvalue = res$omnibus_p, stringsAsFactors = FALSE),
    test = data.frame(name = "permutation", pvalue = res$omnibus_p, omnibus_T = res$omnibus_T,
                      omnibus_p = res$omnibus_p, G = G, n_kept = n_kept, n_perm = as.integer(n_perm),
                      correction = correction, random_state = random_state, stringsAsFactors = FALSE),
    pairs = pairs,
    units = units,
    texts = list(.result_text(s, keep, x, rows$n_no_vector)),
    fit_args = fit_args
  ), class = c("ssd_groups", "ssd_result"))
}

# Python's correct_pvalues() by p.adjust(). Equal to Python's numbers, except
# that "fdr_bh" multiplies in another order (n / i * p, not p * n / i), so it
# can differ in the last bit.
.correct_pvalues <- function(p, correction) {
  .check_correction(correction)
  stats::p.adjust(p, method = .CORRECTIONS[[correction]])
}

.check_correction <- function(correction) {
  if (!is.character(correction) || length(correction) != 1L || !correction %in% names(.CORRECTIONS)) {
    stop(sprintf("Unknown correction %s. Choose from ('holm', 'bonferroni', 'fdr_bh', 'none')",
                 if (is.character(correction) && length(correction) == 1L) paste0("'", correction, "'")
                 else deparse1(correction)), call. = FALSE)
  }
}

# "low" / "high" at the median of y; values at the median are shared out in
# a random order (shared generator) so that "low" gets floor(n / 2) if it can.
.median_split <- function(y, random_state = 2137) {
  if (max(y) - min(y) == 0) {
    stop(sprintf("Cannot median-split: all y values are identical (%s)", .py_repr(y[1L])), call. = FALSE)
  }
  n <- length(y)
  s <- sort(y)
  # numpy's median: the middle value, or the mean of the two middle values as (a + b) / 2.
  med <- if (n %% 2L) s[(n + 1L) %/% 2L] else (s[n %/% 2L] + s[n %/% 2L + 1L]) / 2
  labels <- ifelse(y < med, "low", "high")
  tied <- which(y == med)
  if (length(tied)) {
    n_low <- max(0L, n %/% 2L - sum(y < med))
    tied <- tied[.rng_permutation(.rng(random_state), length(tied))]
    labels[tied] <- rep(c("low", "high"), c(min(n_low, length(tied)), max(0L, length(tied) - n_low)))
  }
  labels
}

# Unit-length group means of x for groups 1..G of group_idx.
.group_centroids <- function(x, group_idx, G) {
  .l2_normalize_rows(.label_means(x, group_idx, G)$centers)
}

# T = 1 - cos for the pairs (a[k], b[k]) of rows of the centroid matrix.
.pair_T <- function(centroids, a, b) {
  dot <- rowSums(centroids[a, , drop = FALSE] * centroids[b, , drop = FALSE])
  1 - pmin(pmax(dot, -1), 1)
}

# Mean of each row as numpy's mean: its pairwise sum (.np_sum), which is
# left to right only below 8 values (fewer than 5 groups).
.row_means_seq <- function(m) apply(m, 1L, .np_sum) / ncol(m)

# The omnibus and pairwise permutation test (unified_permutation_test):
# `omnibus_T`, `omnibus_p`, and `pairs`, a data frame with the group numbers
# `a`, `b` and T, p_raw, p_corrected, cohens_d, n_g1, n_g2, contrast_norm.
.permutation_test <- function(x, group_idx, G, n_perm, correction, random_state, verbose) {
  .check_correction(correction)
  pr <- utils::combn(G, 2L)
  a <- pr[1L, ]
  b <- pr[2L, ]
  n <- length(group_idx)
  centroids <- .group_centroids(x, group_idx, G)
  t_obs <- .pair_T(centroids, a, b)
  omnibus_t <- .row_means_seq(matrix(t_obs, 1L))

  r <- .rng(random_state)
  null_t <- matrix(0, n_perm, length(a))
  bar <- if (verbose && n_perm > 0) utils::txtProgressBar(min = 0, max = n_perm, style = 3) else NULL
  for (p in seq_len(n_perm)) {
    perm <- group_idx[.rng_permutation(r, n)]
    null_t[p, ] <- .pair_T(.group_centroids(x, perm, G), a, b)
    if (!is.null(bar)) utils::setTxtProgressBar(bar, p)
  }
  if (!is.null(bar)) close(bar)
  null_omnibus <- .row_means_seq(null_t)

  omnibus_p <- (sum(null_omnibus >= omnibus_t) + 1) / (n_perm + 1)
  if (verbose) message(sprintf("[groups] omnibus p=%s (%d groups, %d perms)", .fmt_p(omnibus_p), G, n_perm))
  p_raw <- (colSums(null_t >= rep(t_obs, each = n_perm)) + 1) / (n_perm + 1)
  p_corrected <- .correct_pvalues(p_raw, correction)

  stat <- lapply(seq_along(a), function(k) {
    contrast <- centroids[a[k], ] - centroids[b[k], ]
    proj <- drop(x %*% .unit_vector(contrast))
    p1 <- proj[group_idx == a[k]]
    p2 <- proj[group_idx == b[k]]
    dof <- length(p1) + length(p2) - 2L
    pooled <- if (dof > 0) {
      sqrt(((length(p1) - 1) * .var1(p1) + (length(p2) - 1) * .var1(p2)) / dof)
    } else 0
    c(cohens_d = (mean(p1) - mean(p2)) / max(pooled, 1e-12), n_g1 = length(p1), n_g2 = length(p2),
      contrast_norm = sqrt(sum(contrast^2)))
  })
  stat <- do.call(rbind, stat)
  list(omnibus_T = omnibus_t, omnibus_p = omnibus_p,
       pairs = data.frame(a = a, b = b, T = t_obs, p_raw = p_raw, p_corrected = p_corrected,
                          cohens_d = stat[, "cohens_d"], n_g1 = as.integer(stat[, "n_g1"]),
                          n_g2 = as.integer(stat[, "n_g2"]), contrast_norm = stat[, "contrast_norm"],
                          row.names = NULL))
}

# Sample variance (n - 1); NaN for one value, as numpy's var(ddof = 1).
.var1 <- function(v) if (length(v) > 1L) stats::var(v) else NaN

# A p-value as Python's fmt_p writes it: "<.001", else 3 decimals, no leading zero.
.fmt_p <- function(p) {
  if (!is.finite(p)) return("nan")
  if (p < 0.001) return("<.001")
  sub("^0", "", sprintf("%.3f", p))
}
