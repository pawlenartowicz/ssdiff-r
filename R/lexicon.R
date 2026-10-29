# Port of ssdiff/utils/lexicon.py and the lexicon helpers of the SSD class in
# ssdiff/ssd.py (suggest_lexicon, token_stats, coverage_summary,
# evaluate_lexicon). Change them together.

#' Suggest lexicon words
#'
#' Ranks the tokens of a text as candidate lexicon words, or lists the tokens
#' close to seed words. Only the tokens are needed, so `s` may be a
#' tokens-only text (`ssd(texts, use_full_doc = FALSE)`).
#'
#' Without `seeds`, every token in at least `min_docs` documents is ranked by
#' balanced coverage (the mean coverage over quantile bins of `y`, or over
#' groups for categorical `y`) times `1 - min(1, |corr| / corr_cap)`, where
#' `corr` is the token's correlation with `y` (Cramér's V for categorical
#' `y`). Tokens that cover the documents evenly and do not track the outcome
#' rank highest. `pvalue` is the point-biserial test (continuous) or the
#' chi-square test (categorical, with Yates' correction for 2 groups).
#'
#' With `seeds`, each seed lists itself, its `top_k` closest candidate tokens
#' by cosine in the embeddings `s` was built on, and its forms with an English
#' ending (`-s`, `-es`, `-d`, `-ed`, `-ing`, and `-ies`, `-ied` after `y`). A
#' token appears once, under its closest seed, and rows are sorted by cosine.
#' A lexicon matches exact tokens, so a plural counts only when it is listed.
#'
#' @param s A text from [ssd()].
#' @param y The outcome: one value per document (numbers for continuous,
#'   labels for categorical), or a column name of `data`. Rows where it is
#'   missing are dropped. `NULL` is allowed with `seeds`: the outcome columns
#'   are then `NA`.
#' @param data A data frame with one row per document, holding the column
#'   `y` names.
#' @param seeds Character vector of seed words, each a key of the embeddings.
#'   A seed need not occur in the text.
#' @param top_k The number of tokens returned; with `seeds`, the number of
#'   neighbours listed per seed.
#' @param min_docs Minimum number of documents a token must be in.
#' @param n_bins Number of quantile bins of a continuous `y`.
#' @param corr_cap Correlation at which the rank penalty is complete.
#' @param var_type `"continuous"` or `"categorical"`.
#' @return An object of class `ssd_lexicon`, a list with `stats` (one-row data
#'   frame: `var_type`, `n_docs`, `n_tokens`) and `suggestions` (data frame:
#'   `token`, `freq`, `cov_all`, `cov_bal`, `corr`, `pvalue`, `direction`,
#'   `rank`; with seeds also `seed` first and `cos` after `token`).
#' @seealso [evaluate_lexicon()], [lexicon_tokens()]
#' @examples
#' \dontrun{
#' s <- ssd(df$essay, use_full_doc = FALSE)
#' lex <- suggest_lexicon(s, "score", data = df)
#' lex <- suggest_lexicon(s, seeds = c("work", "job"))
#' s2 <- with_lexicon(s, lexicon_tokens(lex))
#' }
#' @export
suggest_lexicon <- function(s, y = NULL, data = NULL, seeds = NULL, top_k = 30, min_docs = 5, n_bins = 4,
                            corr_cap = 0.30, var_type = "continuous") {
  .check_ssd(s)
  .check_var_type(var_type)
  if (!is.numeric(top_k) || length(top_k) != 1L || is.na(top_k) || top_k < 0) {
    stop("top_k must be >= 0, got ", format(top_k), call. = FALSE)
  }
  if (is.null(seeds)) {
    if (is.null(y)) stop("Ranking needs an outcome. To list neighbours of words without one, pass seeds.", call. = FALSE)
  } else {
    if (!is.character(seeds)) stop("seeds must be a character vector of words", call. = FALSE)
    if (!length(seeds)) stop("seeds is empty; pass at least one seed word.", call. = FALSE)
    emb <- .embeddings_from_ref(s$pack)
    missing <- seeds[wordpack_lookup(emb$pack, seeds) < 0L]
    if (length(missing)) {
      stop("These seeds are not in the embeddings: ", paste0("\"", missing, "\"", collapse = ", "), call. = FALSE)
    }
  }
  o <- .lexicon_outcome(s, y, data, "suggest_lexicon", var_type)
  pres <- .presence(o$docs)
  if (is.null(seeds)) {
    vocab <- names(pres)[lengths(pres) >= min_docs]
    sug <- .token_stats(pres, vocab, o$y, length(o$docs), n_bins, corr_cap, var_type)
    sug <- sug[seq_len(min(nrow(sug), top_k)), , drop = FALSE]
  } else {
    listed <- .seed_neighbours(pres, emb, seeds, top_k, min_docs)
    stats <- .token_stats(pres, listed$token, o$y, length(o$docs), n_bins, corr_cap, var_type)
    stats <- stats[match(listed$token, stats$token), names(stats) != "token", drop = FALSE]
    sug <- cbind(listed, stats)
  }
  rownames(sug) <- NULL
  .new_lexicon(var_type, length(o$docs), sug)
}

#' Evaluate a lexicon
#'
#' Per-token statistics of a lexicon (as in [suggest_lexicon()]) and its
#' coverage of the documents.
#'
#' @inheritParams suggest_lexicon
#' @param y The outcome, as in [suggest_lexicon()]; `NULL` gives only the
#'   statistics that need no outcome, the others are `NA`.
#' @param lexicon Character vector of lexicon words; a repeated word is
#'   evaluated once.
#' @return An object of class `ssd_lexicon` with `stats`, `suggestions` (one
#'   row per word, sorted by rank, or by `freq` without `y`) and `summary`, a
#'   one-row data frame: `docs_any` (documents with any lexicon word),
#'   `cov_all` (their share), `q1` and `q4` (that share in the lowest and
#'   highest bin of `y`, or the lowest and highest group), `corr_any` (its
#'   correlation with `y`, or Cramér's V), `hits_mean`, `hits_median`
#'   (occurrences per document), `types_mean`, `types_median` (distinct words
#'   per document); for categorical `y` also `group_cov`, a list column
#'   holding the share per group as a named numeric vector.
#' @examples
#' \dontrun{
#' ev <- evaluate_lexicon(s, "score", c("work", "job", "career"), data = df)
#' ev$summary
#' }
#' @export
evaluate_lexicon <- function(s, y, lexicon, data = NULL, n_bins = 4, corr_cap = 0.30,
                             var_type = "continuous") {
  .check_ssd(s)
  .check_var_type(var_type)
  if (!is.character(lexicon)) stop("lexicon must be a character vector of words", call. = FALSE)
  lexicon <- unique(lexicon)
  o <- .lexicon_outcome(s, y, data, "evaluate_lexicon", var_type)
  pres <- .presence(o$docs)
  sug <- .token_stats(pres, lexicon, o$y, length(o$docs), n_bins, corr_cap, var_type)
  summary <- .coverage_summary(o$docs, pres, lexicon, o$y, n_bins, var_type)
  .new_lexicon(var_type, length(o$docs), sug, summary)
}

#' Tokens of a lexicon suggestion
#'
#' @param lex A result of [suggest_lexicon()] or [evaluate_lexicon()].
#' @return The `token` column of `lex$suggestions`, ready to pass to [ssd()]
#'   or [with_lexicon()].
#' @export
lexicon_tokens <- function(lex) {
  if (!inherits(lex, "ssd_lexicon")) stop("lex must be a result of suggest_lexicon() or evaluate_lexicon()", call. = FALSE)
  lex$suggestions$token
}

#' Print a lexicon suggestion
#'
#' Shows the stats, the first 20 suggestions and, from [evaluate_lexicon()],
#' the coverage summary.
#'
#' @param x A result of [suggest_lexicon()] or [evaluate_lexicon()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.ssd_lexicon <- function(x, ...) {
  st <- x$stats
  cat(sprintf("<ssd_lexicon> %d tokens, %s documents, %s\n", st$n_tokens, format(st$n_docs, big.mark = ","),
              st$var_type))
  if (nrow(x$suggestions)) {
    shown <- utils::head(x$suggestions, 20L)
    if (nrow(x$suggestions) > 20L) cat("First 20 of", nrow(x$suggestions), "suggestions:\n")
    print(shown, row.names = FALSE, digits = 3)
  }
  if (!is.null(x$summary)) {
    cat("Coverage summary:\n")
    sm <- x$summary
    sm$group_cov <- NULL
    print(sm, row.names = FALSE, digits = 3)
    if (!is.null(x$summary$group_cov)) {
      cat("Coverage per group:\n")
      print(round(x$summary$group_cov[[1]], 3))
    }
  }
  invisible(x)
}

.check_var_type <- function(var_type) {
  if (!is.character(var_type) || length(var_type) != 1L || !var_type %in% c("continuous", "categorical")) {
    stop("var_type must be 'continuous' or 'categorical', got ", deparse1(var_type), call. = FALSE)
  }
}

.new_lexicon <- function(var_type, n_docs, suggestions, summary = NULL) {
  structure(list(
    stats = data.frame(var_type = var_type, n_docs = n_docs, n_tokens = nrow(suggestions)),
    suggestions = suggestions,
    summary = summary
  ), class = "ssd_lexicon")
}

# The documents and outcome a lexicon helper uses (_lexicon_outcome and
# _filter_y): `y = NULL` drops nothing; otherwise rows with a missing outcome
# are dropped. Categorical labels are kept as given.
.lexicon_outcome <- function(s, y, data, method, var_type) {
  if (is.null(y)) return(list(docs = s$docs, y = NULL))
  o <- .outcome_values(y, data, s$n_docs, var_type, method = method, check_data = FALSE)
  list(docs = s$docs[o$valid], y = o$values[o$valid])
}

# For each distinct token, the documents it is in, in order of first
# appearance.
.presence <- function(docs) {
  sets <- lapply(docs, unique)
  tokens <- unlist(sets, use.names = FALSE)
  keys <- unique(tokens)
  split(rep.int(seq_along(sets), lengths(sets)), factor(tokens, levels = keys))
}

# Per-token statistics (token_stats): one row per token with `freq`,
# `cov_all`, and the outcome columns (NA without `y`), sorted by rank, or by
# freq without `y`.
.token_stats <- function(pres, tokens, y, n, n_bins, corr_cap, var_type) {
  tokens <- as.character(tokens)
  docs_of <- unname(pres[match(tokens, names(pres))])
  freq <- lengths(docs_of)
  m <- length(tokens)
  out <- data.frame(token = tokens, freq = freq, cov_all = if (n) freq / n else rep(0, m),
                    cov_bal = rep(NA_real_, m), corr = rep(NA_real_, m), pvalue = rep(NA_real_, m),
                    direction = rep(NA_character_, m), rank = rep(NA_real_, m))
  if (is.null(y)) {
    out <- out[order(-out$freq, out$token, method = "radix"), , drop = FALSE]
  } else {
    categorical <- var_type == "categorical"
    group <- if (categorical) match(y, sort(unique(y), method = "radix")) else {
      bins <- .quantile_bins(y, n_bins)
      match(bins, sort(unique(bins)))
    }
    size <- tabulate(group)
    rows <- lapply(docs_of, function(d) {
      p <- numeric(n)
      p[d] <- 1
      if (categorical) .stats_categorical(p, group, size) else .stats_continuous(p, y, group, size)
    })
    out$cov_bal <- vapply(rows, `[[`, 0, "cov_bal")
    out$corr <- vapply(rows, `[[`, 0, "corr")
    out$pvalue <- vapply(rows, `[[`, 0, "pvalue")
    out$direction <- vapply(rows, `[[`, "", "direction")
    pen <- abs(out$corr) / corr_cap
    pen[is.na(pen) | pen > 1] <- 1   # Python's min(1.0, nan) is 1.0
    out$rank <- out$cov_bal * (1 - pen)
    out <- out[order(-out$rank, -out$cov_bal, -out$freq, out$token, method = "radix"), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

# Mean of p within each group 1..G (all groups present), in group order.
.group_means <- function(p, group, size) if (length(p)) as.vector(rowsum(p, group)) / size else numeric(0)

# _rank_for_token_stats, _pointbiserial_pvalue and _effect_direction for a
# continuous outcome; `group` numbers its quantile bins from 1.
.stats_continuous <- function(p, y, group, size) {
  covs <- .group_means(p, group, size)
  n <- length(p)
  constant <- length(unique(p)) < 2L
  r <- if (constant) NaN else suppressWarnings(stats::cor(p, y))
  pvalue <- if (constant || n < 3L || !is.finite(r)) NaN else if (abs(r) >= 1) 0 else {
    df <- n - 2
    2 * .t_sf(abs(r * sqrt(df / (1 - r * r))), df)
  }
  direction <- if (constant || !is.finite(r) || abs(r) < 1e-9) "none" else if (r > 0) "positive" else "negative"
  list(cov_bal = .np_mean(covs), corr = if (is.finite(r)) r else 0, pvalue = pvalue,
       direction = direction)
}

# The same for a categorical outcome; `group` indexes the sorted labels.
.stats_categorical <- function(p, group, size) {
  covs <- .group_means(p, group, size)
  direction <- if (length(covs) < 2L || max(covs) - min(covs) < 1e-9) "none" else
    if (covs[length(covs)] > covs[1]) "positive" else "negative"
  list(cov_bal = .np_mean(covs), corr = .cramers_v(p, group),
       pvalue = .chi2_pvalue(p, group), direction = direction)
}

# np.mean(): numpy's pairwise double sum over the length; 0 when empty. R's
# mean() sums in long double and corrects, which can differ in the last bit
# and reorder tied tokens.
.np_mean <- function(v) if (length(v)) .np_sum(v) / length(v) else 0

# numpy's pairwise sum: one by one below 8 values, eight running sums up to
# 128, halves (at a multiple of 8) above.
.np_sum <- function(a) {
  n <- length(a)
  if (n < 8L) {
    s <- 0
    for (x in a) s <- s + x
    return(s)
  }
  if (n > 128L) {
    half <- n %/% 2L
    half <- half - half %% 8L
    return(.np_sum(a[seq_len(half)]) + .np_sum(a[(half + 1L):n]))
  }
  r <- a[1:8]
  i <- 9L
  while (i + 7L <= n) {
    r <- r + a[i:(i + 7L)]
    i <- i + 8L
  }
  s <- ((r[1] + r[2]) + (r[3] + r[4])) + ((r[5] + r[6]) + (r[7] + r[8]))
  while (i <= n) {
    s <- s + a[i]
    i <- i + 1L
  }
  s
}

# Contingency table of presence p (0/1) by group 1..G (every group present):
# a row for each presence value that occurs, as table() gives. tabulate()
# instead of table(), which costs milliseconds per token.
.crosstab <- function(p, group) {
  hits <- tabulate(group[p > 0], max(group))
  ct <- rbind(tabulate(group, max(group)) - hits, hits)
  ct[rowSums(ct) > 0, , drop = FALSE]
}

.cramers_v <- function(p, group) {
  ct <- .crosstab(p, group)
  if (nrow(ct) < 2L || ncol(ct) < 2L) return(0)
  n <- sum(ct)
  e <- outer(rowSums(ct), colSums(ct)) / n
  nz <- e > 0
  chi2 <- sum((ct[nz] - e[nz])^2 / e[nz])
  k <- min(dim(ct)) - 1
  if (n * k > 0) sqrt(chi2 / (n * k)) else 0
}

# Chi-square test of independence, with Yates' correction for a 2 x 2 table.
.chi2_pvalue <- function(p, group) {
  ct <- .crosstab(p, group)
  if (nrow(ct) < 2L || ncol(ct) < 2L) return(NaN)
  n <- sum(ct)
  if (n == 0) return(NaN)
  e <- outer(rowSums(ct), colSums(ct)) / n
  nz <- e > 0
  d <- if (all(dim(ct) == 2L)) pmax(abs(ct - e) - 0.5, 0) else ct - e
  chi2 <- sum(d[nz]^2 / e[nz])
  df <- (nrow(ct) - 1) * (ncol(ct) - 1)
  if (df < 1 || chi2 < 0) return(NaN)
  if (chi2 <= 0) 1 else 1 - stats::pchisq(chi2, df)
}

# P(T > t) for Student's t, as Python's t_sf computes it.
.t_sf <- function(t, df) {
  if (t == 0) return(0.5)
  half <- 0.5 * stats::pbeta(df / (df + t * t), df / 2, 0.5)
  if (t > 0) half else 1 - half
}

# Bins 0..k-1 of y at its quantiles (_quantile_bins): numpy's percentile with
# linear interpolation at linspace(0, 100, n_bins + 1), duplicate edges
# dropped, then searchsorted on the inner edges. Falls back to a median split.
.quantile_bins <- function(y, n_bins = 4) {
  valid <- y[is.finite(y)]
  if (!length(valid)) return(integer(length(y)))
  q <- (0:n_bins) * (100 / n_bins)
  q[length(q)] <- 100
  edges <- unique(sort(.np_percentile(valid, q / 100)))
  if (length(edges) < 2L) return(as.integer(y > stats::median(y)))
  inner <- edges[-c(1L, length(edges))]
  if (!length(inner)) return(integer(length(y)))
  findInterval(y, inner)
}

# numpy.percentile(x, 100 * q) with method "linear", including numpy's lerp.
.np_percentile <- function(x, q) {
  x <- sort(x)
  n <- length(x)
  vi <- (n - 1) * q
  lo <- floor(vi)
  above <- vi >= n - 1
  a <- x[pmin(lo, n - 1) + 1]
  b <- x[pmin(lo + 1, n - 1) + 1]
  t <- vi - lo
  out <- a + (b - a) * t
  hi <- t >= 0.5
  out[hi] <- b[hi] - (b[hi] - a[hi]) * (1 - t[hi])
  out[above] <- x[n]
  out
}

# Coverage of a lexicon (coverage_summary): a one-row data frame, with a
# group_cov list column for categorical y.
.coverage_summary <- function(docs, pres, lexicon, y, n_bins, var_type) {
  categorical <- !is.null(y) && var_type == "categorical"
  if (!is.null(y) && !length(docs)) {
    out <- data.frame(docs_any = 0L, cov_all = 0, q1 = 0, q4 = 0, corr_any = 0, hits_mean = 0,
                      hits_median = 0, types_mean = 0, types_median = 0)
    if (categorical) out$group_cov <- list(stats::setNames(numeric(0), character(0)))
    return(out)
  }
  n <- length(docs)
  p <- numeric(n)
  p[unique(unlist(pres[intersect(lexicon, names(pres))], use.names = FALSE))] <- 1
  q1 <- q4 <- corr_any <- NA_real_
  group_cov <- NULL
  if (categorical) {
    labels <- sort(unique(y), method = "radix")
    group <- match(y, labels)
    group_cov <- stats::setNames(as.vector(rowsum(p, group)) / tabulate(group),
                                 if (is.numeric(labels)) .number_labels(labels) else labels)
    q1 <- min(group_cov)
    q4 <- max(group_cov)
    corr_any <- .cramers_v(p, group)
  } else if (!is.null(y)) {
    bins <- .quantile_bins(y, n_bins)
    q1 <- mean(p[bins == min(bins)])
    q4 <- mean(p[bins == max(bins)])
    r <- if (length(unique(p)) < 2L) 0 else suppressWarnings(stats::cor(p, y))
    corr_any <- if (is.finite(r)) r else 0
  }
  # Doubles, so the medians are doubles as numpy's (median() keeps integers).
  hits <- vapply(docs, function(d) as.double(sum(d %in% lexicon)), 0)
  types <- vapply(docs, function(d) as.double(length(intersect(d, lexicon))), 0)
  avg <- function(v, f) if (length(v)) f(v) else 0
  out <- data.frame(docs_any = as.integer(sum(p)), cov_all = avg(p, mean), q1 = q1, q4 = q4,
                    corr_any = corr_any, hits_mean = avg(hits, mean), hits_median = avg(hits, stats::median),
                    types_mean = avg(types, mean), types_median = avg(types, stats::median))
  if (categorical) out$group_cov <- list(group_cov)
  out
}

# Tokens to list for seeds (_seed_neighbours): a data frame with `seed`,
# `token` and `cos`. Candidates are tokens in at least `min_docs` documents
# that are keys of the pack. Each seed lists itself, its `top_k` closest
# other candidates and its candidate _seed_endings forms; a token keeps the
# row of the seed it is closest to (the first such seed on a tie), a seed
# its own row. Rows are sorted by cosine, then seed order, the seed's own row
# first; other ties keep the order in which tokens were first listed.
.seed_neighbours <- function(pres, emb, seeds, top_k, min_docs) {
  empty <- data.frame(seed = character(0), token = character(0), cos = numeric(0))
  frequent <- names(pres)[lengths(pres) >= min_docs]
  candidates <- sort(frequent[wordpack_lookup(emb$pack, frequent) >= 0L], method = "radix")
  if (!length(candidates)) return(empty)
  unit <- .l2_normalize_rows(.word_vectors(emb, candidates)[candidates, , drop = FALSE])
  seed_unit <- .l2_normalize_rows(.word_vectors(emb, seeds)[seeds, , drop = FALSE])

  tok <- character(0)
  best_seed <- character(0)
  best_cos <- numeric(0)
  for (s in seq_along(seeds)) {
    seed <- seeds[s]
    cos <- as.vector(unit %*% seed_unit[s, ])
    closest <- order(-cos, method = "radix")
    closest <- closest[candidates[closest] != seed]
    endings <- match(.seed_endings(seed), candidates)
    for (i in c(utils::head(closest, top_k), endings[!is.na(endings)])) {
      k <- match(candidates[i], tok)
      if (is.na(k)) {
        tok <- c(tok, candidates[i])
        best_seed <- c(best_seed, seed)
        best_cos <- c(best_cos, cos[i])
      } else if (cos[i] > best_cos[k]) {
        best_seed[k] <- seed
        best_cos[k] <- cos[i]
      }
    }
  }
  for (seed in intersect(seeds, candidates)) {
    k <- match(seed, tok)
    if (is.na(k)) {
      tok <- c(tok, seed)
      k <- length(tok)
    }
    best_seed[k] <- seed
    best_cos[k] <- 1
  }
  o <- order(-pmin(best_cos, 1), match(best_seed, seeds), tok != best_seed, method = "radix")
  data.frame(seed = best_seed[o], token = tok[o], cos = best_cos[o])
}

# Forms of a seed with an English ending, tried for every language.
.seed_endings <- function(seed) {
  forms <- paste0(seed, c("s", "es", "d", "ed", "ing"))
  stem <- substr(seed, 1L, nchar(seed) - 1L)
  if (endsWith(seed, "y")) forms <- c(forms, paste0(stem, c("ies", "ied")))
  if (endsWith(seed, "e")) forms <- c(forms, paste0(stem, "ing"))
  forms
}
