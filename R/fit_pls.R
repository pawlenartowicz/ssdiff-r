# Port of SSD.fit_pls, SSD.fit_multipls, SSD._orient_beta and
# _regression_gradients in ssdiff/ssd.py, and of ssdiff/backends/pls.py (with
# the fields of PLSResult and MultiPLSResult in ssdiff/results). Change them
# together.

#' Fit a PLS gradient
#'
#' Fits a partial least squares regression (PLS1) of a numeric outcome on the
#' text's document vectors, both standardized, and turns its coefficients
#' into a semantic gradient (`beta`): the direction in embedding space along
#' which the outcome grows. `beta` is oriented so that a higher alignment
#' means a higher outcome.
#'
#' With `k = "auto"`, plskit's `pls1_find_k_optimal()` (selector `"r2_se"`,
#' five-fold cross-validation) picks the number of components, and the
#' p-value is plskit's confirmatory test at k = 1, independent of that
#' choice. With a whole number `k`, the model has `k` components and the test
#' runs at `k`. The test is plskit's `test_method = "auto"` (Lenartowicz and
#' Plisiecki, 2026), which runs `split_exact` or `split_nb`, chosen from the
#' vectors alone; `test$name` says which ran.
#'
#' @section Covariates:
#' With a formula such as `score ~ age + gender`, the outcome and every
#' column of the document vectors are first regressed on the formula's fixed
#' effects by least squares, and PLS runs on the residuals. `r2` is then the
#' residualized outcome's, and the residuals the fit reports have the
#' covariates taken out. The covariates' own coefficients are not reported:
#' use [fit_regression()] for them, and for random effects, which `fit_pls()`
#' does not take.
#'
#' @section Several gradients:
#' With `ssd_gradients = N > 1`, the result also holds `N` gradients,
#' `units[["dim-1"]]` ... `units[["dim-N"]]`, from the weight vectors of an
#' N-component PLS fit on the same data. With `rotate = "words"` they are
#' varimax-rotated towards the projection of the 50,000 most frequent words
#' of the pack, ordered by the correlation of their scores with the outcome,
#' and signed so that it is positive. `rotate = NULL` keeps them in PLS
#' order. The standard fit does not change.
#'
#' @param s A text from [ssd()], with document vectors.
#' @param y The numeric outcome: values, one per document, a column name of
#'   `data`, or a formula whose left side names the outcome column and whose
#'   right side holds fixed covariates. Documents with a missing outcome drop
#'   out.
#' @param data A data frame with one row per document, or `NULL`.
#' @param k The number of PLS components, or `"auto"`.
#' @param k_max The largest k that `k = "auto"` tries, capped at `n - 1` and
#'   the vector dimension.
#' @param ssd_gradients The number of gradients, from 1 to `min(n - 1, D)`.
#' @param rotate `"words"` or `NULL`: the rotation of the gradients when
#'   `ssd_gradients > 1`.
#' @param n_splits Random splits of the split-half tests.
#' @param random_state Seed of plskit's random splits, a whole number from 0
#'   to 2^53 - 1.
#' @param test `FALSE` skips the confirmatory test; `pvalue` is then `NA`.
#' @param verbose `TRUE` shows the steps.
#' @return An object of class `c("ssd_pls", "ssd_result")`, a list with:
#'   \describe{
#'   \item{`stats`}{a one-row data frame: `backend` (`"PLS"`), `r2`,
#'     `pvalue`, `n_raw`, `n_kept`, `n_dropped`, `y_mean`, `y_std`,
#'     `beta_norm`, `delta`, `iqr_effect`, `y_corr_pred`;}
#'   \item{`fit_info`}{a one-row data frame of the fit's settings, the number
#'     of components (`n_components`) and the k of the test (`p_at_k`)
#'     included;}
#'   \item{`test`}{a one-row data frame: `name` (the method plskit ran),
#'     `pvalue`, `split_r2`, `n_splits`, `random_state`;}
#'   \item{`beta`, `gradient`}{the gradient, on the standardized-outcome
#'     scale, and its unit vector;}
#'   \item{`y`, `alignment_scores`, `y_hat`, `residuals`}{the outcome of the
#'     rows used, their cosine with the gradient, the predicted outcome and
#'     the residual (with the covariates taken out, given a formula);}
#'   \item{`outcome_name`}{the outcome's name;}
#'   \item{`component_scores`, `component_weights`}{the PLS scores T and
#'     weights W;}
#'   \item{`find_k`}{plskit's k-selection result for `k = "auto"`, else
#'     `NULL`;}
#'   \item{`pls_info`, `units`}{with `ssd_gradients > 1`: the rotation's
#'     diagnostics (`ssd_gradients`, `rotate`, `order`, `signs`,
#'     `kaiser_normalized`, `sweeps`, `V_converged`, `random_state`), and one
#'     unit per gradient, keyed `"dim-1"`, ..., each with `beta`, `gradient`
#'     and `stats` (`cov_share`, `score_corr`); else `NULL`;}
#'   \item{`texts`}{the text the fit used, with its tokens, settings, pack,
#'     and the vectors (`x`) and rows (`keep`) of the documents used;}
#'   \item{`fit_args`}{the formula (as a string) and the other arguments.}
#'   }
#' @examples
#' \dontrun{
#' s <- ssd(df$essay, lexicon = c("work", "job"))
#' r <- fit_pls(s, "score", data = df)
#' r <- fit_pls(s, score ~ age + gender, data = df, k = 2)
#' r <- fit_pls(s, df$score, ssd_gradients = 3)
#' r$units[["dim-1"]]$stats
#' }
#' @export
fit_pls <- function(s, y, data = NULL, k = "auto", k_max = 5, ssd_gradients = 1, rotate = "words", n_splits = 50,
                    random_state = 2137, test = TRUE, verbose = FALSE) {
  need_pkg("plskit", "fit_pls()")
  if (!identical(k, "auto") && !.is_whole(k)) stop(sprintf("k must be a whole number or \"auto\", got %s", deparse1(k)),
                                                    call. = FALSE)
  if (!identical(k, "auto") && k < 1) stop(sprintf("k must be >= 1, got %s", deparse1(k)), call. = FALSE)
  .check_count(k_max, "k_max")
  .check_count(n_splits, "n_splits")
  .check_flag(test, "test")
  .check_flag(verbose, "verbose")
  .rng(random_state)                                   # checks the seed
  if (inherits(y, "formula") && length(.text_names(y))) {
    stop("fit_pls() takes covariates only in its formula, no ssd() term; the text is the SSD itself.", call. = FALSE)
  }
  rows <- .outcome_rows(s, y, data, "pls")
  formula <- rows$formula
  fit_args <- c(list(formula = if (is.null(formula)) NULL else deparse1(formula)),
                mget(c("k", "k_max", "ssd_gradients", "rotate", "n_splits", "random_state", "test")))
  n_kept <- rows$n_kept
  D <- ncol(rows$x)
  g_max <- min(n_kept - 1L, D)
  if (!.is_whole(ssd_gradients) || ssd_gradients < 1 || ssd_gradients > g_max) {
    stop(sprintf("ssd_gradients must be a whole number in [1, %d] (min(n - 1, D) = min(%d, %d)); got %s",
                 g_max, n_kept - 1L, D, deparse1(ssd_gradients)), call. = FALSE)
  }
  if (ssd_gradients > 1 && !(is.null(rotate) || identical(rotate, "words"))) {
    stop(sprintf("fit_pls() rotate must be 'words' or NULL; got %s. One outcome gives no outcome loadings to rotate against.",
                 .shown(rotate)), call. = FALSE)
  }
  y_std <- .sd_pop(rows$y)
  if (y_std <= 1e-12) {
    stop(sprintf("fit_pls() requires y with nonzero variance (std=%s). PLS1 is undefined when y is constant.",
                 sprintf("%.3g", y_std)), call. = FALSE)
  }

  x_fit <- rows$x
  y_fit <- rows$y
  if (!is.null(formula)) {
    res <- .pls_residualize(formula, rows$data, rows$y, rows$x)
    x_fit <- res$x
    y_fit <- res$y
    if (.sd_pop(y_fit) <= 1e-8 * y_std) {
      stop("The formula's covariates explain y fully, so no variance is left for the text.", call. = FALSE)
    }
  }
  y_st <- .standardize(y_fit)
  ys <- drop(y_st$x)
  y_mean <- y_st$mean[[1L]]
  y_scale <- y_st$scale[[1L]]
  x_st <- .standardize(x_fit)
  xs <- x_st$x
  x_mean <- x_st$mean
  if (!is.null(formula)) {
    # Residuals have mean 0: y_hat centres on y's own mean and the raw
    # vectors' mean, and y_scale stays the residuals' SD, the unit beta is in.
    y_mean <- mean(rows$y)
    x_mean <- colMeans(rows$x)
  }

  find_k <- NULL
  test_name <- NA_character_
  pvalue <- NA_real_
  confirm <- function(at) {
    r <- plskit::pls1_confirmatory_test(xs, ys, at, test_method = "auto", args = list(n_splits = as.integer(n_splits)),
                                        pre_standardized = TRUE, seed = random_state, verbose = verbose)
    test_name <<- r$test_method
    pvalue <<- as.numeric(r$pvalue)
  }
  if (identical(k, "auto")) {
    eff_k_max <- max(min(as.integer(k_max), n_kept - 1L, D), 1L)
    if (verbose) message(sprintf("[fit_pls] auto-selecting k via find_k_optimal (selector=r2_se, k_max=%d)", eff_k_max))
    find_k <- plskit::pls1_find_k_optimal(xs, ys, eff_k_max, selector = "r2_se", diagnostic = "split_nb",
                                          args = list(n_folds = 5L, n_splits = as.integer(n_splits)),
                                          pre_standardized = TRUE, seed = random_state, verbose = verbose)
    k_star <- as.integer(find_k$k_star)
    if (k_star < 1L) {
      stop("fit_pls(k = \"auto\"): plskit found no component (k_star = 0); y shows no signal in X.", call. = FALSE)
    }
    if (verbose) {
      message(sprintf("[fit_pls] k_star = %d", k_star))
      if (k_star > 1L) {
        p <- find_k$pvalues[seq_len(k_star)]
        message(sprintf("[fit_pls] diagnostic chain (split_nb, k=1..%d): %s", k_star,
                        paste(sprintf("%.3g", p), collapse = ", ")))
        message(sprintf("        path-max: %.3g  [same-sample]", max(p, na.rm = TRUE)))
      }
    }
    if (test) {
      if (verbose) {
        message(sprintf("[fit_pls] running k=1 confirmatory test (test_method='auto', n_splits=%d)", as.integer(n_splits)))
      }
      confirm(1L)
      if (verbose) message(sprintf("[fit_pls] p-value (k=1, %s) = %.6g", test_name, pvalue))
    }
    n_comp <- k_star
    p_at_k <- 1L
  } else {
    n_comp <- as.integer(min(k, n_kept - 1L, D))
    if (test) confirm(n_comp)
    p_at_k <- n_comp
  }

  m <- plskit::pls1_fit(xs, ys, k = n_comp, pre_standardized = TRUE, seed = random_state)
  y_pred <- drop(xs %*% m$coef)
  ss_tot <- sum((ys - mean(ys))^2)
  r2 <- if (ss_tot > 0) 1 - sum((ys - y_pred)^2) / ss_tot else 0
  beta <- .orient_beta(m$coef / x_st$scale, ys, x_fit)

  grad <- NULL
  if (ssd_gradients > 1) {
    # NIPALS weights do not depend on how many components are fitted, so W's
    # leading columns are reused when the main fit has enough.
    W <- if (ssd_gradients <= ncol(m$W)) m$W[, seq_len(ssd_gradients), drop = FALSE] else
      plskit::pls1_fit(xs, ys, k = as.integer(ssd_gradients), pre_standardized = TRUE, seed = random_state)$W
    if (ncol(W) < ssd_gradients) {
      stop(sprintf(paste0("plskit::pls1_fit returned %d weight vectors but ssd_gradients=%d was requested ",
                          "(NIPALS deflation collapsed). Reduce ssd_gradients."), ncol(W), ssd_gradients), call. = FALSE)
    }
    L <- if (identical(rotate, "words")) .word_target(.embeddings_from_ref(s$pack), W, x_mean, x_st$scale)
    grad <- .pls_gradients(W, xs, ys, x_st$scale, rotate, L, random_state)
  }

  y_hat <- y_mean + y_scale * drop(sweep(rows$x, 2L, x_mean) %*% beta)
  unit <- .regression_unit(beta, rows$x, rows$y, drop(rows$x %*% beta), y_scale)
  tested <- !is.na(test_name)
  structure(list(
    stats = data.frame(backend = "PLS", r2 = r2, pvalue = pvalue, n_raw = rows$n_raw, n_kept = n_kept,
                       n_dropped = rows$n_dropped, y_mean = mean(rows$y), y_std = y_std, unit$stats,
                       stringsAsFactors = FALSE),
    fit_info = data.frame(
      n_components = n_comp, p_at_k = p_at_k, n_splits = as.integer(n_splits), random_state = random_state,
      k_min = NA_integer_, k_max = NA_integer_, k_step = NA_integer_, best_k = NA_integer_,
      pca_k_source = NA_character_, formula = if (is.null(formula)) NA_character_ else deparse1(formula),
      family = NA_character_, ssd_gradients = as.integer(ssd_gradients),
      rotate = if (ssd_gradients > 1 && !is.null(rotate)) rotate else NA_character_, outcome = rows$outcome_name,
      stringsAsFactors = FALSE),
    test = data.frame(name = test_name, pvalue = pvalue, split_r2 = NA_real_,
                      n_splits = if (tested) as.integer(n_splits) else NA_integer_,
                      random_state = if (tested) random_state else NA_real_, stringsAsFactors = FALSE),
    beta = unit$beta,
    gradient = unit$gradient,
    y = rows$y,
    alignment_scores = unit$alignment_scores,
    y_hat = y_hat,
    residuals = if (is.null(formula)) rows$y - y_hat else (ys - y_pred) * y_scale,
    outcome_name = rows$outcome_name,
    component_scores = m$T,
    component_weights = m$W,
    find_k = find_k,
    pls_info = grad$pls_info,
    units = grad$units,
    texts = list(.result_text(s, rows$keep, rows$x, rows$n_no_vector)),
    fit_args = fit_args
  ), class = c("ssd_pls", "ssd_result"))
}

#' Fit PLS gradients for several outcomes
#'
#' Fits PLS3 (PLS correlation): the singular value decomposition of the
#' cross-covariance of the standardized document vectors and the
#' standardized outcomes. Dimension i pairs a direction in embedding space
#' (`U[, i]`) with outcome loadings (`V[, i]`). Each dimension is a unit,
#' `units[["dim-1"]]`, ..., with its own gradient.
#'
#' Significance is plskit's PLS3 test of the first dimension, with
#' `test_method = "auto"` (`split_exact` or `split_nb`, chosen from the
#' vectors alone): one p-value for the fit, the same for every `rotate`.
#' Dimensions have no p-value of their own.
#'
#' @section Rotation:
#' One orthogonal `R` is applied to both sides (`U_rot = U %*% R`, `V_rot = V
#' %*% R`):
#' - `NULL`: the dimensions of the decomposition as they are.
#' - `"words"`: varimax towards the projection of the 50,000 most frequent
#'   words of the pack onto `U`.
#' - `"outcomes"`: varimax towards the outcome loadings `V`, or towards one
#'   row per group of `outcome_groups`. With `ssd_gradients` equal to the
#'   number of outcomes and no groups, `R = t(V)`: each dimension loads 1 on
#'   one outcome.
#'
#' The rotated modes order the dimensions by their squared cross-covariance.
#' Every mode signs each dimension so that its largest outcome loading is
#' positive.
#'
#' @param s A text from [ssd()], with document vectors.
#' @param y The outcomes, numeric: a matrix or data frame with one row per
#'   document and at least two columns, or column names of `data`. A row with
#'   a missing value in any outcome drops the document.
#' @param data A data frame with one row per document, or `NULL`.
#' @param ssd_gradients The number of dimensions, from 1 to `min(q, n - 1,
#'   D)` for `q` outcomes. `NULL` means `q`, or with `outcome_groups` the
#'   number of rotation target rows.
#' @param rotate `NULL`, `"words"` or `"outcomes"`.
#' @param outcome_groups With `rotate = "outcomes"` only: a named list of
#'   character vectors of outcome names. Each group becomes one rotation
#'   target row, the mean of its outcomes' loadings; outcomes in no group keep
#'   their own row. `ssd_gradients` must not exceed the number of rows.
#' @param n_splits Random splits of the test.
#' @param n_perm Permutations of the test; unused when `split_nb` runs.
#' @param random_state Seed of the test, a whole number from 0 to 2^53 - 1.
#' @param verbose `TRUE` shows plskit's progress.
#' @return An object of class `c("ssd_multipls", "ssd_result")`, a list
#'   with:
#'   \describe{
#'   \item{`stats`}{a one-row data frame: `n`, `q`, `ssd_gradients`,
#'     `rotate`, `cov_explained` (the share of the squared cross-covariance
#'     the dimensions hold), `pvalue`;}
#'   \item{`pls_info`}{the rotation's diagnostics: `ssd_gradients`, `rotate`,
#'     `order`, `signs`, `kaiser_normalized`, `sweeps`, `V_converged`,
#'     `random_state`;}
#'   \item{`test`}{a one-row data frame: `name` (the method plskit ran),
#'     `pvalue`, `statistic`, `n_splits`, `n_perm`, `random_state`;}
#'   \item{`y`, `outcome_names`, `outcome_groups`}{the outcomes of the rows
#'     used, their names, and the groups;}
#'   \item{`U`, `V`, `singular_values`, `U_rot`, `V_rot`, `R`, `x_scores`,
#'     `y_scores`}{the decomposition, its rotation, and the scores of the
#'     standardized vectors and outcomes on the rotated dimensions;}
#'   \item{`units`}{one per dimension, keyed `"dim-1"`, ..., each with
#'     `beta` (`U_rot[, i]` in the scale of the document vectors), `gradient`
#'     and `stats` (`cov_share`, `score_corr`);}
#'   \item{`texts`}{the text the fit used, with its tokens, settings, pack,
#'     and the vectors (`x`) and rows (`keep`) of the documents used;}
#'   \item{`fit_args`}{the fit's arguments.}
#'   }
#' @examples
#' \dontrun{
#' s <- ssd(df$essay, lexicon = c("work", "job"))
#' r <- fit_multipls(s, c("valence", "arousal"), data = df)
#' r <- fit_multipls(s, c("v1", "v2", "a1"), data = df, rotate = "outcomes",
#'                   outcome_groups = list(valence = c("v1", "v2")))
#' r$units[["dim-1"]]$stats
#' }
#' @export
fit_multipls <- function(s, y, data = NULL, ssd_gradients = NULL, rotate = NULL, outcome_groups = NULL, n_splits = 50,
                         n_perm = 1000, random_state = 2137, verbose = FALSE) {
  need_pkg("plskit", "fit_multipls()")
  fit_args <- mget(setdiff(names(formals()), c("s", "y", "data", "verbose")))
  .check_count(n_splits, "n_splits")
  .check_count(n_perm, "n_perm")
  .check_flag(verbose, "verbose")
  .rng(random_state)                                   # checks the seed
  rows <- .outcome_rows(s, y, data, "multipls")
  out_names <- rows$outcome_names
  n_kept <- rows$n_kept
  D <- ncol(rows$x)
  q <- ncol(rows$y)

  constant <- out_names[apply(rows$y, 2L, .sd_pop) <= 1e-12]
  if (length(constant)) {
    stop(sprintf("fit_multipls() requires outcomes with nonzero variance; constant: %s.", .quoted(constant)),
         call. = FALSE)
  }
  if (!(is.null(rotate) || (is.character(rotate) && length(rotate) == 1L && rotate %in% c("words", "outcomes")))) {
    stop(sprintf("rotate must be NULL, 'words' or 'outcomes'; got %s", .shown(rotate)), call. = FALSE)
  }

  grouped <- character(0)                              # outcome -> its group
  if (!is.null(outcome_groups)) {
    if (!identical(rotate, "outcomes")) {
      stop(sprintf("outcome_groups requires rotate='outcomes'; got rotate=%s", .shown(rotate)), call. = FALSE)
    }
    if (!is.list(outcome_groups) || (length(outcome_groups) && (is.null(names(outcome_groups)) ||
        any(!nzchar(names(outcome_groups))) || anyDuplicated(names(outcome_groups)))) ||
        !all(vapply(outcome_groups, function(m) is.null(m) || is.character(m), TRUE))) {
      stop("outcome_groups must be a list of character vectors with unique names", call. = FALSE)
    }
    for (g in names(outcome_groups)) {
      members <- outcome_groups[[g]]
      if (!length(members)) stop(sprintf("outcome_groups[['%s']] is empty.", g), call. = FALSE)
      for (name in members) {
        if (!name %in% out_names) {
          stop(sprintf("outcome_groups[['%s']]: unknown outcome '%s'; outcomes are %s.", g, name, .quoted(out_names)),
               call. = FALSE)
        }
        if (name %in% names(grouped)) {
          stop(sprintf("outcome '%s' is in two groups: '%s' and '%s'.", name, grouped[[name]], g), call. = FALSE)
        }
        grouped[[name]] <- g
      }
    }
  }

  # One rotation target row per group plus one per ungrouped outcome.
  n_rows <- if (length(outcome_groups)) length(outcome_groups) + q - length(grouped) else q
  default_note <- if (is.null(ssd_gradients)) " (the default: one per outcome, or per group with outcome_groups)" else ""
  if (is.null(ssd_gradients)) ssd_gradients <- n_rows
  k_max <- min(q, n_kept - 1L, D)
  if (!.is_whole(ssd_gradients) || ssd_gradients < 1 || ssd_gradients > k_max) {
    stop(sprintf("ssd_gradients must be a whole number in [1, %d] (min(q, n - 1, D) = min(%d, %d, %d)); got %s%s",
                 k_max, q, n_kept - 1L, D, deparse1(ssd_gradients), default_note), call. = FALSE)
  }
  k <- as.integer(ssd_gradients)
  if (k > n_rows) {
    stop(sprintf(paste0("rotate='outcomes' needs ssd_gradients <= the number of rotation target rows (%d: one per ",
                        "group plus one per ungrouped outcome); got ssd_gradients=%d"), n_rows, k), call. = FALSE)
  }

  x_st <- .standardize(rows$x)
  xs <- x_st$x
  ys <- unname(.standardize(rows$y)$x)
  m <- tryCatch(
    plskit::pls3_fit(xs, ys, k = k, pre_standardized_X = TRUE, pre_standardized_Y = TRUE),
    plskit_error = function(e) stop(sprintf("plskit::pls3_fit rejected ssd_gradients=%d: %s", k, conditionMessage(e)),
                                    call. = FALSE))
  if (m$k_used < k) {
    stop(sprintf(paste0("plskit::pls3_fit kept k_used=%d components of ssd_gradients=%d%s; the rest fall below its ",
                        "truncation floor (e.g. collinear outcomes). Refit with ssd_gradients <= %d."),
                 m$k_used, k, default_note, m$k_used), call. = FALSE)
  }
  U <- m$U
  V <- m$V
  sigma <- m$singular_values

  # rotate = "outcomes" without groups leaves L = NULL: .pls3_rotate() targets V itself.
  L <- NULL
  if (identical(rotate, "words")) {
    L <- .word_target(.embeddings_from_ref(s$pack), U, x_st$mean, x_st$scale)
  } else if (identical(rotate, "outcomes") && !is.null(outcome_groups)) {
    L <- rbind(do.call(rbind, lapply(outcome_groups, function(m) colMeans(V[match(m, out_names), , drop = FALSE]))),
               V[!out_names %in% names(grouped), , drop = FALSE])
  }
  rot <- .pls3_rotate(U, V, sigma, rotate, L)
  cov_explained <- sum(sigma^2) / sum(crossprod(xs, ys)^2)

  if (verbose) {
    message(sprintf("[fit_multipls] running LV1 test (test_method='auto', n_splits=%d, n_perm=%d)",
                    as.integer(n_splits), as.integer(n_perm)))
  }
  r <- plskit::pls3_confirmatory_test(xs, ys, 1L, test_method = "auto",
                                      args = list(n_perm = as.integer(n_perm), n_splits = as.integer(n_splits)),
                                      seed = random_state, verbose = verbose)

  x_scores <- xs %*% rot$U_rot
  y_scores <- ys %*% rot$V_rot
  sigma2 <- sigma^2
  units <- lapply(seq_len(k), function(i) {
    beta <- rot$U_rot[, i] / x_st$scale
    list(beta = beta, gradient = .unit_vector(beta),
         stats = data.frame(cov_share = sum(rot$R[, i]^2 * sigma2) / sum(sigma2),
                            score_corr = stats::cor(x_scores[, i], y_scores[, i])))
  })
  names(units) <- paste0("dim-", seq_len(k))
  structure(list(
    stats = data.frame(n = n_kept, q = q, ssd_gradients = k, rotate = if (is.null(rotate)) NA_character_ else rotate,
                       cov_explained = cov_explained, pvalue = as.numeric(r$pvalue), stringsAsFactors = FALSE),
    pls_info = .pls_info(k, rotate, rot, random_state),
    test = data.frame(name = r$test_method, pvalue = as.numeric(r$pvalue), statistic = as.numeric(r$statistic),
                      n_splits = as.integer(n_splits), n_perm = if (is.null(r$n_perm)) NA_integer_ else as.integer(r$n_perm),
                      random_state = random_state, stringsAsFactors = FALSE),
    y = rows$y,
    outcome_names = out_names,
    outcome_groups = outcome_groups,
    U = U,
    V = V,
    singular_values = sigma,
    U_rot = rot$U_rot,
    V_rot = rot$V_rot,
    R = rot$R,
    x_scores = x_scores,
    y_scores = y_scores,
    units = units,
    texts = list(.result_text(s, rows$keep, rows$x, rows$n_no_vector)),
    fit_args = fit_args
  ), class = c("ssd_multipls", "ssd_result"))
}

# ---------------------------------------------------------------------------
# Helpers of the fits
# ---------------------------------------------------------------------------

# y and every column of x with the formula's fixed effects regressed out by
# OLS: one QR decomposition for all D + 1 columns. `data` holds the kept rows.
.pls_residualize <- function(formula, data, y, x) {
  rhs <- formula[[3L]]
  if (.has_bar(rhs)) {
    stop("fit_pls() takes fixed covariates only, but the formula has random effects. Use fit_regression() for random effects.",
         call. = FALSE)
  }
  unknown <- setdiff(all.vars(rhs), names(data))
  if (length(unknown)) {
    stop(sprintf("fit_pls(): the formula uses %s, which data does not have.", .quoted(unknown)), call. = FALSE)
  }
  .check_nobs(formula, data)
  X <- stats::model.matrix(stats::as.formula(call("~", rhs), env = environment(formula)), data)
  r <- unname(qr.resid(qr(X), cbind(y, x)))
  list(y = r[, 1L], x = r[, -1L, drop = FALSE])
}

# beta flipped so that a higher alignment (x %*% beta) means a higher y.
.orient_beta <- function(beta, ys, x) {
  yhat <- drop(x %*% beta)
  corr <- if (.sd_pop(yhat) > 0) stats::cor(ys, yhat) else 0
  if (is.finite(corr) && corr < 0) -beta else beta
}

# ---------------------------------------------------------------------------
# Rotation (ssdiff/backends/pls.py)
# ---------------------------------------------------------------------------

.NO_ROTATION_META <- list(sweeps = 0L, V_converged = 0, kaiser_normalized = FALSE)

# The varimax target of rotate = "words": the leading 50,000 rows of the
# pack, centred and scaled as the document vectors were, projected onto W
# (D x k). Python computes it in float32, R in double.
.word_target <- function(emb, W, x_mean, x_scale) {
  M <- W / x_scale
  L <- emb$vectors[seq_len(min(.RAM_TOP_N, nrow(emb$vectors))), , drop = FALSE] %*% M
  sweep(L, 2L, drop(x_mean %*% M))
}

# Orthogonal R (k x k) from plskit's varimax of W against the target L, and
# its diagnostics.
.varimax <- function(W, L) {
  spec <- plskit::rotate(W, method = "varimax", L = L)$spec
  list(R = spec$R, meta = list(sweeps = as.integer(spec$sweeps), V_converged = spec$V_converged,
                               kaiser_normalized = spec$args$kaiser_normalize))
}

# Orthonormal weights W (D x N) rotated and scored (pls1_rotate):
# rotate = "words" runs varimax against L, orders the dimensions by
# |corr(xs %*% W_rot[, i], ys)| and flips each to a positive correlation;
# rotate = NULL keeps PLS order and signs. cov_share[i] is
# (W_rot[, i]' xs' ys)^2 over its sum.
.pls1_rotate <- function(W, xs, ys, rotate, L = NULL) {
  k <- ncol(W)
  v <- if (is.null(rotate)) list(R = diag(k), meta = .NO_ROTATION_META) else .varimax(W, L)
  W_rot <- W %*% v$R
  score_corr <- vapply(seq_len(k), function(i) stats::cor(drop(xs %*% W_rot[, i]), ys), 0)
  if (is.null(rotate)) {
    ord <- seq_len(k)
    signs <- rep(1, k)
  } else {
    ord <- order(-abs(score_corr), method = "radix")
    signs <- ifelse(score_corr[ord] < 0, -1, 1)
  }
  W_rot <- sweep(W_rot[, ord, drop = FALSE], 2L, signs, "*")
  c2 <- drop(crossprod(W_rot, crossprod(xs, ys)))^2
  list(W_rot = W_rot, R = sweep(v$R[, ord, drop = FALSE], 2L, signs, "*"), order = ord, signs = signs,
       score_corr = score_corr[ord] * signs, cov_share = c2 / sum(c2), meta = v$meta)
}

# A PLS3 pair (U, V) rotated by one orthogonal R (pls3_rotate). c[i], the
# squared cross-covariance of dimension i, orders the rotated modes; every
# mode signs each column so that its largest |V_rot| entry is positive.
# rotate = "outcomes" with L = NULL targets V, and takes R = t(V) when V is square.
.pls3_rotate <- function(U, V, sigma, rotate, L = NULL) {
  k <- ncol(U)
  v <- if (is.null(rotate)) {
    list(R = diag(k), meta = .NO_ROTATION_META)
  } else if (rotate == "outcomes" && is.null(L) && nrow(V) == k) {
    list(R = t(V), meta = .NO_ROTATION_META)
  } else {
    .varimax(U, if (is.null(L)) V else L)
  }
  R <- v$R
  c2 <- colSums(R^2 * sigma^2)
  ord <- if (is.null(rotate)) seq_len(k) else order(-c2, method = "radix")
  U_rot <- (U %*% R)[, ord, drop = FALSE]
  V_rot <- (V %*% R)[, ord, drop = FALSE]
  R <- R[, ord, drop = FALSE]
  idx <- apply(abs(V_rot), 2L, which.max)
  signs <- ifelse(V_rot[cbind(idx, seq_len(k))] < 0, -1, 1)
  list(U_rot = sweep(U_rot, 2L, signs, "*"), V_rot = sweep(V_rot, 2L, signs, "*"), R = sweep(R, 2L, signs, "*"),
       c = c2[ord], order = ord, signs = signs, meta = v$meta)
}

# The rotation diagnostics as a one-row data frame; `order` (1-based) and
# `signs` are list columns.
.pls_info <- function(n, rotate, rot, random_state) {
  data.frame(ssd_gradients = as.integer(n), rotate = if (is.null(rotate)) NA_character_ else rotate,
             order = I(list(rot$order)), signs = I(list(rot$signs)),
             kaiser_normalized = rot$meta$kaiser_normalized, sweeps = rot$meta$sweeps,
             V_converged = rot$meta$V_converged, random_state = random_state, stringsAsFactors = FALSE)
}

# The gradient units "dim-1" ... of orthonormal weights W, rotated by
# .pls1_rotate(), and their pls_info. Each unit's beta is W_rot[, i] / x_scale.
.pls_gradients <- function(W, xs, ys, x_scale, rotate, L, random_state) {
  rot <- .pls1_rotate(W, xs, ys, rotate, L)
  units <- lapply(seq_len(ncol(W)), function(i) {
    beta <- rot$W_rot[, i] / x_scale
    list(beta = beta, gradient = .unit_vector(beta),
         stats = data.frame(cov_share = rot$cov_share[i], score_corr = rot$score_corr[i]))
  })
  names(units) <- paste0("dim-", seq_along(units))
  list(units = units, pls_info = .pls_info(ncol(W), rotate, rot, random_state))
}

# ---------------------------------------------------------------------------
# fit_regression(ssd_gradients > 1)
# ---------------------------------------------------------------------------

# Stops unless the gradients can be built: plskit loads, and the formula has
# one text with one plain ssd() term.
.check_regression_gradients <- function(several, labels, text_names, full) {
  need_pkg("plskit", "fit_regression(ssd_gradients > 1)")
  ph <- .placeholder(text_names[1L])
  if (several || !identical(.text_labels(labels, ph), ph)) {
    stop(sprintf(paste0("fit_regression(ssd_gradients > 1) needs the text as one plain ssd() term, but the formula ",
                        "is '%s'. A moderated term's beta is a slope difference, which corr(score, y) does not ",
                        "describe."), deparse1(full)), call. = FALSE)
  }
}

# The n rotated gradients of one text term from its fitted component weights
# `w` (K): the K components (p$components, orthonormal rows in standardized
# space) are split into n consecutive blocks, the first ones a component
# longer when n does not divide K. Block j's piece t(components[b, ]) %*%
# w[b] is the term's beta restricted to those components; the pieces sum to
# it and are mutually orthogonal, so their unit vectors are orthonormal
# weights, rotated by words as in fit_pls(). `p` is fit_regression's part of
# the text, `ys` the standardized y.
.regression_gradients <- function(p, w, ys, n) {
  K <- length(w)
  if (n > K) {
    stop(sprintf("ssd_gradients=%d is more than the K=%d PCA components of the fit; use ssd_gradients <= K.", n, K),
         call. = FALSE)
  }
  sizes <- rep(K %/% n, n) + (seq_len(n) <= K %% n)
  blocks <- split(seq_len(K), rep(seq_len(n), sizes))
  pieces <- vapply(blocks, function(b) drop(crossprod(p$components[b, , drop = FALSE], w[b])),
                   numeric(ncol(p$components)))
  total <- sqrt(sum(rowSums(pieces)^2))
  norms <- sqrt(colSums(pieces^2))
  for (j in seq_len(n)) {
    b <- blocks[[j]]
    comps <- if (length(b) > 1L) sprintf("components %d-%d", b[1L], b[length(b)]) else sprintf("component %d", b)
    if (!all(is.finite(pieces[, j]))) {
      stop(sprintf("ssd_gradients: gradient %d (%s) has a coefficient the model could not estimate (aliased).",
                   j, comps), call. = FALSE)
    }
    if (isTRUE(norms[j] <= 1e-12 * total)) {
      stop(sprintf(paste0("ssd_gradients: gradient %d (%s) is zero: the model gives those components no weight. ",
                          "Use a smaller ssd_gradients."), j, comps), call. = FALSE)
    }
  }
  W <- unname(sweep(pieces, 2L, norms, "/"))
  L <- .word_target(.embeddings_from_ref(p$s$pack), W, p$x_mean, p$x_scale)
  .pls_gradients(W, .standardize(p$x)$x, ys, p$x_scale, "words", L, NA)
}
