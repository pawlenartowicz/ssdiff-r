# Port of _fit_regression, _wald_row, _compare_row, _whole_model_test,
# SSD.fit_regression and the module-level fit_regression in ssdiff/ssd.py, and
# of the K sweep in ssdiff/backends/pca_sweep.py and _sweep_math.py (with the
# fields of PCAOLSResult in ssdiff/results/continuous_result.py). Change them
# together.

#' Regress an outcome on texts, with covariates
#'
#' Fits a regression of the outcome on the first K principal components of the
#' text's document vectors, together with the covariates of a formula, and
#' turns the components' coefficients back into a semantic gradient (`beta`):
#' the direction in embedding space along which the outcome changes.
#'
#' @section The `ssd()` term:
#' The text is a term of the formula, written `ssd()`: `score ~ age +
#' ssd()`. A formula without an `ssd()` term gets `+ ssd()` appended, and
#' without a formula the model is the intercept plus the text. The term can
#' interact with covariates: in `score ~ age + sex * ssd()` (treatment
#' coding), `ssd()` is the gradient at the reference level of `sex`, and
#' `sexM:ssd()` the difference between level `M` and the reference, so the
#' gradient at `M` is the sum of the two. `score ~ sex:ssd()` has one gradient
#' per level, and in `score ~ age * ssd()` `ssd()` is the gradient at `age =
#' 0`. Each term is then a unit of the result (`res$units`, keyed
#' `"sexM:ssd()"`, ...). A formula with `ssd(...)` is used as written, so
#' `~ sex:ssd()` has no main text term. `ssd()` cannot appear inside a
#' random-effect term.
#'
#' The term is never evaluated: `fit_regression()` reads it from the formula
#' and replaces it with `(.ssd_pca_1 + ... + .ssd_pca_K)` before a model
#' function sees the formula. Names starting with `.ssd_` are reserved.
#'
#' @section Several texts:
#' `s` can be a named list of texts, each named in the formula as
#' `ssd(name)`: `fit_regression(list(essay = s1, diary = s2), score ~ age +
#' ssd(essay) + sex * ssd(diary), data = df)`. A name that is not an element
#' of `s` is built from the column of that name in `data` as
#' `ssd(data[[name]])`, with every default and one [load_wordpack()] for all
#' such texts. When every text comes from a column, the formula can come
#' first: `fit_regression(score ~ ssd(essay), data = df)`. A row stays when it
#' has a valid outcome and a vector in every text. Each text has its own
#' standardization, PCA and K (`fixed_k` applies to all).
#'
#' @section Choosing K:
#' Unless `fixed_k` is given, K is chosen once per text by a sweep over
#' `seq(k_min, k_max, k_step)`: OLS of the standardized outcome on the
#' standardized components, whatever the family, scoring each K by the
#' interpretability of the clusters of neighbours of both poles and by the
#' stability of the gradient, and picking the best joint score. The per-K
#' table is `res$texts[[i]]$sweep`.
#'
#' @section Two model paths:
#' A model without random effects, Gaussian with the identity link, is fitted
#' with [stats::lm()]; binomial and Poisson models without random effects with
#' [stats::glm()]. A formula with random effects (`(1 | g)`) or another family
#' or link is fitted with the fastglmm package, which must be installed. Both
#' paths give the numbers of Python's ssdiff.
#'
#' fastglmm takes column names of letters, digits and underscores only,
#' starting with a letter: a formula column such as `age.years` raises on
#' that path (rename it), and `data` columns starting with `ssd__` are
#' reserved. There the text's columns are passed as `ssd__pca_1`, ..., so the
#' coefficient names of `res$model` use `ssd__` where everything else in the
#' result shows `.ssd_`.
#'
#' @section Tests:
#' `res$test` holds a refit test of all text terms (an F test for OLS, else a
#' likelihood-ratio test; on fastglmm as fastglmm reports it), with several
#' texts one refit test per text, and with a text that has several terms a
#' Wald test per term label (`sex:ssd()` tests whether the gradient differs by
#' `sex`). `res$terms` has a Wald test per unit, and `res$model_test` the
#' whole model against `outcome ~ 1`.
#'
#' @param s A text from [ssd()]; a named list of texts; or the formula, when
#'   every text comes from a column of `data` (then leave out `y`).
#' @param y The outcome: numeric values, one per document, a column name of
#'   `data`, or a formula whose left side names the outcome column, such as
#'   `score ~ age + sex * ssd() + (1 | school)`. Documents with a missing
#'   outcome drop out.
#' @param data A data frame with one row per document, or `NULL`.
#' @param family A family as [stats::glm()] takes it: a family object
#'   (`binomial(link = "probit")`), a family function or its name.
#' @param fixed_k A fixed number of components, or `NULL` for the sweep.
#' @param k_min,k_max,k_step The K values the sweep tries.
#' @param ssd_gradients The number of gradients, from 1 to K. With `N > 1`
#'   (one text, one plain `ssd()` term), the K components are split into `N`
#'   consecutive blocks; each block's part of `beta` is one gradient, and the
#'   `N` gradients are varimax-rotated towards the pack's most frequent words
#'   as in [fit_pls()]. They are the units `"dim-1"`, ..., with `beta`,
#'   `gradient` and `stats` (`cov_share`, `score_corr`), and `pls_info` holds
#'   the rotation's diagnostics. The standard fit does not change. Needs plskit.
#' @param verbose `TRUE` shows the sweep's progress.
#' @param sweep_args A named list of sweep settings: `cluster_topn` (100),
#'   `cluster_k_min` (2), `cluster_k_max` (5), `weight_by_size` (`TRUE`),
#'   `auck_radius` (3).
#' @return An object of class `c("ssd_pca_ols", "ssd_result")`, a list with:
#'   \describe{
#'   \item{`stats`}{a one-row data frame: `backend` (`"PCA+OLS"` for a
#'     Gaussian identity model without random effects, else `"PCA+GLMM"`),
#'     `r2`, `r2_adj` (the whole model's; `NA` for GLMs and mixed models),
#'     `pvalue` (the test of all text terms), `n_raw`, `n_kept`, `n_dropped`,
#'     `y_mean`, `y_std`, `beta_norm`, `delta`, `iqr_effect`, `y_corr_pred`;}
#'   \item{`fit_info`}{a one-row data frame of the fit's settings, K
#'     (`n_components`) included;}
#'   \item{`test`}{one row per test: `label`, `name`, `pvalue`, `statistic`,
#'     `df`, `df_resid`, `r2_delta`;}
#'   \item{`beta`, `gradient`}{the gradient (on the standardized-outcome
#'     scale for a Gaussian identity model, else on the link scale) and its
#'     unit vector; `NULL` with several units;}
#'   \item{`y`, `alignment_scores`, `y_hat`}{the outcome of the rows used,
#'     their cosine with the gradient, and the predicted outcome;}
#'   \item{`outcome_name`, `formula`, `pca_weights`}{the outcome's name, the
#'     formula as written (`NULL` without one), and the component weights;}
#'   \item{`model`}{the `lm`, `glm` or `fastglmm` fit;}
#'   \item{`covariates`}{the covariates' estimates with the text held fixed:
#'     `term`, `estimate`, `se`, `statistic`, `p_value`, `ci_low`, `ci_high`;}
#'   \item{`terms`}{a Wald test per unit: `term`, `name`, `statistic`, `df`,
#'     `df_resid`, `pvalue`;}
#'   \item{`model_test`}{`fvalue`, `f_pvalue` (OLS) or `llr`, `llr_pvalue`:
#'     the whole model against `outcome ~ 1`;}
#'   \item{`units`}{with several text terms, one entry per term, keyed as
#'     `"sexM:ssd(essay)"`, with its own `beta`, `gradient`, `pca_weights`,
#'     `alignment_scores`, `stats`, `fit_info`, `text` (its position in
#'     `texts`) and `test` (the label of the `test` row of its term); with
#'     `ssd_gradients > 1`, one entry per gradient, keyed `"dim-1"`, ...,
#'     with `beta`, `gradient` and `stats` (`cov_share`, `score_corr`); else
#'     `NULL`;}
#'   \item{`pls_info`}{with `ssd_gradients > 1`, the rotation's diagnostics
#'     (`ssd_gradients`, `rotate`, `order`, `signs`, `kaiser_normalized`,
#'     `sweeps`, `V_converged`, `random_state`); else `NULL`;}
#'   \item{`texts`}{one entry per text: its tokens, settings, the vectors and
#'     rows used, `k`, `pca_components` and `sweep`;}
#'   \item{`obs_names`}{the names of the rows used: the row names of `data`,
#'     or the row numbers without `data`;}
#'   \item{`fit_args`}{the formula (as a string) and the other arguments.}
#'   }
#'   The model generics [coef()], [vcov()], [confint()], [nobs()],
#'   [logLik()], [residuals()], [fitted()], [df.residual()] and
#'   [ranef()][ranef.ssd_pca_ols] read the model, without the text's
#'   coefficients.
#' @examples
#' \dontrun{
#' s <- ssd(df$essay, lexicon = c("work", "job"))
#' r <- fit_regression(s, "score", data = df)
#' r <- fit_regression(s, score ~ age + sex * ssd(), data = df)
#' r$test
#' coef(r)
#' r <- fit_regression(score ~ ssd(essay) + ssd(diary), data = df)
#' if (requireNamespace("fastglmm", quietly = TRUE)) {
#'   r <- fit_regression(s, score ~ age + (1 | school), data = df)
#' }
#' }
#' @importFrom stats gaussian
#' @export
fit_regression <- function(s, y, data = NULL, family = gaussian(), fixed_k = NULL, k_min = 2, k_max = 120,
                           k_step = 2, ssd_gradients = 1, verbose = FALSE, sweep_args = NULL) {
  if (inherits(s, "formula")) {
    if (!missing(y) && is.data.frame(y)) {
      stop(sprintf("fit_regression(): pass data = df by name: fit_regression(%s, data = df).", deparse1(s)),
           call. = FALSE)
    }
    if (!missing(y)) {
      stop("fit_regression(): s is a formula, so leave out y, or pass the texts first: ",
           "fit_regression(list(), formula, data = df).", call. = FALSE)
    }
    y <- s
    s <- list()
  }
  if (!is.numeric(ssd_gradients) || length(ssd_gradients) != 1L || is.na(ssd_gradients) ||
      ssd_gradients < 1 || ssd_gradients != round(ssd_gradients)) {
    stop(sprintf("ssd_gradients must be a whole number >= 1; got %s", deparse1(ssd_gradients)), call. = FALSE)
  }
  .check_flag(verbose, "verbose")
  family <- .regression_family(family)
  if (!is.null(fixed_k)) .check_count(fixed_k, "fixed_k")
  for (a in c("k_min", "k_max", "k_step")) .check_count(get(a), a)
  if (k_max < k_min) stop(sprintf("k_max (%d) must be >= k_min (%d)", k_max, k_min), call. = FALSE)
  sweep_args <- .check_sweep_args(sweep_args)

  if (inherits(s, "ssd")) {
    out <- .outcome_values(y, data, s$n_docs, "regression")
    full <- if (is.null(out$formula)) {
      stats::as.formula(call("~", as.name(out$outcome_name), quote(ssd())), env = baseenv())
    } else {
      .with_text_term(out$formula)
    }
    text_names <- .text_names(full)
    if (length(text_names) > 1L) {
      shown <- paste(ifelse(is.na(text_names), "ssd()", sprintf("ssd(%s)", text_names)), collapse = ", ")
      stop(sprintf(paste0("fit_regression() got one text, but the formula has %s. ",
                          "Pass several texts as a named list, e.g. fit_regression(list(essay = s1, diary = s2), ...)."),
                   shown), call. = FALSE)
    }
    texts <- list(s)
  } else if (is.list(s)) {
    if (!inherits(y, "formula") || length(y) != 3L) {
      stop("fit_regression() with a list of texts needs a formula with the outcome on its left side, ",
           "e.g. score ~ ssd(essay).", call. = FALSE)
    }
    text_names <- .text_names(y)
    if (!length(text_names)) {
      stop(sprintf(paste0("fit_regression(): the formula %s has no ssd(name) term; write one for each text, ",
                          "e.g. y ~ ssd(essay) + ssd(diary)."), deparse1(y)), call. = FALSE)
    }
    if (anyNA(text_names)) {
      stop("fit_regression(): bare ssd() needs one text from ssd() as s; with a list of texts, ",
           "write ssd(name) and pass list(name = text).", call. = FALSE)
    }
    if (length(s) && (is.null(names(s)) || any(!nzchar(names(s))))) {
      stop("fit_regression(): every element of the list s needs a name, the name its ssd(name) term uses.",
           call. = FALSE)
    }
    for (nm in names(s)) .check_ssd(s[[nm]], sprintf("s$%s", nm))
    unused <- setdiff(names(s), text_names)
    if (length(unused)) {
      stop(sprintf("fit_regression(): the texts %s were passed but the formula names only %s.",
                   .quoted(unused), .quoted(text_names)), call. = FALSE)
    }
    if (is.null(data)) stop("fit_regression() needs data = df.", call. = FALSE)
    if (!is.data.frame(data)) stop("data must be a data frame", call. = FALSE)
    absent <- setdiff(text_names, c(names(s), names(data)))
    if (length(absent)) {
      stop(sprintf(paste0("fit_regression(): the formula names the texts %s, which are neither columns of data ",
                          "nor elements of s. Add a text column to data, or pass s = list(%s = ssd(...))."),
                   .quoted(absent), absent[1L]), call. = FALSE)
    }
    columns <- setdiff(text_names, names(s))
    if (length(columns)) {
      emb <- load_wordpack()
      for (nm in columns) {
        col <- data[[nm]]
        s[[nm]] <- ssd(if (is.factor(col)) as.character(col) else col, embeddings = emb)
      }
    }
    sizes <- vapply(s, function(t) as.integer(t$n_docs), 1L)
    if (length(unique(sizes)) != 1L) {
      stop(sprintf("fit_regression(): the texts differ in length: %s.",
                   paste(sprintf("%s = %d", names(sizes), sizes), collapse = ", ")), call. = FALSE)
    }
    texts <- unname(s[text_names])
    out <- .outcome_values(y, data, sizes[[1L]], "regression")
    full <- y
  } else {
    stop("s must be a text from ssd(), a named list of texts, or the formula", call. = FALSE)
  }
  fit_args <- list(formula = if (is.null(out$formula)) NULL else deparse1(out$formula),
                   family = family$family, link = family$link, fixed_k = fixed_k, k_min = k_min, k_max = k_max,
                   k_step = k_step, ssd_gradients = ssd_gradients, sweep_args = sweep_args)
  .fit_regression(out, full, data, texts, text_names, family, fixed_k, k_min, k_max, k_step, verbose,
                  sweep_args, fit_args)
}

# The family object of a family, a family function or its name (as glm()).
.regression_family <- function(family) {
  if (is.character(family) && length(family) == 1L) {
    family <- tryCatch(get(family, mode = "function"), error = function(e) {
      stop(sprintf("family '%s' is not a family function", family), call. = FALSE)
    })
  }
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) {
    stop("family must be a family object such as gaussian() or binomial(link = \"probit\"), a family function, or its name",
         call. = FALSE)
  }
  family
}

# The sweep's tuning keys (Python's pca_sweep_args) and their defaults.
.SWEEP_DEFAULTS <- list(cluster_topn = 100L, cluster_k_min = 2L, cluster_k_max = 5L, weight_by_size = TRUE,
                        auck_radius = 3L)

.check_sweep_args <- function(sweep_args) {
  if (is.null(sweep_args)) return(NULL)
  if (!is.list(sweep_args) || (length(sweep_args) && (is.null(names(sweep_args)) || any(!nzchar(names(sweep_args)))))) {
    stop("sweep_args must be a named list", call. = FALSE)
  }
  bad <- setdiff(names(sweep_args), names(.SWEEP_DEFAULTS))
  if (length(bad)) {
    stop(sprintf("sweep_args takes only the keys %s; got %s.", .quoted(names(.SWEEP_DEFAULTS)), .quoted(bad)),
         call. = FALSE)
  }
  sweep_args
}

# ---------------------------------------------------------------------------
# The fit
# ---------------------------------------------------------------------------

# `out` is .outcome_values()'s list over all documents, `full` the formula with
# its ssd() terms, `texts` the ssd objects in formula order and `text_names`
# the names inside their ssd(...) (NA for ssd()).
.fit_regression <- function(out, full, data, texts, text_names, family, fixed_k, k_min, k_max, k_step, verbose,
                            sweep_args, fit_args) {
  if ("y" %in% text_names) {
    stop("fit_regression(): the text name 'y' is reserved for the outcome; pick another name.", call. = FALSE)
  }
  .check_text_terms(full)
  several <- length(texts) > 1L
  f_ph <- .to_placeholders(full)
  f_ph[[2L]] <- quote(.ssd_y)
  labels <- .formula_terms(f_ph)$labels
  n_grad <- as.integer(fit_args$ssd_gradients)
  if (n_grad > 1L) .check_regression_gradients(several, labels, text_names, full)

  valid <- out$valid
  keep <- valid
  for (t in texts) {
    .require_vectors(t)
    keep <- keep & t$keep
  }
  n_raw <- length(keep)
  n_kept <- sum(keep)
  if (n_kept < 3L) {
    stop(sprintf("fit_regression(): only %d documents have a valid outcome and a vector in every text.", n_kept),
         call. = FALSE)
  }
  y <- out$values[keep]
  data_kept <- NULL
  obs_names <- as.character(which(keep))
  if (!is.null(data)) {
    obs_names <- row.names(data)[keep]
    data_kept <- data[keep, , drop = FALSE]
    rownames(data_kept) <- NULL
    unknown <- setdiff(.vars_outside_ssd(full[[3L]]), names(data))
    if (length(unknown)) {
      stop(sprintf("fit_regression(): the formula uses %s, which data does not have.", .quoted(unknown)),
           call. = FALSE)
    }
    .check_nobs(full, data_kept)
  } else if (length(.vars_outside_ssd(full[[3L]]))) {
    stop(sprintf("fit_regression(): the formula uses %s, but there is no data.",
                 .quoted(.vars_outside_ssd(full[[3L]]))), call. = FALSE)
  }
  # The path needs only the formula and family: decided, and fastglmm checked,
  # before the sweep.
  path <- .regression_path(full, family)
  if (path == "fastglmm") .check_fastglmm_names(full, data)

  ys <- drop(.standardize(y)$x)
  parts <- lapply(seq_along(texts), function(i) {
    t <- texts[[i]]
    name <- text_names[i]
    x <- .x_rows(t, keep)
    st <- .standardize(x)
    sweep <- NULL
    if (is.null(fixed_k)) {
      sweep <- .pca_sweep(st$x, st$scale, x, ys, .embeddings_from_ref(t$pack), seq(k_min, k_max, k_step),
                          sweep_args, verbose)
      k <- sweep$best_k
    } else {
      k <- as.integer(fixed_k)
    }
    k <- min(k, nrow(st$x), ncol(st$x))
    pca <- .pca_fit_transform(st$x, k)
    prefix <- if (several) paste0(".ssd_", name, "_pca_") else ".ssd_pca_"
    list(name = name, s = t, x = x, x_mean = st$mean, x_scale = st$scale, sweep = sweep, z = pca$z,
         components = pca$components, pcs = paste0(prefix, seq_len(k)), ph = .placeholder(name),
         shown = if (is.na(name)) "ssd()" else sprintf("ssd(%s)", name), n_no_vector = sum(valid & !t$keep))
  })

  cols <- if (is.null(data_kept)) data.frame(row.names = seq_len(n_kept)) else data_kept
  cols$.ssd_y <- y
  for (p in parts) for (i in seq_along(p$pcs)) cols[[p$pcs[i]]] <- p$z[, i]

  phs <- vapply(parts, `[[`, "", "ph")
  for (lab in labels) {
    vars <- strsplit(lab, ":", fixed = TRUE)[[1L]]
    if (sum(phs %in% vars) > 1L) {
      shown <- Reduce(function(l, p) .label_shown(l, p$ph, p$shown), parts, lab)
      stop(sprintf("The term '%s' has two texts; two texts cannot interact in one term.", shown), call. = FALSE)
    }
  }
  own_labels <- lapply(parts, function(p) .text_labels(labels, p$ph))
  all_text_labels <- unlist(own_labels)
  columns <- stats::setNames(lapply(parts, `[[`, "pcs"), phs)
  h0 <- .null_formula(f_ph, all_text_labels)
  h1 <- .expand(f_ph, columns)

  cmp <- .reg_compare(path, cols, h0, h1, family)
  fit <- cmp$fit1
  coefs <- .reg_coef(fit)
  coef_names <- names(coefs)

  linear <- family$family == "gaussian" && family$link == "identity"
  y_st <- .standardize(y)
  y_mean <- if (linear) y_st$mean[[1L]] else 0
  y_scale <- if (linear) y_st$scale[[1L]] else 1

  for (i in seq_along(parts)) {
    p <- parts[[i]]
    comps <- lapply(seq_along(p$pcs), function(j) list(p$shown, j))
    parts[[i]]$groups <- .term_groups(coef_names, stats::setNames(comps, p$pcs))
  }
  n_units <- sum(vapply(parts, function(p) length(p$groups), 1L))

  # The all-text row, one row per text, then one Wald row per label of a text
  # that has several.
  tests <- list(c(list(label = NA_character_), cmp$row))
  wald <- list()
  term_label <- list()                       # key -> label of its Wald row
  for (i in seq_along(parts)) {
    p <- parts[[i]]
    if (several) {
      others <- columns[setdiff(names(columns), p$ph)]
      sub <- .reg_compare(path, cols, .expand(.null_formula(f_ph, own_labels[[i]]), others), h1, family)
      tests[[length(tests) + 1L]] <- c(list(label = p$name), sub$row)
    }
    shown_labels <- vapply(own_labels[[i]], .label_shown, "", ph = p$ph, shown = p$shown, USE.NAMES = FALSE)
    if (length(shown_labels) > 1L) {
      of_key <- vapply(names(p$groups), .label_of, "", labels_shown = shown_labels)
      for (lab in shown_labels) {
        pos <- unlist(lapply(p$groups[of_key %in% lab], `[[`, "pos"))
        wald[[length(wald) + 1L]] <- c(list(label = lab), .reg_wald(fit, pos))
        for (key in names(p$groups)[of_key %in% lab]) term_label[[key]] <- lab
      }
    }
  }
  test <- .rows_df(c(tests, wald), c("label", "name", "pvalue", "statistic", "df", "df_resid", "r2_delta"))
  terms <- .rows_df(unlist(lapply(parts, function(p) lapply(names(p$groups), function(key) {
    c(list(term = key), .reg_wald(fit, p$groups[[key]]$pos))
  })), recursive = FALSE), c("term", "name", "statistic", "df", "df_resid", "pvalue"))
  model_test <- .whole_model_test(path, cols, f_ph, labels, all_text_labels, cmp, family)

  r2 <- .reg_r2(fit)
  backend <- if (linear && !.has_bar(full[[3L]])) "PCA+OLS" else "PCA+GLMM"
  fitted_mean <- unname(.reg_fitted(fit))
  formula_str <- if (is.null(out$formula)) NA_character_ else deparse1(out$formula)

  units <- list()
  for (i in seq_along(parts)) {
    p <- parts[[i]]
    fit_info <- data.frame(
      n_components = length(p$pcs), p_at_k = NA_integer_, n_splits = NA_integer_, random_state = NA_integer_,
      k_min = if (is.null(p$sweep)) NA_integer_ else as.integer(k_min),
      k_max = if (is.null(p$sweep)) NA_integer_ else as.integer(k_max),
      k_step = if (is.null(p$sweep)) NA_integer_ else as.integer(k_step),
      best_k = if (is.null(p$sweep)) NA_integer_ else as.integer(p$sweep$best_k),
      pca_k_source = if (is.null(p$sweep)) "fixed" else "sweep", formula = formula_str,
      family = family$family, ssd_gradients = n_grad, rotate = if (n_grad > 1L) "words" else NA_character_,
      outcome = out$outcome_name,
      stringsAsFactors = FALSE)
    for (key in names(p$groups)) {
      if (anyNA(coefs[p$groups[[key]]$pos])) {
        stop(sprintf(paste0("fit_regression(): the coefficients of '%s' cannot all be estimated (aliased), ",
                            "e.g. a factor level with fewer documents than components. ",
                            "Merge or drop that level, or set a smaller fixed_k."), key), call. = FALSE)
      }
      w <- unname(coefs[p$groups[[key]]$pos]) / y_scale
      beta <- drop(t(p$components) %*% w) / p$x_scale
      shared <- n_units > 1L || grepl(":", key, fixed = TRUE)
      y_hat <- if (shared || !linear) fitted_mean else y_mean + y_scale * drop(sweep(p$x, 2L, p$x_mean) %*% beta)
      unit <- .regression_unit(beta, p$x, y, if (shared) fitted_mean else drop(p$x %*% beta), y_scale)
      unit$stats <- data.frame(backend = backend, r2 = r2[[1L]], r2_adj = r2[[2L]], pvalue = cmp$row$pvalue,
                               n_raw = n_raw, n_kept = n_kept, n_dropped = n_raw - n_kept, y_mean = mean(y),
                               y_std = y_st$scale[[1L]], unit$stats, stringsAsFactors = FALSE)
      units[[key]] <- c(unit[c("beta", "gradient")], list(pca_weights = w), unit["alignment_scores"],
                        list(y_hat = y_hat, stats = unit$stats, fit_info = fit_info, text = i,
                             test = if (is.null(term_label[[key]])) NA_character_ else term_label[[key]]))
    }
    parts[[i]]$fit_info <- fit_info
  }

  first <- units[[1L]]
  one <- n_units == 1L
  stats <- first$stats
  if (!one) stats[c("beta_norm", "delta", "iqr_effect")] <- NA_real_
  texts_out <- lapply(parts, function(p) {
    c(.result_text(p$s, keep, p$x, p$n_no_vector, if (is.na(p$name)) NULL else p$name),
      list(k = length(p$pcs), pca_components = p$components, sweep = p$sweep$table))
  })
  grad <- if (n_grad > 1L) .regression_gradients(parts[[1L]], first$pca_weights, ys, n_grad)
  structure(list(
    stats = stats,
    fit_info = first$fit_info,
    test = test,
    beta = if (one) first$beta,
    gradient = if (one) first$gradient,
    y = y,
    alignment_scores = if (one) first$alignment_scores,
    y_hat = first$y_hat,
    outcome_name = out$outcome_name,
    formula = .at_global(out$formula),
    pca_weights = if (one) first$pca_weights,
    model = .release_env(fit),
    covariates = .reg_covariates(fit, !is.null(out$formula)),
    terms = terms,
    model_test = model_test,
    pls_info = grad$pls_info,
    units = if (!one) lapply(units, function(u) u[setdiff(names(u), "y_hat")]) else grad$units,
    texts = texts_out,
    obs_names = obs_names,
    fit_args = fit_args
  ), class = c("ssd_pca_ols", "ssd_result"))
}

# beta, gradient, alignment_scores and the effect sizes (beta_norm, delta,
# iqr_effect, y_corr_pred) of one unit. `pred` is the prediction y_corr_pred
# correlates with y; `y_scale` the unit beta was fitted in.
.regression_unit <- function(beta, x, y, pred, y_scale) {
  gradient <- .unit_vector(beta)
  beta_norm <- sqrt(sum(beta^2))
  align <- .alignment(x, gradient)
  corr <- if (.sd_pop(y) * .sd_pop(pred) > 0) stats::cor(y, pred) else 0
  if (!is.finite(corr)) corr <- 0
  q <- .np_percentile(align, c(0.75, 0.25))
  list(beta = beta, gradient = gradient, alignment_scores = align,
       stats = data.frame(beta_norm = beta_norm, delta = 0.1 * beta_norm * y_scale,
                          iqr_effect = (q[1L] - q[2L]) * beta_norm * y_scale, y_corr_pred = abs(corr)))
}

# ---------------------------------------------------------------------------
# The two model paths: "lm", "glm" (base R) and "fastglmm"
# ---------------------------------------------------------------------------

.GLM_LINKS <- list(binomial = c("logit", "probit", "cloglog"), poisson = "log")

.regression_path <- function(full, family) {
  bar <- .has_bar(full[[3L]])
  fam <- family$family
  if (!bar && fam == "gaussian" && family$link == "identity") return("lm")
  if (!bar && fam %in% names(.GLM_LINKS) && family$link %in% .GLM_LINKS[[fam]]) return("glm")
  why <- if (bar) "a formula with random effects" else {
    default_link <- tryCatch(get(fam, mode = "function")()$link, error = function(e) NA)
    if (identical(default_link, family$link)) sprintf("family = %s", fam)
    else sprintf("family = %s(link = \"%s\")", fam, family$link)
  }
  need_pkg("fastglmm", "fit_regression()", why = why)
  "fastglmm"
}

# fastglmm's formula parser takes names of letters, digits and underscores,
# starting with a letter.
.fastglmm_ok <- function(x) grepl("(*UCP)^\\p{L}[\\p{L}\\p{N}_]*$", x, perl = TRUE)

# Stops before any fitting when a column the formula uses has a name fastglmm
# cannot parse, or when data has a column starting with `ssd__`, the prefix
# of fit_regression()'s own columns on the fastglmm path.
.check_fastglmm_names <- function(full, data) {
  bad <- setdiff(.vars_outside_ssd(full[[3L]]), "")
  bad <- bad[!.fastglmm_ok(bad)]
  if (length(bad)) {
    stop(sprintf(paste0("fit_regression() fits this model with fastglmm, which takes column names of letters, ",
                        "digits and underscores only, starting with a letter. Rename %s in data and the formula."),
                 .quoted(bad)), call. = FALSE)
  }
  cols <- if (is.null(data)) character(0) else names(data)
  clash <- sort(cols[startsWith(cols, "ssd__")], method = "radix")
  if (length(clash)) {
    stop(sprintf(paste0("data columns %s start with 'ssd__', which fit_regression() reserves for its own ",
                        "columns on the fastglmm path. Rename them."), .quoted(clash)), call. = FALSE)
  }
  invisible()
}

# The names fit_regression()'s own columns `.ssd_*` get on the fastglmm path:
# `ssd__` and the rest of the name when fastglmm parses it (`ssd__pca_1`),
# else `ssd__c<k>` (a text name with a dot).
.fastglmm_map <- function(names) {
  ours <- names[startsWith(names, ".ssd_")]
  rest <- substring(ours, 6L)
  stats::setNames(ifelse(.fastglmm_ok(rest), paste0("ssd__", rest), paste0("ssd__c", seq_along(ours))), ours)
}

# Coefficient names of a fastglmm fit with the `.ssd_*` names put back.
.from_fastglmm_names <- function(x, fit) {
  map <- attr(fit, "ssd_names")
  back <- stats::setNames(names(map), map)
  vapply(strsplit(x, ":", fixed = TRUE), function(p) {
    hit <- p %in% names(back)
    p[hit] <- back[p[hit]]
    paste(p, collapse = ":")
  }, "")
}

.fastglmm_formula <- function(f, map) .replace_symbols(f, lapply(map, as.name))

# fastglmm on `cols` with its `.ssd_*` columns renamed: compare() for two
# formulas, one maximum-likelihood fit for one. The fits keep the names in
# their attribute "ssd_names".
.fastglmm_call <- function(cols, formulas, family, ...) {
  map <- .fastglmm_map(names(cols))
  names(cols)[match(names(map), names(cols))] <- map
  fs <- lapply(formulas, .fastglmm_formula, map = map)
  out <- if (length(fs) == 2L) {
    fastglmm::compare(cols, fs[[1L]], fs[[2L]], family = family)
  } else {
    fastglmm::fastglmm(fs[[1L]], cols, family = family, REML = FALSE, anova_test = "Wald")
  }
  if (inherits(out, "fastglmm")) return(structure(out, ssd_names = map))
  out$fit0 <- structure(out$fit0, ssd_names = map)
  out$fit1 <- structure(out$fit1, ssd_names = map)
  out
}

# Fits f0 and f1 on `cols` and compares them as glmm's compare(): `fit0`,
# `fit1` and `row`, the test as a list of name, pvalue, statistic, df,
# df_resid and r2_delta.
.reg_compare <- function(path, cols, f0, f1, family) {
  if (path == "lm") {
    m0 <- .base_fit(path, f0, cols, family)
    m1 <- .base_fit(path, f1, cols, family)
    a <- stats::anova(m0, m1)
    row <- list(name = "f_test", pvalue = a[2L, "Pr(>F)"], statistic = a[2L, "F"], df = as.numeric(a[2L, "Df"]),
                df_resid = as.numeric(a[2L, "Res.Df"]),
                r2_delta = summary(m1)$r.squared - summary(m0)$r.squared)
  } else if (path == "glm") {
    m0 <- .base_fit(path, f0, cols, family)
    m1 <- .base_fit(path, f1, cols, family)
    a <- stats::anova(m0, m1, test = "LRT")
    row <- list(name = "lrt", pvalue = a[2L, "Pr(>Chi)"], statistic = a[2L, "Deviance"],
                df = as.numeric(a[2L, "Df"]), df_resid = NA_real_, r2_delta = NA_real_)
  } else {
    cmp <- .fastglmm_call(cols, list(f0, f1), family)
    m0 <- cmp$fit0
    m1 <- cmp$fit1
    tab <- cmp$table
    f_test <- "Pr(>F)" %in% names(tab)
    r2_delta <- m1$r_squared - m0$r_squared
    row <- list(name = if (f_test) "f_test" else "lrt", pvalue = tab[2L, if (f_test) "Pr(>F)" else "Pr(>Chisq)"],
                statistic = tab[2L, if (f_test) "F" else "Chisq"], df = as.numeric(tab[2L, "Df"]),
                df_resid = if (f_test) as.numeric(tab[2L, "Res.Df"]) else NA_real_,
                r2_delta = if (is.finite(r2_delta)) r2_delta else NA_real_)
  }
  list(fit0 = m0, fit1 = m1, row = row)
}

# One lm() or glm() fit, its call naming the formula.
.base_fit <- function(path, f, cols, family) {
  m <- if (path == "lm") stats::lm(f, data = cols) else stats::glm(f, family = family, data = cols)
  m$call$formula <- f
  m$call$data <- quote(data)
  m
}

# The formula `f` (or NULL) with the global environment as its own.
.at_global <- function(f) {
  if (!is.null(f)) environment(f) <- globalenv()
  f
}

# The fit with its formulas and terms pointing at the global environment, so
# that a saved result does not carry the frame it was fitted from. Every
# column the model reads is in its data.
.release_env <- function(fit) {
  if (!inherits(fit, "lm")) return(fit)
  attr(fit$terms, ".Environment") <- globalenv()
  attr(attr(fit$model, "terms"), ".Environment") <- globalenv()
  fit$call$formula <- .at_global(fit$call$formula)
  if (!is.null(fit$formula)) fit$formula <- .at_global(fit$formula)
  fit
}

# One maximum-likelihood fit of f on the path.
.reg_fit_ml <- function(path, f, cols, family) {
  if (path != "fastglmm") return(.base_fit(path, f, cols, family))
  .fastglmm_call(cols, list(f), family)
}

# Coefficients (NA where aliased), named with `.ssd_` column names.
.reg_coef <- function(fit) {
  .model_pkg(fit)
  if (inherits(fit, "fastglmm")) {
    b <- fastglmm::fixef(fit)
    b[as.logical(fit$aliased)] <- NA
    names(b) <- .from_fastglmm_names(names(b), fit)
    return(b)
  }
  stats::coef(fit)
}

.reg_vcov <- function(fit) {
  .model_pkg(fit)
  v <- stats::vcov(fit)
  if (inherits(fit, "fastglmm")) dimnames(v) <- rep(list(.from_fastglmm_names(rownames(v), fit)), 2L)
  v
}

.reg_fitted <- function(fit) {
  .model_pkg(fit)
  stats::fitted(fit)
}

.reg_resid <- function(fit) {
  .model_pkg(fit)
  if (inherits(fit, "lm") && !inherits(fit, "glm")) stats::residuals(fit) else stats::residuals(fit, type = "response")
}

# The whole model's r2 and r2_adj; NA for GLMs and mixed models.
.reg_r2 <- function(fit) {
  if (inherits(fit, "glm")) return(c(NA_real_, NA_real_))
  if (inherits(fit, "lm")) {
    s <- summary(fit)
    return(c(s$r.squared, s$adj.r.squared))
  }
  r <- c(fit$r_squared, fit$adj_r_squared)
  r[!is.finite(r)] <- NA_real_
  r
}

# The Wald test of the coefficients at positions `pos` on the fit, as a
# list of name, statistic, df, df_resid and pvalue; NA when one is aliased.
.reg_wald <- function(fit, pos) {
  q <- length(pos)
  if (inherits(fit, "fastglmm")) {
    L <- matrix(0, q, length(fit$aliased))
    L[cbind(seq_len(q), pos)] <- 1
    ct <- fastglmm::contrast_test(fit, L)
    if ("F value" %in% names(ct)) {
      return(list(name = "wald_f", statistic = .num(ct[1L, "F value"]), df = .num(ct[1L, "NumDF"]),
                  df_resid = .num(ct[1L, "DenDF"]), pvalue = .num(ct[1L, "Pr(>F)"])))
    }
    return(list(name = "wald_chisq", statistic = .num(ct[1L, "Chisq"]), df = .num(ct[1L, "Df"]),
                df_resid = NA_real_, pvalue = .num(ct[1L, "Pr(>Chisq)"])))
  }
  is_lm <- !inherits(fit, "glm")
  b <- stats::coef(fit)[pos]
  na <- list(name = if (is_lm) "wald_f" else "wald_chisq", statistic = NA_real_, df = NA_real_,
             df_resid = NA_real_, pvalue = NA_real_)
  if (anyNA(b)) return(na)
  v <- stats::vcov(fit)[pos, pos, drop = FALSE]
  w <- tryCatch(drop(crossprod(b, solve(v, b))), error = function(e) NA_real_)
  if (!is.finite(w)) return(na)
  if (is_lm) {
    dfr <- stats::df.residual(fit)
    return(list(name = "wald_f", statistic = w / q, df = q, df_resid = dfr,
                pvalue = stats::pf(w / q, q, dfr, lower.tail = FALSE)))
  }
  list(name = "wald_chisq", statistic = w, df = q, df_resid = NA_real_,
       pvalue = stats::pchisq(w, q, lower.tail = FALSE))
}

.num <- function(v) if (length(v) && !is.null(v)) as.numeric(v) else NA_real_

# The covariates' rows of the fit (its coefficients without the text's), with
# Wald intervals: `term`, `estimate`, `se`, `statistic`, `p_value`, `ci_low`,
# `ci_high`. Empty without a formula and when only the intercept is left.
.reg_covariates <- function(fit, has_formula, level = 0.95) {
  empty <- data.frame(term = character(0), estimate = numeric(0), se = numeric(0), statistic = numeric(0),
                      p_value = numeric(0), ci_low = numeric(0), ci_high = numeric(0), stringsAsFactors = FALSE)
  if (!has_formula) return(empty)
  tab <- .reg_coef_table(fit, level)
  tab <- tab[!.is_text_coef(tab$term), , drop = FALSE]
  if (identical(tab$term, "(Intercept)")) return(empty)
  rownames(tab) <- NULL
  tab
}

# Every coefficient's estimate, se, statistic, p-value and Wald interval: t
# on the residual df for lm, z for glm, as fastglmm's summary() for fastglmm.
.reg_coef_table <- function(fit, level = 0.95) {
  .model_pkg(fit)
  if (inherits(fit, "fastglmm")) {
    # By name: a mixed model's table has a df column before the t value.
    sm <- summary(fit)$coefficients
    stat <- .coef_column(sm, "^[tz] value$", "statistic")
    p <- .coef_column(sm, "^Pr\\(", "p-value")
    ci <- stats::confint(fit, level = level)
    return(data.frame(term = .from_fastglmm_names(rownames(sm), fit), estimate = sm[, "Estimate"],
                      se = sm[, "Std. Error"], statistic = sm[, stat], p_value = sm[, p], ci_low = ci[, 1L],
                      ci_high = ci[, 2L],
                      stringsAsFactors = FALSE, row.names = NULL))
  }
  b <- stats::coef(fit)
  se <- rep(NA_real_, length(b))
  ok <- !is.na(b)
  se[ok] <- sqrt(diag(stats::vcov(fit))[names(b)[ok]])
  stat <- b / se
  a <- (1 - level) / 2
  if (inherits(fit, "glm")) {
    p <- 2 * stats::pnorm(-abs(stat))
    crit <- stats::qnorm(1 - a)
  } else {
    dfr <- stats::df.residual(fit)
    p <- 2 * stats::pt(-abs(stat), dfr)
    crit <- stats::qt(1 - a, dfr)
  }
  data.frame(term = names(b), estimate = unname(b), se = se, statistic = unname(stat), p_value = unname(p),
             ci_low = unname(b - crit * se), ci_high = unname(b + crit * se), stringsAsFactors = FALSE)
}

# The one column of a coefficient table whose name matches `pattern`.
.coef_column <- function(sm, pattern, what) {
  j <- grep(pattern, colnames(sm))
  if (length(j) != 1L) {
    stop(sprintf("fastglmm's coefficient table has no single %s column (its columns: %s). Please report this.",
                 what, paste(colnames(sm), collapse = ", ")), call. = FALSE)
  }
  j
}

# The whole model against `outcome ~ 1` with the same random effects and
# offsets: fvalue and f_pvalue for OLS, llr and llr_pvalue for every other
# model; all NA without an intercept.
.whole_model_test <- function(path, cols, f_ph, labels, text_labels, cmp, family) {
  out <- data.frame(fvalue = NA_real_, f_pvalue = NA_real_, llr = NA_real_, llr_pvalue = NA_real_)
  if (!.formula_terms(f_ph)$intercept) return(out)
  fit1 <- cmp$fit1
  every_text <- setequal(labels, text_labels)
  fit0 <- if (every_text) cmp$fit0 else .reg_fit_ml(path, .null_formula(f_ph, labels), cols, family)
  # As Python: the model decides (OLS: Gaussian, identity link, no random
  # effects), not the path that fitted it.
  ols <- family$family == "gaussian" && family$link == "identity" && !.has_bar(f_ph[[3L]])
  if (ols) {
    if (every_text) {
      out$fvalue <- cmp$row$statistic
      out$f_pvalue <- cmp$row$pvalue
    } else {
      a <- stats::anova(fit0, fit1)
      out$fvalue <- a[2L, "F"]
      out$f_pvalue <- a[2L, "Pr(>F)"]
    }
  } else {
    l1 <- stats::logLik(fit1)
    l0 <- stats::logLik(fit0)
    out$llr <- 2 * (as.numeric(l1) - as.numeric(l0))
    out$llr_pvalue <- stats::pchisq(out$llr, attr(l1, "df") - attr(l0, "df"), lower.tail = FALSE)
  }
  out
}

# ---------------------------------------------------------------------------
# The K sweep
# ---------------------------------------------------------------------------

.SWEEP_COLUMNS <- c("k", "var_explained", "mean_coherence", "mean_abs_cosb", "aggregate", "n_clusters",
                    "total_size", "beta_delta_1_minus_cos", "interp_hat", "interp_resid", "interp_resid_z",
                    "interp_auck", "stab_good_raw", "stab_z_raw", "stab_auck_raw", "joint_score")

# The sweep over the K values `ks` on standardized vectors `xs` (column scales
# `x_scale`), raw vectors `x` and standardized outcome `ys`. Returns `best_k`
# and `table` (one row per K, .SWEEP_COLUMNS).
.pca_sweep <- function(xs, x_scale, x, ys, emb, ks, sweep_args = NULL, verbose = FALSE) {
  a <- utils::modifyList(.SWEEP_DEFAULTS, as.list(sweep_args))
  n <- nrow(xs)
  D <- ncol(xs)
  sv <- svd(xs, nu = 0L)
  ev <- sv$d^2 / (n - 1)
  total <- sum(ev)

  # Per K: OLS on the first K components, beta oriented to correlate
  # positively with y, and its change from the previous K.
  recs <- vector("list", length(ks))
  prev <- NULL
  for (i in seq_along(ks)) {
    K <- ks[i]
    rec <- tryCatch({
      max_k <- min(K, n - 1L, D)
      if (max_k < 1L) stop(sprintf("K = %d is too large for the data (n = %d, D = %d)", K, n, D))
      comp <- sv$v[, seq_len(max_k), drop = FALSE]
      z <- xs %*% comp
      w <- solve(crossprod(z), crossprod(z, ys))
      beta <- drop(comp %*% w) / x_scale
      yhat <- drop(x %*% beta)
      corr <- if (.sd_pop(ys) * .sd_pop(yhat) > 0) stats::cor(ys, yhat) else 0
      if (!is.finite(corr)) corr <- 0
      if (corr < 0) beta <- -beta
      list(ok = TRUE, var = if (total > 0) sum(ev[seq_len(max_k)]) / total * 100 else 0,
           gradient = .unit_vector(beta))
    }, error = function(e) {
      if (verbose) message(sprintf("[sweep] K = %d skipped: %s", K, conditionMessage(e)))
      list(ok = FALSE)
    })
    if (rec$ok) {
      rec$delta <- if (is.null(prev)) NaN else 1 - .cosine(prev, rec$gradient)
      prev <- rec$gradient
    } else {
      prev <- NULL
    }
    recs[[i]] <- rec
  }

  # Interpretability: clusters of the neighbours of both poles.
  bar <- if (verbose) utils::txtProgressBar(min = 0, max = length(ks), style = 3) else NULL
  rows <- lapply(seq_along(ks), function(i) {
    rec <- recs[[i]]
    if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
    if (!rec$ok) {
      return(c(k = ks[i], var_explained = NaN, mean_coherence = NaN, mean_abs_cosb = NaN, aggregate = NaN,
               n_clusters = 0, total_size = 0, beta_delta_1_minus_cos = NaN))
    }
    clusters <- list()
    for (side in c("pos", "neg")) {
      cl <- tryCatch(
        .cluster_top_neighbors(emb, rec$gradient, topn = a$cluster_topn, k_min = a$cluster_k_min,
                               k_max = a$cluster_k_max, side = side),
        error = function(e) {
          if (identical(conditionMessage(e), "Not enough neighbors to cluster.")) list() else stop(e)
        })
      clusters <- c(clusters, cl)
    }
    ov <- .overall_interpretability(clusters, a$weight_by_size)
    c(k = ks[i], var_explained = rec$var, ov, beta_delta_1_minus_cos = rec$delta)
  })
  if (!is.null(bar)) close(bar)
  tab <- as.data.frame(do.call(rbind, rows))
  tab <- tab[order(tab$k), , drop = FALSE]

  det <- .detrend_by_variance(tab$var_explained, tab$aggregate)
  tab$interp_hat <- det$hat
  tab$interp_resid <- det$resid
  tab$interp_resid_z <- .zscore_ignore_nan(det$resid)
  tab$interp_auck <- .compute_auck(tab$interp_resid_z, a$auck_radius)
  tab$stab_good_raw <- -tab$beta_delta_1_minus_cos
  tab$stab_z_raw <- .zscore_ignore_nan(tab$stab_good_raw)
  tab$stab_auck_raw <- .compute_auck(tab$stab_z_raw, a$auck_radius)
  interp_any <- any(is.finite(tab$interp_auck))
  stab_any <- any(is.finite(tab$stab_auck_raw))
  tab$joint_score <- if (interp_any && stab_any) 0.5 * (tab$interp_auck + tab$stab_auck_raw)
                     else if (interp_any) tab$interp_auck else NA_real_

  finite <- is.finite(tab$joint_score)
  if (any(finite)) {
    best <- max(tab$joint_score[finite])
    best_k <- min(tab$k[finite][abs(tab$joint_score[finite] - best) <= 1e-12])
    if (verbose) message(sprintf("[sweep] best K = %d (joint score %.4f)", as.integer(best_k), best))
  } else {
    best_k <- min(tab$k)
    if (verbose) message(sprintf("[sweep] no finite joint score; using k_min = %d", as.integer(best_k)))
  }
  for (col in names(tab)) tab[[col]][is.nan(tab[[col]])] <- NA
  for (col in c("k", "n_clusters", "total_size")) tab[[col]] <- as.integer(tab[[col]])
  rownames(tab) <- NULL
  list(best_k = as.integer(best_k), table = tab[.SWEEP_COLUMNS])
}

# mean_coherence, mean_abs_cosb, aggregate, n_clusters, total_size of a list
# of clusters (each with size, coherence, centroid_cos_beta).
.overall_interpretability <- function(clusters, weight_by_size = TRUE) {
  if (!length(clusters)) {
    return(c(mean_coherence = NaN, mean_abs_cosb = NaN, aggregate = NaN, n_clusters = 0, total_size = 0))
  }
  sizes <- vapply(clusters, function(c) as.numeric(c$size), 0)
  coherence <- vapply(clusters, `[[`, 0, "coherence")
  abs_cosb <- abs(vapply(clusters, `[[`, 0, "centroid_cos_beta"))
  if (weight_by_size && sum(sizes) > 0) {
    w <- sizes / sum(sizes)
    coh <- sum(coherence * w)
    cosb <- sum(abs_cosb * w)
  } else {
    coh <- mean(coherence)
    cosb <- mean(abs_cosb)
  }
  c(mean_coherence = coh, mean_abs_cosb = cosb, aggregate = coh * cosb, n_clusters = length(clusters),
    total_size = sum(sizes))
}

# The fit and residuals of y ~ a + b * log(v) (v in percent), over the rows
# where both are finite and v > 0; NA with fewer than 3 such rows.
.detrend_by_variance <- function(v, y) {
  hat <- rep(NA_real_, length(y))
  ok <- is.finite(v) & is.finite(y) & v > 0
  if (sum(ok) < 3L) return(list(hat = hat, resid = hat))
  X <- cbind(1, log(v[ok]))
  coef <- qr.coef(qr(X), y[ok])
  hat[ok] <- drop(X %*% coef)
  list(hat = hat, resid = y - hat)
}

# z-scores over the finite values (population SD; an SD under 1e-12 counts as 1).
.zscore_ignore_nan <- function(x) {
  f <- x[is.finite(x)]
  if (!length(f)) return(rep(NA_real_, length(x)))
  m <- mean(f)
  s <- .sd_pop(f)
  if (!is.finite(s) || s < 1e-12) s <- 1
  (x - m) / s
}

# The mean of the finite values within `radius` positions of each position.
.compute_auck <- function(z, radius) {
  n <- length(z)
  vapply(seq_len(n), function(i) {
    w <- z[max(1L, i - radius):min(n, i + radius)]
    w <- w[is.finite(w)]
    if (length(w)) mean(w) else NA_real_
  }, 0)
}
