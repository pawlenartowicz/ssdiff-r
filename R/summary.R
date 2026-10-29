# Port of the reprs, summary() and predict() of ssdiff/results/*
# (continuous_result.py, group_result.py, multi_pls_result.py,
# contrast_result.py). Change them together. R prints its own tables, so the
# lm/glm and fastglmm paths print the same thing.

# A bounded number as Python's fmt_r writes it: ".42", "-.18", "1.00".
.fmt_r <- function(r, digits = 2L) {
  if (is.null(r) || length(r) != 1L || !is.finite(r)) return("nan")
  s <- sprintf("%.*f", digits, abs(r))
  if (abs(r) < 1) s <- sub("^0", "", s)
  if (r < 0) paste0("\u2212", s) else s
}

.fmt_count <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

# Items joined by two spaces, wrapped at the console width after `label`.
.print_items <- function(label, items) {
  width <- max(40L, getOption("width", 80L))
  indent <- strrep(" ", nchar(label) + 2L)
  lines <- character(0)
  line <- paste0("  ", label)
  first <- TRUE
  for (it in items) {
    if (!first && nchar(line) + 2L + nchar(it) > width) {
      lines <- c(lines, line)
      line <- paste0(indent, it)
    } else {
      line <- if (first) paste0(line, it) else paste0(line, "  ", it)
    }
    first <- FALSE
  }
  cat(c(lines, line), sep = "\n")
}

# The repr's reminder of what to call next.
.print_access <- function(accessors, fields) {
  .print_items("accessors: ", accessors)
  .print_items("fields:    ", paste0("res$", fields))
}

.DOC_ACCESSORS <- c("words(res)", "clusters(res)", "cluster_words(res)", "snippets(res)", "cluster_snippets(res)",
                    "docs(res)", "doc_text(res, doc_id)")

.dims_line <- function(result) {
  if (!is.null(result$units)) .print_items("units:     ", paste0("\"", names(result$units), "\""))
}

#' Print and summarise a result
#'
#' `print()` shows a result's headline (r², p-value, n, backend; for groups
#' the omnibus test) and lists the accessors and fields to read next.
#' `summary()` prints a longer account and returns its tables invisibly:
#'
#' - [fit_pls()]: the number of components, the split-half test and the R²
#'   (the partial R² with a formula); returns `stats` and `test`.
#' - [fit_regression()]: the formula, the covariates' coefficients
#'   (estimate, standard error, t or z, p), one Wald test per text term,
#'   with several text terms the joint tests (all text terms; each text)
#'   with the R² they add (OLS only), the fit statistics, and the test of the
#'   whole model against the intercept-only model: the F test for OLS, as
#'   [lm()] prints it, else the likelihood-ratio test. It returns a list
#'   with `coefficients` (the covariates' matrix, as [coef()] and [vcov()]
#'   give them), `terms` and `model_test`. `summary(res$model)` gives the
#'   model package's own summary.
#' - [fit_groups()]: the omnibus test and the pairs; returns `test` and
#'   `pairs`.
#' - [fit_multipls()]: the stats and each dimension's `cov_share` and
#'   `score_corr`; returns `stats` and `dims`.
#' - [contrast()]: the two labels, `cos_g1_g2` and `diff_norm`; returns
#'   `stats`.
#'
#' The tables of [snippets()] and [cluster_snippets()] print with
#' `text_window` and `text_surface` clipped to 40 characters; the data frames
#' keep the full text.
#'
#' @param x,object A result, or for `print()` a snippet table.
#' @param ... Unused.
#' @return `print()`: `x`, invisibly. `summary()`: a list, invisibly.
#' @examples
#' \dontrun{
#' r <- fit_regression(s, score ~ age + ssd(), data = df)
#' r
#' summary(r)
#' }
#' @name result-print
NULL

#' @rdname result-print
#' @export
print.ssd_pls <- function(x, ...) {
  s <- x$stats
  cat(sprintf("ssd_pls  r\u00b2=%s  p=%s  n=%s  backend=%s\n", .fmt_r(s$r2), .fmt_p(s$pvalue),
              .fmt_count(s$n_kept), s$backend))
  .dims_line(x)
  .print_access(c(.DOC_ACCESSORS, "predict(res, newdata)", "summary(res)"),
                c("stats", "fit_info", "test", "beta", "gradient", "y", "y_hat", "residuals", "alignment_scores",
                  "component_scores", "component_weights", if (!is.null(x$units)) c("pls_info", "units")))
  invisible(x)
}

#' @rdname result-print
#' @export
print.ssd_pca_ols <- function(x, ...) {
  s <- x$stats
  r2 <- if (is.na(s$r2)) "" else sprintf("r\u00b2=%s  ", .fmt_r(s$r2))
  cat(sprintf("ssd_pca_ols  %sp=%s  n=%s  backend=%s\n", r2, .fmt_p(s$pvalue), .fmt_count(s$n_kept), s$backend))
  cat("  formula:   ", .formula_shown(x), "\n", sep = "")
  .dims_line(x)
  .print_access(c(.DOC_ACCESSORS, "sweep_table(res)", "plot_sweep(res)", "coef(res)", "confint(res)",
                  "predict(res, newdata)", "summary(res)"),
                c("stats", "fit_info", "test", "terms", "covariates", "model_test", "model",
                  if (is.null(x$units) || !is.null(x$pls_info)) c("beta", "gradient", "alignment_scores"),
                  "y", "y_hat", if (!is.null(x$units)) c(if (!is.null(x$pls_info)) "pls_info", "units")))
  invisible(x)
}

#' @rdname result-print
#' @export
print.ssd_groups <- function(x, ...) {
  cat(sprintf("ssd_groups  G=%d  n=%s  omnibus p=%s\n", x$G, .fmt_count(x$n_kept), .fmt_p(x$test$omnibus_p)))
  lab <- x$group_labels
  .print_items("pairs:     ", sprintf("%s (%s vs %s)", x$pairs$contrast, lab[x$pairs$g1], lab[x$pairs$g2]))
  .print_access(c(.DOC_ACCESSORS, "summary(res)"), c("stats", "test", "pairs", "group_labels", "groups", "units"))
  invisible(x)
}

#' @rdname result-print
#' @export
print.ssd_multipls <- function(x, ...) {
  s <- x$stats
  cat(sprintf("ssd_multipls  ssd_gradients=%d  q=%d  rotate=%s  cov_explained=%s  p=%s  n=%s\n", s$ssd_gradients,
              s$q, if (is.na(s$rotate)) "NULL" else s$rotate, .fmt_r(s$cov_explained), .fmt_p(s$pvalue),
              .fmt_count(s$n)))
  .dims_line(x)
  .print_access(c(.DOC_ACCESSORS, "outcomes(res)", "summary(res)"),
                c("stats", "pls_info", "test", "outcome_names", "outcome_groups", "U", "V", "U_rot", "V_rot", "R",
                  "singular_values", "x_scores", "y_scores", "units"))
  invisible(x)
}

#' @rdname result-print
#' @export
print.ssd_contrast <- function(x, ...) {
  cat(sprintf("ssd_contrast  %s \u2212 %s  cos(g1, g2)=%s  diff_norm=%.3f\n", x$labels[1L], x$labels[2L],
              .fmt_r(x$stats$cos_g1_g2), x$stats$diff_norm))
  .print_access(c("words(res)", "clusters(res)", "cluster_words(res)", "summary(res)"),
                c("stats", "labels", "beta", "gradient"))
  invisible(x)
}

# ---------------------------------------------------------------------------
# summary()
# ---------------------------------------------------------------------------

#' @rdname result-print
#' @export
summary.ssd_pls <- function(object, ...) {
  fi <- object$fit_info
  t <- object$test
  formula <- !is.na(fi$formula)
  head <- paste0("PLS, k = ", fi$n_components, if (formula) " (on the residuals of the formula)")
  test <- if (is.na(t$pvalue)) "test not run" else {
    at_k <- if (!is.na(fi$p_at_k) && fi$p_at_k != fi$n_components) paste0(", k = ", fi$p_at_k) else ""
    sprintf("split-half test (%s%s): p = %s", t$name, at_k, .fmt_p(t$pvalue))
  }
  r2 <- sprintf("%s = %s", if (formula) "partial R\u00b2" else "R\u00b2", .fmt_r(object$stats$r2, 3L))
  cat(sprintf("%s fit [ssdiff]\nNumber of obs: %d\n\nText:  %s\n       %s\n       %s\n", object$stats$backend,
              as.integer(object$stats$n_kept), head, test, r2))
  invisible(list(stats = object$stats, test = object$test))
}

#' @rdname result-print
#' @export
summary.ssd_groups <- function(object, ...) {
  t <- object$test
  cat(sprintf("Group comparison [ssdiff]: G = %d, n = %s\n", object$G, .fmt_count(object$n_kept)))
  cat(sprintf("Omnibus permutation test: T = %.4g, p = %s (%s permutations)\n", t$omnibus_T, .fmt_p(t$omnibus_p),
              .fmt_count(t$n_perm)))
  cat("Groups: ", paste(sprintf("%s = %s", names(object$group_labels), object$group_labels), collapse = ", "),
      "\n\n", sep = "")
  cat(sprintf("Pairs (p corrected: %s):\n", t$correction))
  print(object$pairs, row.names = FALSE, digits = 4)
  invisible(list(test = object$test, pairs = object$pairs))
}

#' @rdname result-print
#' @export
summary.ssd_multipls <- function(object, ...) {
  s <- object$stats
  cat(sprintf("PLS3 fit [ssdiff]: q = %d outcomes, n = %s, ssd_gradients = %d, rotate = %s\n", s$q,
              .fmt_count(s$n), s$ssd_gradients, if (is.na(s$rotate)) "NULL" else s$rotate))
  cat(sprintf("cov_explained = %s, LV1 test (%s): p = %s\n\n", .fmt_r(s$cov_explained, 3L), object$test$name,
              .fmt_p(s$pvalue)))
  dims <- data.frame(dim = names(object$units),
                     cov_share = vapply(object$units, function(u) u$stats$cov_share, 0, USE.NAMES = FALSE),
                     score_corr = vapply(object$units, function(u) u$stats$score_corr, 0, USE.NAMES = FALSE),
                     stringsAsFactors = FALSE)
  cat("Dimensions:\n")
  print(dims, row.names = FALSE, digits = 4)
  invisible(list(stats = object$stats, dims = dims))
}

#' @rdname result-print
#' @export
summary.ssd_contrast <- function(object, ...) {
  cat(sprintf("Contrast [ssdiff]: %s \u2212 %s\ncos(g1, g2) = %s   diff_norm = %.3f\n", object$labels[1L],
              object$labels[2L], .fmt_r(object$stats$cos_g1_g2), object$stats$diff_norm))
  invisible(list(stats = object$stats))
}

# The formula as Python's summary shows it: as written, with `+ ssd()` when it
# has no ssd() term, or `<outcome> ~ ssd()` without a formula.
.formula_shown <- function(result) {
  if (is.null(result$formula)) return(paste(result$outcome_name, "~ ssd()"))
  deparse1(.with_text_term(result$formula))
}

# The covariates' lm-style matrix: estimate, standard error, t or z, p.
.coef_matrix <- function(object) {
  m <- object$model
  tab <- .reg_coef_table(m)
  tab <- tab[!.is_text_coef(tab$term), , drop = FALSE]
  stat <- .stat_letter(m)
  out <- cbind(tab$estimate, tab$se, tab$statistic, tab$p_value)
  dimnames(out) <- list(tab$term, c("Estimate", "Std. Error", paste(stat, "value"), sprintf("Pr(>|%s|)", stat)))
  out
}

# "t" or "z": the statistic the model's coefficient tests use.
.stat_letter <- function(m) {
  .model_pkg(m)
  if (inherits(m, "fastglmm")) {
    return(if ("t value" %in% colnames(summary(m)$coefficients)) "t" else "z")
  }
  if (inherits(m, "glm")) "z" else "t"
}

#' @rdname result-print
#' @export
summary.ssd_pca_ols <- function(object, ...) {
  m <- object$model
  digits <- max(3L, getOption("digits") - 3L)
  coefs <- .coef_matrix(object)
  terms <- object$terms
  fam <- object$fit_args
  cat("Formula: ", .formula_shown(object), "\n", sep = "")
  cat(sprintf("Family: %s (%s link)   backend: %s\n", fam$family, fam$link, object$stats$backend))

  cat("\nCovariates:\n")
  if (nrow(coefs)) {
    stats::printCoefmat(coefs, digits = digits, na.print = "NA")
  } else {
    cat("  (none)\n")
  }

  wald_f <- identical(terms$name[1L], "wald_f")
  cat("\nText terms (Wald tests):\n")
  tt <- cbind(Df = terms$df, terms$statistic, terms$pvalue)
  dimnames(tt) <- list(terms$term, c("Df", if (wald_f) "F value" else "Chisq", if (wald_f) "Pr(>F)" else "Pr(>Chisq)"))
  stats::printCoefmat(tt, digits = digits, cs.ind = NULL, tst.ind = 2L, has.Pvalue = TRUE, P.values = TRUE,
                      na.print = "NA")

  if (nrow(terms) > 1L) {
    joint <- object$test[object$test$name %in% c("f_test", "lrt"), , drop = FALSE]
    has_r2 <- any(is.finite(joint$r2_delta))
    f_test <- identical(joint$name[1L], "f_test")
    jt <- cbind(joint$r2_delta, joint$statistic, joint$pvalue)
    dimnames(jt) <- list(ifelse(is.na(joint$label), "all text terms", joint$label),
                         c("\u0394R\u00b2", if (f_test) "F value" else "Chisq", if (f_test) "Pr(>F)" else "Pr(>Chisq)"))
    if (!has_r2) jt <- jt[, -1L, drop = FALSE]
    cat("\nJoint tests:\n")
    stats::printCoefmat(jt, digits = digits, cs.ind = NULL, tst.ind = ncol(jt) - 1L, has.Pvalue = TRUE,
                        P.values = TRUE, na.print = "NA")
  }

  n <- stats::nobs(object)
  dfr <- stats::df.residual(object)
  intercept <- "(Intercept)" %in% names(.reg_coef(m))
  df_model <- n - dfr - intercept
  ks <- vapply(object$texts, `[[`, 0L, "k")
  names_k <- vapply(object$texts, function(t) if (is.null(t$name)) "" else sprintf(" (%s)", t$name), "")
  ll <- stats::logLik(object)
  cat(sprintf("\nFit: n = %d, K = %s", n, paste0(ks, names_k, collapse = ", ")))
  if (!is.na(object$stats$r2)) {
    cat(sprintf(", R\u00b2 = %s, adjusted R\u00b2 = %s", formatC(object$stats$r2, digits = digits),
                formatC(object$stats$r2_adj, digits = digits)))
  }
  cat(sprintf("\n     logLik = %s (df = %s), AIC = %s, BIC = %s\n", format(signif(as.numeric(ll), digits + 2L)),
              format(attr(ll, "df")), format(signif(stats::AIC(object), digits + 2L)),
              format(signif(stats::BIC(object), digits + 2L))))

  mt <- object$model_test
  if (!is.na(mt$fvalue)) {
    cat("F-statistic:", formatC(mt$fvalue, digits = digits), "on", df_model, "and", dfr, "DF,  p-value:",
        format.pval(mt$f_pvalue, digits = digits), "\n")
  } else if (!is.na(mt$llr)) {
    cat("Likelihood-ratio test against the intercept-only model:", formatC(mt$llr, digits = digits), "on", df_model,
        "DF,  p-value:", format.pval(mt$llr_pvalue, digits = digits), "\n")
  } else {
    cat("No test of the whole model: the formula has no intercept.\n")
  }
  invisible(list(coefficients = coefs, terms = terms, model_test = mt))
}

# ---------------------------------------------------------------------------
# predict()
# ---------------------------------------------------------------------------

#' Predict the outcome of new texts
#'
#' Builds the new documents' vectors as the fit built its own: the same pack,
#' lexicon, window, `sif_a` and mode, with SIF weights from the training
#' documents' word counts. The prediction is the text's part on the
#' outcome's scale, `y_mean + y_std * (x - x_mean) %*% beta`, so on the
#' training texts it equals `res$y_hat`.
#'
#' Only fits without a formula, and for [fit_regression()] only the
#' Gaussian family with the identity link: another model's prediction needs
#' the new documents' covariates and random effects, or the model's
#' intercept and inverse link.
#'
#' @param object A result of [fit_pls()] or [fit_regression()].
#' @param newdata Texts in the form of the fit's text: a character vector
#'   (tokenized with the result's pack) for a fit on strings, a list of
#'   token lists for a fit on token lists.
#' @param ... Unused.
#' @return A numeric vector with one prediction per text; `NA` for a text
#'   that gets no vector (no word in the pack, or no lexicon word).
#' @examples
#' \dontrun{
#' r <- fit_pls(s, "score", data = df)
#' predict(r, c("A new essay about work.", "Another one."))
#' }
#' @name predict.ssd_result
NULL

#' @rdname predict.ssd_result
#' @importFrom stats predict
#' @export
predict.ssd_pls <- function(object, newdata, ...) .predict_texts(object, newdata)

#' @rdname predict.ssd_result
#' @export
predict.ssd_pca_ols <- function(object, newdata, ...) .predict_texts(object, newdata)

.predict_texts <- function(object, newdata) {
  if (!is.null(object$fit_args$formula) || identical(object$stats$backend, "PCA+GLMM")) {
    stop("predict() works only for fits without a formula, and for fit_regression() only with the Gaussian family.",
         call. = FALSE)
  }
  if (missing(newdata)) stop("predict() needs newdata: the texts to predict.", call. = FALSE)
  t <- object$texts[[1L]]
  token_lists <- is.list(newdata) && !is.object(newdata)
  if (is.null(t$pre_docs) && !token_lists) {
    stop("This result was fit on token lists; pass token lists, one list of tokens per document.", call. = FALSE)
  }
  if (!is.null(t$pre_docs) && !is.character(newdata)) {
    stop("This result was fit on strings; pass strings, not token lists.", call. = FALSE)
  }
  emb <- .embeddings_from_ref(t$pack)
  tok <- .tokenize_texts(newdata, emb)
  v <- .doc_vectors(tok$docs, emb, t$lexicon, t$window, t$sif_a, t$use_full_doc,
                    if (!is.null(t$sif_a)) .sif_counts(t$docs))
  y_st <- .standardize(object$y)
  out <- rep(NA_real_, length(tok$docs))
  if (nrow(v$x)) out[v$keep] <- y_st$mean + y_st$scale * drop(sweep(v$x, 2L, colMeans(t$x)) %*% object$beta)
  out
}
