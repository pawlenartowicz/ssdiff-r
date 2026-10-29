# Port of _RegressionModelView in ssdiff/results/continuous_result.py: R's
# model generics on a fit_regression() result. Change them together.

# Positions of the covariates (coefficients without a text component).
.covariate_names <- function(object) {
  b <- .reg_coef(object$model)
  names(b)[!.is_text_coef(names(b))]
}

#' Model generics of a regression result
#'
#' A [fit_regression()] result reads like the model it holds (`res$model`).
#' A text term has K component coefficients and no single estimate, so these
#' hold the covariates only: the model's coefficients without the `.ssd_*`
#' ones. The text terms are tested in `res$terms` and `res$test`, and the
#' whole model in `res$model_test`.
#'
#' `confint()` gives Wald intervals on every model path (t on the residual df
#' for `lm`, normal for `glm`, as fastglmm reports them), the intervals of
#' `res$covariates`. `df.residual()` is `nobs()` minus the fixed coefficients
#' that are not aliased, for every model (for a mixed model this is not a
#' Satterthwaite df). `logLik()` carries `df` and `nobs`, so [AIC()] and
#' [BIC()] work. `residuals()` (response residuals) and `fitted()` (fitted
#' means on the response scale) are named by the row names of `data` of the
#' rows used, or by the row numbers without `data`. `ranef()` gives one data
#' frame per grouping factor, one row per level and one column per
#' random-effect term; an empty list without random effects.
#'
#' @param object A result of [fit_regression()].
#' @param parm Covariates (names or positions among them); all by default.
#' @param level The confidence level.
#' @param ... Unused.
#' @return `coef()`: named estimates. `vcov()`: their covariance matrix.
#'   `confint()`: a matrix with the lower and upper limits. `nobs()`,
#'   `df.residual()`: numbers. `logLik()`: a `logLik` object. `residuals()`,
#'   `fitted()`: named vectors. `ranef()`: a named list of data frames.
#' @examples
#' \dontrun{
#' r <- fit_regression(s, score ~ age + sex * ssd(), data = df)
#' coef(r)
#' confint(r)
#' AIC(r)
#' }
#' @importFrom stats coef confint df.residual fitted logLik nobs residuals vcov
#' @name regression-generics
NULL

#' @rdname regression-generics
#' @export
coef.ssd_pca_ols <- function(object, ...) {
  b <- .reg_coef(object$model)
  b[.covariate_names(object)]
}

#' @rdname regression-generics
#' @export
vcov.ssd_pca_ols <- function(object, ...) {
  keep <- .covariate_names(object)
  .reg_vcov(object$model)[keep, keep, drop = FALSE]
}

#' @rdname regression-generics
#' @export
confint.ssd_pca_ols <- function(object, parm, level = 0.95, ...) {
  tab <- .reg_coef_table(object$model, level)
  tab <- tab[!.is_text_coef(tab$term), , drop = FALSE]
  if (!missing(parm)) tab <- tab[if (is.numeric(parm)) parm else match(parm, tab$term), , drop = FALSE]
  a <- (1 - level) / 2
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
  matrix(c(tab$ci_low, tab$ci_high), ncol = 2L, dimnames = list(tab$term, pct))
}

#' @rdname regression-generics
#' @export
nobs.ssd_pca_ols <- function(object, ...) length(object$y)

#' @rdname regression-generics
#' @export
df.residual.ssd_pca_ols <- function(object, ...) {
  length(object$y) - sum(!is.na(.reg_coef(object$model)))
}

#' @rdname regression-generics
#' @export
logLik.ssd_pca_ols <- function(object, ...) {
  .model_pkg(object$model)
  ll <- stats::logLik(object$model)
  attr(ll, "nobs") <- length(object$y)
  ll
}

#' @rdname regression-generics
#' @export
residuals.ssd_pca_ols <- function(object, ...) {
  stats::setNames(unname(.reg_resid(object$model)), object$obs_names)
}

#' @rdname regression-generics
#' @export
fitted.ssd_pca_ols <- function(object, ...) {
  stats::setNames(unname(.reg_fitted(object$model)), object$obs_names)
}

#' @rdname regression-generics
#' @importFrom nlme ranef
#' @export
ranef.ssd_pca_ols <- function(object, ...) {
  if (!inherits(object$model, "fastglmm")) return(list())
  .model_pkg(object$model)
  nlme::ranef(object$model)
}

#' @importFrom nlme ranef
#' @export
nlme::ranef
