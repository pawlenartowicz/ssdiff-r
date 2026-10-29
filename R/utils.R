# Small helpers shared by the fits and the results. No Python counterpart.

# Stops unless x is a whole number >= 1.
.check_count <- function(x, arg) {
  if (!.is_whole(x) || x < 1) {
    stop(sprintf("%s must be a whole number >= 1, got %s", arg, deparse1(x)), call. = FALSE)
  }
}

# Values quoted and joined for a message: 'a', 'b'.
.quoted <- function(x) paste0("'", x, "'", collapse = ", ")

# TRUE for one finite whole number.
.is_whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)

# An argument as a message shows it: 'words', or deparsed.
.shown <- function(x) if (is.character(x) && length(x) == 1L) .quoted(x) else deparse1(x)

# Population SD (ddof = 0).
.sd_pop <- function(v) sqrt(mean((v - mean(v))^2))

# A data frame from a list of rows (named lists), columns in `cols`.
.rows_df <- function(rows, cols) {
  out <- lapply(cols, function(c) unlist(lapply(rows, function(r) {
    v <- r[[c]]
    if (is.null(v)) NA else v
  })))
  names(out) <- cols
  if (!length(rows)) out <- lapply(stats::setNames(cols, cols), function(c) logical(0))
  as.data.frame(out, stringsAsFactors = FALSE)
}

# Cosine of two vectors; NaN when either is zero.
.cosine <- function(u, v) {
  nu <- sqrt(sum(u^2))
  nv <- sqrt(sum(v^2))
  if (nu < 1e-12 || nv < 1e-12) return(NaN)
  sum(u * v) / (nu * nv)
}
