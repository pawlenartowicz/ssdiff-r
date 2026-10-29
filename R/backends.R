# The two Suggests backends (plskit, fastglmm): the check at call time and
# the installer. No Python counterpart.

# Stops unless `pkg` loads in at least the version DESCRIPTION's Suggests
# asks for. `what` names the caller, `why` (optional) says why it needs the
# package this time.
need_pkg <- function(pkg, what, why = NULL) {
  reason <- if (is.null(why)) "" else sprintf(" (for %s)", why)
  install <- sprintf(paste0(
    "  ssdiff::install_backends()\n",
    "or:\n",
    "  install.packages(\"%s\", repos = c(\"https://pawlenartowicz.r-universe.dev\", getOption(\"repos\")))"), pkg)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("%s needs the '%s' package%s, which is on r-universe, not CRAN. Install it with:\n%s",
                 what, pkg, reason, install), call. = FALSE)
  }
  min <- .min_version(pkg)
  if (!is.null(min) && utils::packageVersion(pkg) < min) {
    stop(sprintf("%s needs '%s' %s or later%s; version %s is installed. Update it with:\n%s",
                 what, pkg, min, reason, utils::packageVersion(pkg), install), call. = FALSE)
  }
}

# The minimum version of `pkg` in ssdiff's Suggests, or NULL.
.min_version <- function(pkg) {
  suggests <- suppressWarnings(utils::packageDescription("ssdiff", fields = "Suggests"))
  if (!is.character(suggests) || is.na(suggests)) return(NULL)
  m <- regmatches(suggests, regexec(sprintf("(^|,)\\s*%s\\s*\\(>=\\s*([^)[:space:]]+)\\s*\\)", pkg), suggests))[[1L]]
  if (length(m)) m[3L] else NULL
}

# Loads the package of a kept model before its methods run: in a new session
# a result read with readRDS() has a fastglmm fit whose S3 methods are not
# registered until fastglmm loads.
.model_pkg <- function(fit) {
  if (inherits(fit, "fastglmm")) need_pkg("fastglmm", "this fit_regression() result")
  invisible(fit)
}

#' Install the plskit and fastglmm backends
#'
#' Installs the packages that the PLS fits and the mixed-model path of
#' `fit_regression()` need. Both are on r-universe, not CRAN. Run it once.
#'
#' On Windows and macOS the r-universe repository serves binaries. On Linux
#' the prebuilt binary for Ubuntu `resolute` is tried first; a package whose
#' binary installs but does not load (for example on a system with glibc
#' older than 2.39) is then built from source, which needs Rust
#' (<https://rustup.rs>).
#'
#' @param pkgs Packages to install.
#' @param ... Passed to [utils::install.packages()]. When it holds `lib`, the
#'   load check on Linux looks there.
#' @return `NULL`, invisibly, on Linux; elsewhere the value of
#'   [utils::install.packages()].
#' @examples
#' \dontrun{
#' install_backends()
#' }
#' @export
install_backends <- function(pkgs = c("plskit", "fastglmm"), ...) {
  universe <- "https://pawlenartowicz.r-universe.dev"
  repos <- c(universe, getOption("repos"))
  if (Sys.info()[["sysname"]] != "Linux") {
    return(utils::install.packages(pkgs, repos = repos, ...))
  }
  r_minor <- sub("^(\\d+\\.\\d+).*", "\\1", getRversion())
  bin <- sprintf("%s/bin/linux/resolute-%s/%s", universe, R.version$arch, r_minor)
  utils::install.packages(pkgs, repos = c(bin, repos), ...)
  lib <- list(...)$lib
  broken <- pkgs[!vapply(pkgs, requireNamespace, logical(1),
                         lib.loc = lib, quietly = TRUE)]
  if (length(broken)) {
    message("The prebuilt Linux binaries of ", paste(broken, collapse = ", "),
            " do not load on this system. Building from source; this needs ",
            "Rust (https://rustup.rs).")
    utils::install.packages(broken, repos = repos, type = "source", ...)
  }
  invisible(NULL)
}
