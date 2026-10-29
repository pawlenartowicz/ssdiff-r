# Fits of the golden text for the result tests. Call inside with_empty_cache().

# The golden text (240 documents, lexicon) and its data columns, the pack loaded.
results_text <- function() {
  d <- golden_texts()
  emb <- load_wordpack(test_pack())
  list(d = d, emb = emb, s = ssd(d$text, emb, lexicon = d$lexicon),
       data = data.frame(y = d$y, y2 = d$y2, cov = d$cov, group = d$group, stringsAsFactors = FALSE))
}

.results_fits <- new.env()

# A fit made once per test run by `make(results_text())`; the pack is loaded
# again each call, as with_empty_cache() empties the registry.
cached_fit <- function(name, make) {
  if (is.null(.results_fits[[name]])) {
    .results_fits[[name]] <- make(results_text())
  } else {
    load_wordpack(test_pack())
  }
  .results_fits[[name]]
}

groups_fit <- function() {
  cached_fit("groups", function(r) fit_groups(r$s, "group", data = r$data, n_perm = 50, random_state = 42))
}

# The page fixtures' pca_ols fit: the sweep over K = 2..6.
pca_ols_fit <- function() {
  cached_fit("pca_ols", function(r) fit_regression(r$s, "y", data = r$data, k_min = 2, k_max = 6, k_step = 1))
}

interaction_fit <- function() {
  cached_fit("interaction", function(r) fit_regression(r$s, y ~ cov + group * ssd(), data = r$data, fixed_k = 10))
}

pls_fit <- function() {
  cached_fit("pls", function(r) fit_pls(r$s, "y", data = r$data, k = 2, ssd_gradients = 2, test = FALSE))
}

multipls_fit <- function() {
  cached_fit("multipls", function(r) {
    fit_multipls(r$s, c("y", "y2"), data = r$data, ssd_gradients = 2, n_perm = 19, n_splits = 5, random_state = 0)
  })
}

# A page fixture, and one table of its first unit (or of unit `i`) as a data frame.
page_fixture <- function(kind) jsonlite::read_json(test_path("fixtures", "page", paste0(kind, ".json")))
page_columns <- function(page, table, i = 1L) unlist(page$units[[i]]$tables[[table]]$columns)
page_table <- function(page, table, i = 1L) py_table(page$units[[i]]$tables[[table]])

# A result holding only a gradient in a pack, for the accessors on Python's betas.
gradient_result <- function(beta, emb) {
  structure(list(beta = beta, gradient = .unit_vector(beta), texts = list(list(pack = .pack_ref(emb)))),
            class = c("ssd_contrast", "ssd_result"))
}

# Runs the R lines `code` in a new R session with ssdiff loaded (the source
# tree under devtools, else the installed package); returns its output.
in_fresh_session <- function(code) {
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  load <- if (file.exists(file.path(root, "DESCRIPTION")) && "pkgload" %in% loadedNamespaces()) {
    sprintf("pkgload::load_all(%s, quiet = TRUE, helpers = FALSE)", deparse(root))
  } else {
    "library(ssdiff)"
  }
  script <- tempfile(fileext = ".R")
  writeLines(c(load, code), script)
  system2(file.path(R.home("bin"), "Rscript"), shQuote(script), stdout = TRUE, stderr = TRUE)
}
