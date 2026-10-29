# The results page (spec 04, "Page tests" 1 to 5): the format number, the
# shape of page_data() against Python's page fixtures, one-element arrays,
# the serialiser on Python's own output, and the file save_page() writes.

PAGE_KINDS <- c("pls", "pca_ols", "groups", "multipls", "regression", "contrast")

read_page_fixture <- function(kind) {
  jsonlite::read_json(test_path("fixtures", "page", paste0(kind, ".json")), simplifyVector = FALSE)
}

# page_data() as the page reads it: through the serialiser and parsed back.
parsed_page <- function(result, ...) jsonlite::fromJSON(page_json(page_data(result, ...)), simplifyVector = FALSE)

json_type <- function(v) {
  if (is.null(v)) return("null")
  if (is.list(v)) return(if (is.null(names(v))) "array" else "object")
  if (is.logical(v)) return("boolean")
  if (is.character(v)) return("string")
  "number"
}

# Differences in shape between two parsed pages, as messages: keys, table
# names, unit fields, columns in order, and per column the JSON types of the
# non-null cells, which must be among the fixture's.
page_shape_problems <- function(r, py) {
  out <- character(0)
  same <- function(a, b, where) {
    if (!identical(a, b)) out <<- c(out, sprintf("%s: %s, fixture %s", where, deparse1(a), deparse1(b)))
  }
  table <- function(rt, pt, where) {
    cols <- unlist(rt$columns)
    same(cols, unlist(pt$columns), paste(where, "columns"))
    if (!identical(cols, unlist(pt$columns))) return()
    for (j in seq_along(cols)) {
      types <- function(rows) setdiff(unique(vapply(rows, function(row) json_type(row[[j]]), "")), "null")
      bad <- setdiff(types(rt$rows), types(pt$rows))
      if (length(bad)) out <<- c(out, sprintf("%s column %s: %s cells", where, cols[[j]], paste(bad, collapse = ", ")))
    }
  }
  tables <- function(rt, pt, where) {
    same(names(rt), names(pt), paste(where, "table names"))
    for (name in intersect(names(rt), names(pt))) table(rt[[name]], pt[[name]], paste(where, name))
  }

  same(names(r), names(py), "top-level keys")
  same(names(r$software), names(py$software), "software keys")
  same(r$software$language, "r", "software.language")
  same(names(r$settings), names(py$settings), "settings keys")
  same(names(r$settings$clusters), names(py$settings$clusters), "settings.clusters keys")
  same(length(r$texts), length(py$texts), "number of texts")
  for (i in seq_len(min(length(r$texts), length(py$texts)))) {
    rt <- r$texts[[i]]
    pt <- py$texts[[i]]
    where <- sprintf("texts[%d]", i - 1L)
    same(names(rt), names(pt), paste(where, "keys"))
    same(names(rt$embeddings), names(pt$embeddings), paste(where, "embeddings keys"))
    same(names(rt$settings), names(pt$settings), paste(where, "settings keys"))
    same(is.null(rt$sweep), is.null(pt$sweep), paste(where, "has a sweep"))
    if (!is.null(rt$sweep) && !is.null(pt$sweep)) table(rt$sweep, pt$sweep, paste(where, "sweep"))
  }
  tables(r$tables, py$tables, "tables")
  same(length(r$units), length(py$units), "number of units")
  for (i in seq_len(min(length(r$units), length(py$units)))) {
    where <- sprintf("units[%d]", i - 1L)
    same(sort(names(r$units[[i]])), sort(names(py$units[[i]])), paste(where, "fields"))
    tables(r$units[[i]]$tables, py$units[[i]]$tables, where)
  }
  out
}

pls1_fit <- function() {
  cached_fit("pls1", function(r) fit_pls(r$s, "y", data = r$data, k = 2, n_splits = 5, random_state = 42))
}

# The regression.json fit: two texts, the second from token lists, with an interaction.
two_text_fit <- function() {
  cached_fit("two_texts", function(r) {
    diary <- ssd(rev(r$s$docs), r$emb, lexicon = r$d$lexicon)
    fit_regression(list(essay = r$s, diary = diary), y ~ ssd(essay) + group * ssd(diary), data = r$data,
                   k_min = 2, k_max = 6, k_step = 1)
  })
}

contrast_fit <- function() {
  contrast(pca_ols_fit()$gradient, groups_fit()$units[[1]]$gradient, load_wordpack(test_pack()),
           labels = c("pls", "pca_ols"))
}

test_that("each fixture's report_format is REPORT_FORMAT", {
  for (kind in PAGE_KINDS) expect_identical(read_page_fixture(kind)$report_format, REPORT_FORMAT, info = kind)
})

test_that("page_data() has the shape of Python's pages: pca_ols, groups, regression, contrast", with_empty_cache({
  fits <- list(pca_ols = pca_ols_fit, groups = groups_fit, regression = two_text_fit, contrast = contrast_fit)
  for (kind in names(fits)) {
    p <- parsed_page(fits[[kind]](), snippets_per_side = 10, snippets_per_cluster = 5)
    expect_identical(p$kind, if (kind == "regression") "pca_ols" else kind)
    expect_identical(page_shape_problems(p, read_page_fixture(kind)), character(0), info = kind)
  }
}))

test_that("page_data() has the shape of Python's pages: pls and multipls", with_empty_cache({
  skip_if_not_installed("plskit")
  for (kind in c("pls", "multipls")) {
    fit <- if (kind == "pls") pls1_fit() else multipls_fit()
    p <- parsed_page(fit, snippets_per_side = 10, snippets_per_cluster = 5)
    expect_identical(p$kind, kind)
    expect_identical(page_shape_problems(p, read_page_fixture(kind)), character(0), info = kind)
  }
}))

test_that("a token-list text has no snippet tables and its docs' text is the tokens", with_empty_cache({
  r <- two_text_fit()
  d <- page_data(r, snippets_per_side = 5, snippets_per_cluster = 2)
  expect_identical(vapply(d$texts, `[[`, TRUE, "snippets"), c(TRUE, FALSE))
  expect_identical(vapply(d$units, `[[`, 0L, "text"), c(0L, 1L, 1L, 1L))
  expect_true(all(c("snippets", "cluster_snippets") %in% names(d$units[[1]]$tables)))
  docs <- d$units[[2]]$tables$docs$rows
  expect_false(any(c("snippets", "cluster_snippets") %in% names(d$units[[2]]$tables)))
  has <- !is.na(docs$text)
  expect_true(any(has))
  expect_identical(docs$text[has], vapply(r$texts[[2]]$docs[docs$doc_id[has]], paste, "", collapse = " "))
}))

test_that("doc_texts keeps the text of both ends of each ranking", with_empty_cache({
  r <- pca_ols_fit()
  docs_of <- function(n) page_data(r, snippets_per_side = 1, snippets_per_cluster = 1, doc_texts = n)$units[[1]]$tables$docs$rows
  d <- docs_of(1)
  ends <- d$doc_id[c(which.min(d$alignment_score), which.max(d$alignment_score), which.min(d$residual),
                     which.max(d$residual))]
  expect_setequal(d$doc_id[!is.na(d$text)], ends)
  has <- !is.na(d$text)
  expect_identical(d$text[has], vapply(r$texts[[1]]$pre_docs[d$doc_id[has]], `[[`, "", "raw"))
  expect_false(anyNA(docs_of(Inf)$text))
  expect_true(all(is.na(docs_of(0)$text)))
}))

test_that("one-element arrays stay arrays", with_empty_cache({
  t <- results_text()
  one <- fit_regression(ssd(t$d$text, t$emb, lexicon = "people"), "y", data = t$data, fixed_k = 4)
  json <- page_json(page_data(one, snippets_per_side = 1, snippets_per_cluster = 1, doc_texts = 0))
  p <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  expect_identical(json_type(p$texts), "array")
  expect_length(p$texts, 1L)
  expect_identical(p$texts[[1]]$settings$lexicon, list("people"))
  expect_identical(json_type(p$units), "array")
  expect_identical(json_type(p$tables$stats$rows), "array")
  expect_length(p$tables$stats$rows, 1L)
  expect_identical(json_type(p$tables$stats$rows[[1]]), "array")
  expect_true(grepl("\"stats\":{\"columns\":[\"backend\",", json, fixed = TRUE))
  expect_true(grepl("\"rows\":[[\"PCA+OLS\",", json, fixed = TRUE))

  skip_if_not_installed("plskit")
  m <- fit_multipls(t$s, c("y", "y2"), data = t$data, ssd_gradients = 1, n_perm = 19, n_splits = 5, random_state = 0)
  p <- parsed_page(m, snippets_per_side = 1, snippets_per_cluster = 1, doc_texts = 0)
  expect_length(p$units, 1L)
  info <- p$tables$pls_info
  row <- info$rows[[1]]
  cols <- unlist(info$columns)
  expect_identical(json_type(row[[match("order", cols)]]), "array")
  expect_length(row[[match("order", cols)]], 1L)
  expect_identical(json_type(row[[match("signs", cols)]]), "array")
  expect_length(row[[match("signs", cols)]], 1L)
  expect_identical(json_type(p$tables$stats$rows[[1]]), "array")
}))

test_that("page_json() writes each fixture back to the same value", {
  for (kind in PAGE_KINDS) {
    py <- read_page_fixture(kind)
    expect_equal(jsonlite::fromJSON(page_json(py), simplifyVector = FALSE), py, tolerance = 1e-12, info = kind)
  }
})

test_that("save_page() writes the template around the escaped data", with_empty_cache({
  cc <- contrast_fit()
  tpl_path <- system.file("page_template.html", package = "ssdiff")
  tpl <- readBin(tpl_path, "raw", file.size(tpl_path))
  at <- regexpr("__SSD_DATA__", rawToChar(tpl), fixed = TRUE, useBytes = TRUE)
  before <- tpl[seq_len(at - 1L)]
  after <- tpl[(at + attr(at, "match.length")):length(tpl)]
  path <- tempfile(fileext = ".html")
  on.exit(unlink(path))

  # The bytes around the data are the template's; the data is the page's JSON.
  read_back <- function(path) {
    bytes <- readBin(path, "raw", file.size(path))
    expect_identical(bytes[seq_along(before)], before)
    expect_identical(utils::tail(bytes, length(after)), after)
    html <- rawToChar(bytes)
    Encoding(html) <- "UTF-8"
    expect_false(grepl("__SSD_DATA__", html, fixed = TRUE))
    expect_identical(lengths(regmatches(html, gregexpr("id=\"ssd-data\"", html, fixed = TRUE))), 1L)
    json <- rawToChar(bytes[(length(before) + 1L):(length(bytes) - length(after))])
    Encoding(json) <- "UTF-8"
    list(html = html, data = jsonlite::fromJSON(json, simplifyVector = FALSE))
  }

  title <- "a </script><script>alert(1)</script> <!--<script> b"
  expect_invisible(out <- save_page(cc, path, title = title))
  expect_identical(out, path)
  got <- read_back(path)
  expect_false(grepl("</script><script>alert", got$html, fixed = TRUE))
  expect_false(grepl("<!--<script", got$html, fixed = TRUE))
  expect_identical(got$data$title, title)
  want <- jsonlite::fromJSON(page_json(page_data(cc, title = title)), simplifyVector = FALSE)
  got$data$created <- NULL
  want$created <- NULL
  expect_equal(got$data, want)

  polish <- "Za\u017c\u00f3\u0142\u0107 g\u0119\u015bl\u0105 ja\u017a\u0144"
  save_page(cc, path, title = polish, dataset = c(File = "C:\\data\\essays.csv", Note = "\\1 $& \\\\"))
  got <- read_back(path)$data
  expect_identical(got$title, polish)
  expect_identical(got$dataset, list(File = "C:\\data\\essays.csv", Note = "\\1 $& \\\\"))
  expect_identical(got$texts[[1]]$snippets, NULL)
}))

test_that("page_data() checks its arguments", with_empty_cache({
  cc <- contrast_fit()
  expect_error(page_data(list()), "no results page for list", fixed = TRUE)
  expect_error(page_data(cc, title = 1), "title must be a string", fixed = TRUE)
  expect_error(page_data(cc, dataset = c("a", "b")), "dataset must be a named list", fixed = TRUE)
  expect_error(page_data(cc, doc_texts = -1), "doc_texts must be a whole number", fixed = TRUE)
  expect_error(page_data(cc, snippets_per_side = 1.5), "snippets_per_side must be a whole number", fixed = TRUE)
  expect_error(page_data(cc, snippets_per_cluster = "5"), "snippets_per_cluster must be a whole number", fixed = TRUE)
  d <- page_data(cc, title = "", dataset = list())
  expect_identical(d$title, "SSD results")
  expect_identical(as.character(page_json(d["dataset"])), "{\"dataset\":{}}")
  expect_null(page_data(cc)$dataset)
}))
