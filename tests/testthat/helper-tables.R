# Python's {"columns": [...], "rows": [[...], ...]} tables as data frames;
# null becomes NA.
py_table <- function(tab) {
  cols <- unlist(tab$columns)
  out <- lapply(seq_along(cols), function(j) {
    v <- lapply(tab$rows, function(r) if (is.null(r[[j]])) NA else r[[j]])
    if (all(vapply(v, function(x) length(x) == 1L && !is.list(x), TRUE))) unlist(v) else v
  })
  names(out) <- cols
  as.data.frame(lapply(out, function(v) if (is.list(v)) I(v) else v), stringsAsFactors = FALSE)
}

golden_texts <- function() {
  jsonlite::read_json(test_path("fixtures", "golden", "texts.json"), simplifyVector = TRUE)
}

parity_case <- function(name) jsonlite::read_json(file.path(case_dir(), paste0(name, ".json")))
