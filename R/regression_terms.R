# Port of ssdiff/backends/regression.py (the ssd() terms of a fit_regression
# formula). Change the two together. R reads the terms from the formula's call
# tree and its null formulas from stats::terms(), not from strings.

# TRUE for a call to ssd().
.is_ssd_call <- function(e) is.call(e) && identical(e[[1L]], as.name("ssd"))

.is_empty_arg <- function(e) identical(e, quote(expr = ))

# Names inside the ssd(...) terms of an expression, in order of first
# appearance; NA for a bare ssd().
.text_names <- function(e) {
  out <- character(0)
  walk <- function(e) {
    if (.is_ssd_call(e)) {
      if (length(e) == 1L) {
        name <- NA_character_
      } else if (length(e) == 2L && is.name(e[[2L]])) {
        name <- as.character(e[[2L]])
      } else {
        stop(sprintf("A text term is ssd() or ssd(name) with one bare column name, got %s.", deparse1(e)),
             call. = FALSE)
      }
      if (!name %in% out) out <<- c(out, name)
    } else if (is.call(e)) {
      for (a in as.list(e)[-1L]) if (!.is_empty_arg(a)) walk(a)
    }
  }
  walk(e)
  out
}

# The formula with ` + ssd()` appended when it has no ssd(...) term.
.with_text_term <- function(f) {
  if (length(.text_names(f))) return(f)
  f[[3L]] <- call("+", f[[3L]], quote(ssd()))
  f
}

.is_bar_call <- function(e) is.call(e) && (identical(e[[1L]], as.name("|")) || identical(e[[1L]], as.name("||")))

# A random-effect group: `( ... | ... )`.
.is_random_term <- function(e) is.call(e) && identical(e[[1L]], as.name("(")) && .is_bar_call(e[[2L]])

.is_offset_term <- function(e) is.call(e) && identical(e[[1L]], as.name("offset"))

# TRUE when the expression has a random-effect bar anywhere.
.has_bar <- function(e) {
  if (.is_bar_call(e)) return(TRUE)
  is.call(e) && any(vapply(as.list(e)[-1L], function(a) !.is_empty_arg(a) && .has_bar(a), TRUE))
}

# Stops for an ssd(...) inside a random-effect term.
.check_text_terms <- function(f) {
  walk <- function(e, shown) {
    if (!is.call(e)) return(invisible())
    if (.is_bar_call(e)) {
      if (length(.text_names(e))) {
        stop(sprintf("ssd(...) cannot be used inside a random-effect term: %s", deparse1(shown)), call. = FALSE)
      }
      return(invisible())
    }
    for (a in as.list(e)[-1L]) if (!.is_empty_arg(a)) walk(a, if (.is_random_term(e)) e else a)
  }
  walk(f[[length(f)]], f[[length(f)]])
}

# The placeholder symbol name of a text term: `.ssd_<name>`, `.ssd_text` for ssd().
.placeholder <- function(name) if (is.na(name)) ".ssd_text" else paste0(".ssd_", name)

# Every ssd(...) replaced by its placeholder symbol.
.to_placeholders <- function(e) {
  if (.is_ssd_call(e)) return(as.name(.placeholder(.text_names(e))))
  if (!is.call(e)) return(e)
  for (i in seq_along(e)[-1L]) if (!.is_empty_arg(e[[i]])) e[[i]] <- .to_placeholders(e[[i]])
  e
}

# Symbols renamed or replaced: `map` is a named list from symbol name to the
# expression that replaces it.
.replace_symbols <- function(e, map) {
  if (is.name(e)) {
    hit <- map[[as.character(e)]]
    return(if (is.null(hit)) e else hit)
  }
  if (!is.call(e)) return(e)
  for (i in seq_along(e)[-1L]) if (!.is_empty_arg(e[[i]])) e[[i]] <- .replace_symbols(e[[i]], map)
  e
}

# Each placeholder replaced by `(col1 + col2 + ...)`; `columns` maps a
# placeholder to its component column names.
.expand <- function(f, columns) {
  map <- lapply(columns, function(cols) {
    call("(", Reduce(function(a, b) call("+", a, b), lapply(cols, as.name)))
  })
  .replace_symbols(f, map)
}

# An expression split into its random-effect and offset terms (`extras`, as
# written, in order) and the rest (`fixed`, NULL when nothing is left). The
# terms are found anywhere in the `+` / `-` tree, so `- 1` after them keeps
# them, as Python's scan of the whole right side does.
.split_extras <- function(e) {
  if (.is_random_term(e) || .is_offset_term(e)) return(list(fixed = NULL, extras = list(e)))
  # A parenthesised group of fixed terms is looked into, `( ... | ... )` aside.
  if (is.call(e) && identical(e[[1L]], as.name("(")) && length(e) == 2L) {
    inner <- .split_extras(e[[2L]])
    return(list(fixed = if (!is.null(inner$fixed)) call("(", inner$fixed), extras = inner$extras))
  }
  if (!is.call(e) || !(identical(e[[1L]], as.name("+")) || identical(e[[1L]], as.name("-")))) {
    return(list(fixed = e, extras = list()))
  }
  parts <- lapply(as.list(e)[-1L], .split_extras)
  extras <- unlist(lapply(parts, `[[`, "extras"), recursive = FALSE)
  fixed <- lapply(parts, `[[`, "fixed")
  op <- e[[1L]]
  if (length(fixed) == 1L) return(list(fixed = if (!is.null(fixed[[1L]])) call(as.character(op), fixed[[1L]]),
                                       extras = extras))
  a <- fixed[[1L]]
  b <- fixed[[2L]]
  out <- if (is.null(a) && is.null(b)) NULL
         else if (is.null(b)) a
         else if (is.null(a)) (if (identical(op, as.name("+"))) b else call("-", b))
         else call(as.character(op), a, b)
  list(fixed = out, extras = extras)
}

# The fixed part of formula `f` (its term labels and intercept) and the
# random-effect and offset terms of its right side, as written.
.formula_terms <- function(f) {
  sp <- .split_extras(f[[length(f)]])
  fixed <- if (is.null(sp$fixed)) 1 else sp$fixed
  tt <- stats::terms(stats::as.formula(call("~", quote(.ssd_y), fixed), env = baseenv()))
  list(labels = attr(tt, "term.labels"), intercept = attr(tt, "intercept") == 1L, extras = sp$extras)
}

# The formula without the fixed terms in `drop` (labels of .formula_terms).
# The response, a removed intercept, random-effect and offset terms are kept.
.null_formula <- function(f, drop) {
  ft <- .formula_terms(f)
  keep <- setdiff(ft$labels, drop)
  if (!length(keep) && !ft$intercept) {
    stop(sprintf("'%s' has no fixed term left once the text terms are dropped and the intercept is removed",
                 deparse1(f)), call. = FALSE)
  }
  rhs <- c(if (!ft$intercept) "0", if (length(keep)) keep else "1", vapply(ft$extras, deparse1, ""))
  stats::as.formula(paste(deparse1(f[[2L]]), "~", paste(rhs, collapse = " + ")), env = environment(f))
}

# Labels whose `:`-separated variables include the placeholder `ph`.
.text_labels <- function(labels, ph) labels[vapply(strsplit(labels, ":", fixed = TRUE), function(p) ph %in% p, TRUE)]

# A label with its placeholder written as the text term (`sex:ssd(essay)`).
.label_shown <- function(label, ph, shown) {
  parts <- strsplit(label, ":", fixed = TRUE)[[1L]]
  paste(ifelse(parts == ph, shown, parts), collapse = ":")
}

# Coefficient names grouped by text term. `components` maps a component column
# to list(shown text term, component index). A name with a component column
# among its `:` parts gets as key the name with that part replaced by the
# shown term. Returns a named list (by key, in order of first appearance) of
# data frames with `idx` (component index) and `pos` (position in names),
# sorted by idx.
.term_groups <- function(names, components) {
  groups <- list()
  for (pos in seq_along(names)) {
    parts <- strsplit(names[pos], ":", fixed = TRUE)[[1L]]
    for (k in seq_along(parts)) {
      comp <- components[[parts[k]]]
      if (!is.null(comp)) {
        parts[k] <- comp[[1L]]
        key <- paste(parts, collapse = ":")
        groups[[key]] <- rbind(groups[[key]], data.frame(idx = comp[[2L]], pos = pos))
        break
      }
    }
  }
  lapply(groups, function(g) g[order(g$idx), , drop = FALSE])
}

# All orders of a vector.
.permutations <- function(v) {
  if (length(v) <= 1L) return(list(v))
  out <- list()
  for (i in seq_along(v)) for (p in .permutations(v[-i])) out[[length(out) + 1L]] <- c(v[i], p)
  out
}

# The label (shown form) that produced the term key, or NA: the same number of
# `:` parts paired one to one, text parts equal, every other key part starting
# with its label variable (`sexM` with `sex`). When several labels fit, the
# longest wins.
.label_of <- function(key, labels_shown) {
  parts <- strsplit(key, ":", fixed = TRUE)[[1L]]
  fits <- character(0)
  for (lab in labels_shown) {
    lparts <- strsplit(lab, ":", fixed = TRUE)[[1L]]
    if (length(lparts) != length(parts)) next
    for (perm in .permutations(lparts)) {
      ok <- ifelse(startsWith(perm, "ssd("), parts == perm, startsWith(parts, perm))
      if (all(ok)) {
        fits <- c(fits, lab)
        break
      }
    }
  }
  if (!length(fits)) return(NA_character_)
  fits[which.max(nchar(fits))]
}

# TRUE for a coefficient with a text component column among its `:` parts.
.is_text_coef <- function(names) {
  vapply(strsplit(names, ":", fixed = TRUE), function(p) any(startsWith(p, ".ssd_")), TRUE)
}
