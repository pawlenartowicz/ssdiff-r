# Port of ssdiff/wordpack/reader.py. Change the two together.
# Reader for .wordpack files, format version 1: a port of the reader in the
# Python package ssdiff. Every pack carries a self-test, which holds this
# reader and every other reader of the format to the same output.
#
# Conventions on the R side: documents and sentences count from 1, and a
# token's `start`/`end` are 1-based inclusive code-point positions in its
# normalized text, so `substring(x$text[doc], start, end)` is the token's text.
# Rows count from 0, as in the file; -1 means no row.

.wp_format_version <- 1L

# Keys answered from an environment built when a pack opens. Rows are in
# frequency order, so these are the most frequent keys.
.wp_top_keys <- 50000L

# The chunk cache is emptied when it reaches this size, which bounds its memory.
.wp_chunk_cache_max <- 250000L

# R's limit on the bytes of an environment name: a longer chunk, key or prefix
# is computed and not cached.
.wp_name_max <- 10000L

.wp_itemsize <- c(F32 = 4L, I32 = 4L, I8 = 1L, U8 = 1L)
.wp_required <- c("abtt", "norms", "vocab_order", "matrix", "vocab", "stopwords", "masked", "chars", "selftest")

.wp_error <- function(...) {
  stop(structure(class = c("wordpack_error", "error", "condition"),
                 list(message = paste0(...), call = NULL)))
}

# ---------------------------------------------------------------------------
# Opening a pack
# ---------------------------------------------------------------------------

#' Open a .wordpack file
#'
#' Reads the header, the vocabulary, the character table, the stopwords, the
#' masked-word list and the self-test. Matrix rows are read on request by
#' [wordpack_read_rows()].
#'
#' @param path Path to a `.wordpack` file.
#' @return An object of class `wordpack`.
#' @export
wordpack_open <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  fsize <- file.size(path)
  con <- file(path, "rb")
  on.exit(close(con))

  raw8 <- readBin(con, "raw", 8L)
  if (length(raw8) != 8L) .wp_error(path, ": file is truncated (no header)")
  size <- sum(as.numeric(raw8) * 256^(0:7))
  # A file that is not safetensors gives a huge size here; check it before reading.
  if (size > fsize - 8) .wp_error(path, ": not a safetensors file, or truncated (header)")
  header <- rawToChar(readBin(con, "raw", size))
  Encoding(header) <- "UTF-8"
  tensors <- tryCatch(jsonlite::parse_json(header), error = function(e) {
    .wp_error(path, ": not a safetensors file (header is not JSON)")
  })
  if (!is.list(tensors) || (length(tensors) && is.null(names(tensors)))) {
    .wp_error(path, ": not a safetensors file (header is not an object)")
  }
  meta <- tensors[["__metadata__"]]
  tensors[["__metadata__"]] <- NULL
  if (is.null(meta)) meta <- list()
  if (!is.list(meta)) .wp_error(path, ": not a safetensors file (__metadata__ is not an object)")
  .wp_check_meta(meta, path)
  meta <- .wp_decode_meta(meta, path)

  missing <- setdiff(.wp_required, names(tensors))
  if (length(missing)) .wp_error(path, ": missing tensors [", paste0("'", sort(missing), "'", collapse = ", "), "]")
  for (name in names(tensors)) .wp_check_entry(name, tensors[[name]], path)
  base <- 8 + size
  # Tensors may be stored in any order, so the file's end is the largest tensor end.
  data_end <- base + max(vapply(tensors, function(e) as.numeric(e$data_offsets[[2]]), 0))
  if (data_end > fsize) .wp_error(path, ": file is truncated")

  bytes <- function(name) {
    span <- as.numeric(unlist(tensors[[name]]$data_offsets))
    seek(con, base + span[1])
    n <- span[2] - span[1]
    data <- readBin(con, "raw", n)
    if (length(data) != n) .wp_error(path, ": file is truncated (tensor '", name, "')")
    data
  }
  utf8 <- function(raw) {
    s <- rawToChar(raw)
    Encoding(s) <- "UTF-8"
    s
  }
  lines <- function(raw) {
    if (!length(raw)) return(character(0))
    strsplit(utf8(raw), "\n", fixed = TRUE)[[1]]
  }

  shape <- as.numeric(unlist(tensors$matrix$shape))
  abtt_shape <- as.integer(unlist(tensors$abtt$shape))
  abtt <- matrix(readBin(bytes("abtt"), "numeric", prod(abtt_shape), size = 4L, endian = "little"),
                 nrow = abtt_shape[1], byrow = TRUE)
  order <- readBin(bytes("vocab_order"), "integer", prod(as.numeric(unlist(tensors$vocab_order$shape))),
                   size = 4L, endian = "little")

  pack <- new.env(parent = emptyenv())
  pack$path <- path
  pack$meta <- meta
  pack$n <- shape[1]
  pack$dims <- shape[2]
  pack$matrix_offset <- base + as.numeric(tensors$matrix$data_offsets[[1]])
  pack$abtt <- abtt
  pack$vocab <- tryCatch(.wp_vocab(bytes("vocab"), order), wordpack_error = function(e) {
    .wp_error(path, ": ", conditionMessage(e))
  })
  pack$stopwords <- lines(bytes("stopwords"))
  pack$masked <- lines(bytes("masked"))
  pack$chars <- .wp_chars(jsonlite::parse_json(utf8(bytes("chars"))))
  pack$selftest_data <- jsonlite::parse_json(utf8(bytes("selftest")))
  # Chunk text -> list(keys, rows, starts, ends), offsets 0-based within the chunk.
  # A chunk's tokens depend on nothing but its text and the pack, and most chunks repeat.
  pack$chunks <- new.env(hash = TRUE, parent = emptyenv())
  by_len <- split(pack$masked, nchar(pack$masked))
  pack$by_len <- lapply(by_len, function(w) lapply(w, utf8ToInt))
  pack$clitics <- lapply(meta$clitics, utf8ToInt)
  class(pack) <- "wordpack"
  pack
}

#' Print a wordpack
#'
#' One line: name, version, number of keys, dimensions and language.
#'
#' @param x A `wordpack` from [wordpack_open()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.wordpack <- function(x, ...) {
  cat("<wordpack ", x$meta$name, "-", x$meta$version, ": ", format(x$n, big.mark = ","),
      " keys x ", x$dims, " dims, language ", x$meta$language, ">\n", sep = "")
  invisible(x)
}

.wp_check_meta <- function(meta, path) {
  if (!identical(meta$format, "wordpack")) {
    .wp_error(path, " is not a wordpack (format=", deparse(meta$format), ")")
  }
  if (!identical(meta$format_version, as.character(.wp_format_version))) {
    .wp_error(path, " has format version ", meta$format_version,
              "; this reader reads format version ", .wp_format_version)
  }
  if (!identical(meta$keys, "surface")) {
    # A lemma-keyed pack needs a form table this reader does not apply; it would tokenize silently wrong.
    .wp_error(path, " has keys=", deparse(meta$keys), "; this reader reads surface-keyed packs only")
  }
}

.wp_decode_meta <- function(meta, path) {
  for (key in c("case_sensitive", "clitics")) {
    if (is.null(meta[[key]])) .wp_error(path, ": missing metadata '", key, "'")
  }
  for (key in c("clitics", "source")) {
    if (!is.null(meta[[key]])) {
      meta[[key]] <- tryCatch(jsonlite::parse_json(meta[[key]], simplifyVector = TRUE), error = function(e) {
        .wp_error(path, ": metadata '", key, "' is not valid JSON")
      })
    }
  }
  meta$clitics <- as.character(meta$clitics)
  meta$case_sensitive <- identical(meta$case_sensitive, "true")
  meta
}

.wp_check_entry <- function(name, entry, path) {
  if (!is.list(entry) || !all(c("dtype", "shape", "data_offsets") %in% names(entry))) {
    .wp_error(path, ": tensor '", name, "' needs dtype, shape and data_offsets")
  }
  if (!is.character(entry$dtype) || !(entry$dtype %in% names(.wp_itemsize))) {
    .wp_error(path, ": tensor '", name, "' has unknown dtype ", deparse(entry$dtype))
  }
  whole <- function(x) is.list(x) && all(vapply(x, function(v) is.numeric(v) && length(v) == 1L && v >= 0 && v == floor(v), TRUE))
  shape <- entry$shape
  if (!(whole(shape) && (name != "matrix" || length(shape) == 2L))) {
    .wp_error(path, ": tensor '", name, "' has invalid shape")
  }
  span <- entry$data_offsets
  if (!(whole(span) && length(span) == 2L && span[[1]] <= span[[2]])) {
    .wp_error(path, ": tensor '", name, "' has invalid data_offsets")
  }
  nbytes <- prod(as.numeric(unlist(shape))) * .wp_itemsize[[entry$dtype]]
  if (nbytes != as.numeric(span[[2]]) - as.numeric(span[[1]])) {
    .wp_error(path, ": tensor '", name, "' shape does not match its ",
              as.numeric(span[[2]]) - as.numeric(span[[1]]), " bytes")
  }
}

# ---------------------------------------------------------------------------
# Vocabulary: exact lookup and prefix test by binary search over UTF-8 bytes
# ---------------------------------------------------------------------------

# No R vector of all keys is built: for 1.29M keys it takes 84 MB and 1.2 s,
# against 11 MB for the key bytes. The key bytes, the key offsets and
# `vocab_order` are kept, and any key past the first .wp_top_keys is found by
# binary search. The comparisons need byte order, so every caller of the
# search runs under LC_COLLATE=C (.wp_with_c_collate).
.wp_vocab <- function(blob, order) {
  nl <- which(blob == as.raw(10L))
  n <- length(nl)
  if (length(order) != n) .wp_error("vocab has ", n, " keys but vocab_order has ", length(order), " entries")
  v <- new.env(parent = emptyenv())
  v$blob <- blob
  v$start <- c(1L, nl[-n] + 1L)
  v$end <- nl - 1L
  v$order <- order
  v$n <- n
  top <- min(n, .wp_top_keys)
  v$rows <- new.env(hash = TRUE, parent = emptyenv(), size = 2L * top)
  if (top > 0L) {
    s <- rawToChar(blob[seq_len(nl[top] - 1L)])
    Encoding(s) <- "UTF-8"
    if (!validUTF8(s)) .wp_error("vocab is not valid UTF-8")
    keys <- strsplit(s, "\n", fixed = TRUE)[[1]]
    keep <- nzchar(keys)
    list2env(stats::setNames(as.list(seq_along(keys)[keep] - 1L), keys[keep]), envir = v$rows)
  }
  v$prefixes <- new.env(hash = TRUE, parent = emptyenv())
  v
}

.wp_key <- function(v, row) {
  a <- v$start[row + 1L]
  b <- v$end[row + 1L]
  s <- if (b >= a) rawToChar(v$blob[a:b]) else ""
  Encoding(s) <- "UTF-8"
  s
}

# Every .wp_block-th key in sorted order, as R strings: the search compares
# these first, so it extracts keys from the bytes only inside one block.
.wp_block <- 32L

# 0-based position of the first sorted key >= key.
.wp_search <- function(v, key) {
  if (is.null(v$sample)) {
    v$sample <- vapply(v$order[seq.int(1L, v$n, by = .wp_block)], function(r) .wp_key(v, r), "")
  }
  sample <- v$sample
  lo <- 0L
  hi <- length(sample)
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    if (sample[mid + 1L] < key) lo <- mid + 1L else hi <- mid
  }
  # Sample lo - 1 sorts before key and sample lo does not, so the answer lies
  # after the first and at most at the second.
  hi <- min(lo * .wp_block, v$n)
  lo <- if (lo == 0L) 0L else (lo - 1L) * .wp_block + 1L
  order <- v$order
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    if (.wp_key(v, order[mid + 1L]) < key) lo <- mid + 1L else hi <- mid
  }
  lo
}

.wp_row <- function(v, key) {
  if (!nzchar(key)) return(-1L)
  cached <- nchar(key, "bytes") <= .wp_name_max
  row <- if (cached) v$rows[[key]]
  if (is.null(row)) {
    i <- .wp_search(v, key)
    row <- if (i < v$n && .wp_key(v, v$order[i + 1L]) == key) v$order[i + 1L] else -1L
    if (cached) assign(key, row, envir = v$rows)
  }
  row
}

.wp_has_prefix <- function(v, prefix) {
  cached <- nchar(prefix, "bytes") <= .wp_name_max
  # A key is its own prefix, so a known key needs no search.
  row <- if (cached) v$rows[[prefix]]
  if (!is.null(row) && row >= 0L) return(TRUE)
  hit <- if (cached) v$prefixes[[prefix]]
  if (is.null(hit)) {
    i <- .wp_search(v, prefix)
    hit <- i < v$n && startsWith(.wp_key(v, v$order[i + 1L]), prefix)
    if (cached) assign(prefix, hit, envir = v$prefixes)
  }
  hit
}

.wp_with_c_collate <- function(expr) {
  old <- Sys.getlocale("LC_COLLATE")
  if (old != "C") {
    invisible(Sys.setlocale("LC_COLLATE", "C"))
    on.exit(invisible(Sys.setlocale("LC_COLLATE", old)))
  }
  expr
}

#' Rows of exact keys
#'
#' @param pack A `wordpack` from [wordpack_open()].
#' @param keys Character vector of keys, matched exactly (no case change).
#' @return Integer vector of 0-based rows, -1 where a key is absent.
#' @export
wordpack_lookup <- function(pack, keys) {
  keys <- enc2utf8(as.character(keys))
  if (anyNA(keys)) stop("keys has missing values.")
  .wp_with_c_collate(vapply(keys, function(k) .wp_row(pack$vocab, k), 1L, USE.NAMES = FALSE))
}

# ---------------------------------------------------------------------------
# Characters
# ---------------------------------------------------------------------------

# Character classes and maps from the pack's `chars` table. Taken from the
# pack, never from R's own Unicode data, so readers built on different
# Unicode versions tokenize the same text alike.
.wp_chars <- function(tab) {
  ranges <- function(x) {
    m <- matrix(as.integer(unlist(x)), ncol = 2L, byrow = TRUE)
    list(lo = m[, 1], hi = m[, 2])
  }
  letter <- ranges(tab$letter)
  digit <- ranges(tab$digit)
  space <- ranges(tab$space)
  # Class per code point, indexed by code point + 1: 1 letter, 2 digit, 3 space, 0 other.
  top <- max(letter$hi, digit$hi, space$hi)
  kind <- integer(top + 1L)
  for (k in 1:3) {
    r <- list(letter, digit, space)[[k]]
    kind[unlist(Map(seq.int, r$lo, r$hi)) + 1L] <- k
  }
  map <- function(m) {
    from <- as.integer(names(m))
    to <- lapply(unname(m), function(s) utf8ToInt(enc2utf8(s)))
    index <- rep(NA_integer_, max(from) + 1L)
    index[from + 1L] <- seq_along(from)
    # Mapped code point per code point + 1: itself when unmapped, NA when the
    # mapping is empty or longer than one code point.
    one <- lengths(to) == 1L
    single <- seq.int(0L, max(from))
    single[from[one] + 1L] <- unlist(to[one])
    single[from[!one] + 1L] <- NA_integer_
    list(from = from, to = to, index = index, single = single)
  }
  cls <- function(r) {
    paste0(ifelse(r$lo == r$hi, sprintf("\\x{%x}", r$lo), sprintf("\\x{%x}-\\x{%x}", r$lo, r$hi)), collapse = "")
  }
  S <- paste0("[", cls(space), "]")
  NS <- paste0("[^", cls(space), "]")
  fold <- map(tab$fold)
  lower <- map(tab$lower)
  skip <- paste0("https?://", NS, "+|(?<!", NS, ")@", NS, "+")
  list(
    kind = kind,
    lower = lower,
    apostrophe = c(39L, lower$from[vapply(lower$to, function(x) length(x) > 0L && x[1] == 39L, TRUE)]),
    upper = map(tab$upper),
    fold = fold,
    fold_any = paste0("(*UTF)[", cls(list(lo = fold$from, hi = fold$from)), "]"),
    # Only a whole run of . ! ? can match, so the lookbehind changes no match.
    # Without it a run that no whitespace follows is tried again from each of
    # its positions, which takes quadratic time.
    sent_end = paste0("(*UTF)(?<![.!?])[.!?]+(?=", S, "|\\z)|\\n"),
    skip = paste0("(*UTF)", skip),
    # Skipped spans come first, so a star inside one is never matched as a star run.
    stars = paste0("(*UTF)", skip, "|(\\*+)|\\n")
  )
}

.wp_kind <- function(chars, cp) {
  k <- chars$kind[cp + 1L]
  k[is.na(k)] <- 0L
  k
}

# Replaces each code point by its mapping; a mapping may be empty or longer
# than one code point.
.wp_map <- function(cp, m) {
  out <- m$single[cp + 1L]
  if (!anyNA(out)) return(out)
  # Code points past the table, which map to themselves, and mappings that are
  # empty or longer than one code point.
  j <- m$index[cp + 1L]
  hit <- which(!is.na(j))
  out <- as.list(cp)
  out[hit] <- m$to[j[hit]]
  as.integer(unlist(out, use.names = FALSE))
}

# ---------------------------------------------------------------------------
# Tokenizing
# ---------------------------------------------------------------------------

#' Tokenize texts with a wordpack
#'
#' Applies the pack's tokenizer: NFC, the fold map, emphasis stars,
#' sentences, skipped URL and handle spans, vocabulary matches, masked words,
#' clitics and hyphen parts.
#'
#' @param pack A `wordpack` from [wordpack_open()].
#' @param texts Character vector of texts.
#' @return A list with
#'   * `text`: the texts after NFC, fold and emphasis blanking, which every
#'     position indexes;
#'   * `sents`: a data frame `doc`, `sent`, `start`, `end`, one row per sentence;
#'   * `tokens`: a data frame `doc`, `sent`, `key`, `row`, `start`, `end`, one
#'     row per token, in text order.
#' @export
wordpack_tokenize <- function(pack, texts) {
  texts <- as.character(texts)
  if (anyNA(texts)) stop("texts has missing values; pass \"\" for an empty text.")
  .wp_with_c_collate(.wp_tokenize(pack, texts))
}

.wp_tokenize <- function(pack, texts) {
  chars <- pack$chars
  text <- utf8::utf8_normalize(enc2utf8(texts))
  folded <- which(grepl(chars$fold_any, text, perl = TRUE))
  for (i in folded) text[i] <- intToUtf8(.wp_map(utf8ToInt(text[i]), chars$fold))
  starred <- which(grepl("*", text, fixed = TRUE))
  for (i in starred) text[i] <- .wp_blank_emphasis(text[i], chars)
  Encoding(text) <- "UTF-8"

  ends <- gregexpr(chars$sent_end, text, perl = TRUE)
  n_docs <- length(text)
  # Every skipped span contains "@" or "://", so only those texts need the regex.
  skips <- vector("list", n_docs)
  linked <- grepl("@", text, fixed = TRUE) | grepl("://", text, fixed = TRUE)
  skips[linked] <- gregexpr(chars$skip, text[linked], perl = TRUE)
  sents <- vector("list", n_docs)
  chunks <- vector("list", n_docs)
  for (d in seq_len(n_docs)) {
    cp <- utf8ToInt(text[d])
    n <- length(cp)
    if (!n) next
    space <- .wp_kind(chars, cp) == 3L
    # Pieces end after each sentence end; a sentence runs from a piece's first
    # non-space character to its last one, and an all-space piece is dropped.
    m <- ends[[d]]
    cuts <- if (m[1] > 0L) m + attr(m, "match.length") - 1L else integer(0)
    nonspace <- which(!space)
    if (!length(nonspace)) next
    piece <- findInterval(nonspace - 1L, cuts) # 0-based piece of each non-space position
    first <- nonspace[!duplicated(piece)]
    last <- nonspace[!duplicated(piece, fromLast = TRUE)]
    sents[[d]] <- cbind(first, last)
    # Chunks: maximal runs of non-space characters outside skipped spans. A
    # sentence end is always followed by whitespace or the text's end, so no
    # chunk or skipped span crosses a sentence.
    inside <- !space
    s <- skips[[d]]
    if (!is.null(s) && s[1] > 0L) {
      blank <- unlist(Map(seq.int, s, s + attr(s, "match.length") - 1L))
      inside[blank] <- FALSE
    }
    prev <- c(FALSE, inside[-n])
    nxt <- c(inside[-1], FALSE)
    cs <- which(inside & !prev)
    ce <- which(inside & !nxt)
    chunks[[d]] <- cbind(cs, ce, findInterval(cs, first))
  }

  n_sents <- vapply(sents, function(x) if (is.null(x)) 0L else nrow(x), 1L)
  sent_mat <- do.call(rbind, sents)
  sent_df <- data.frame(
    doc = rep(seq_len(n_docs), n_sents),
    sent = sequence(n_sents),
    start = if (is.null(sent_mat)) integer(0) else sent_mat[, 1],
    end = if (is.null(sent_mat)) integer(0) else sent_mat[, 2]
  )

  n_chunks <- vapply(chunks, function(x) if (is.null(x)) 0L else nrow(x), 1L)
  chunk_mat <- do.call(rbind, chunks)
  if (is.null(chunk_mat)) chunk_mat <- matrix(integer(0), 0L, 3L)
  chunk_doc <- rep(seq_len(n_docs), n_chunks)
  chunk_text <- substring(text[chunk_doc], chunk_mat[, 1], chunk_mat[, 2])
  distinct <- unique(chunk_text)
  res <- .wp_chunk_results(pack, distinct)
  at <- match(chunk_text, distinct)
  field <- function(k) lapply(res, `[[`, k)[at]
  n_tok <- lengths(lapply(res, `[[`, 1L))[at]
  base <- rep(chunk_mat[, 1], n_tok)
  tokens <- data.frame(
    doc = rep(chunk_doc, n_tok),
    sent = rep(chunk_mat[, 3], n_tok),
    key = as.character(unlist(field(1L), use.names = FALSE)),
    row = as.integer(unlist(field(2L), use.names = FALSE)),
    start = base + as.integer(unlist(field(3L), use.names = FALSE)),
    end = base + as.integer(unlist(field(4L), use.names = FALSE)) - 1L,
    stringsAsFactors = FALSE
  )
  list(text = text, sents = sent_df, tokens = tokens)
}

# Tokens of each distinct chunk, from the pack's cache or by walking the chunk.
.wp_chunk_results <- function(pack, distinct) {
  cache <- pack$chunks
  short <- nchar(distinct, "bytes") <= .wp_name_max
  got <- vector("list", length(distinct))
  if (any(short)) got[short] <- mget(distinct[short], envir = cache, ifnotfound = list(NULL))
  todo <- which(vapply(got, is.null, TRUE))
  if (length(todo)) {
    new <- lapply(distinct[todo], function(c) .wp_chunk_tokens(pack, c))
    got[todo] <- new
    todo <- todo[short[todo]]
    if (length(cache) + length(todo) > .wp_chunk_cache_max) {
      rm(list = ls(cache, all.names = TRUE, sorted = FALSE), envir = cache)
    }
    keep <- todo[seq_len(min(length(todo), .wp_chunk_cache_max))]
    list2env(stats::setNames(got[keep], distinct[keep]), envir = cache)
  }
  unname(got)
}

# `text` with paired emphasis stars replaced by spaces. Left in, the closing
# stars of `**bit**` would join the masked run and be read as masked letters.
# Spaces, not deletion, keep every offset in place.
.wp_blank_emphasis <- function(text, chars) {
  cp <- utf8ToInt(text)
  n <- length(cp)
  m <- gregexpr(chars$stars, text, perl = TRUE)[[1]]
  if (m[1] < 0L) return(text)
  len <- attr(m, "match.length")
  star <- attr(m, "capture.start")[, 1] > 0L
  kind <- .wp_kind(chars, cp)
  # One-character tests on the neighbours; outside the text both fail, so the
  # text's ends count as whitespace.
  nonspace <- function(i) i >= 1L && i <= n && kind[i] != 3L
  alnum <- function(i) i >= 1L && i <= n && (kind[i] == 1L || kind[i] == 2L)
  open <- integer(0) # names: run length; value: start of the open run of that length on this line
  out <- cp
  for (j in seq_along(m)) {
    a <- m[j]
    b <- a + len[j] - 1L
    if (!star[j]) {
      if (cp[a] == 10L) open <- integer(0)
      next
    }
    key <- as.character(b - a + 1L)
    if (!is.na(open[key]) && nonspace(a - 1L) && !alnum(b + 1L)) {
      for (st in c(open[[key]], a)) out[st:(st + b - a)] <- 32L
      open <- open[names(open) != key]
    } else if (nonspace(b + 1L) && !alnum(a - 1L)) {
      open[key] <- a
    }
  }
  intToUtf8(out)
}

# The walk over one chunk. Offsets inside are 0-based, as in the Python
# reader, so each step reads like its counterpart there; `cp[i + 1L]` is the
# code point at offset i. Returns list(keys, rows, starts, ends) with
# end-exclusive offsets.
.wp_chunk_tokens <- function(pack, chunk) {
  cp <- utf8ToInt(chunk)
  n <- length(cp)
  chars <- pack$chars
  kind <- .wp_kind(chars, cp)
  w <- new.env(parent = emptyenv())
  w$keys <- character(0)
  w$rows <- integer(0)
  w$s <- integer(0)
  w$e <- integer(0)
  emit <- function(key, row, s, e) {
    w$keys[length(w$keys) + 1L] <- key
    w$rows[length(w$rows) + 1L] <- row
    w$s[length(w$s) + 1L] <- s
    w$e[length(w$e) + 1L] <- e
  }
  # Units: a run of letters, a run of digits, or any single other character.
  newunit <- c(TRUE, kind[-1] != kind[-n] | kind[-1] == 0L)
  us <- which(newunit) - 1L
  ue <- c(us[-1], n)
  uk <- kind[us + 1L]
  nu <- length(us)
  ctx <- list(pack = pack, cp = cp, n = n, kind = kind, emit = emit, us = us, ue = ue, uk = uk)
  # A masked run needs a *, so a chunk without one skips the masked-word step.
  masking <- length(pack$by_len) > 0L && any(cp == 42L)
  # Without a character other than a letter, ' and -, no vocabulary match is possible.
  special <- any(kind != 1L & cp != 39L & cp != 45L)
  i <- 1L
  while (i <= nu) {
    u0 <- us[i]
    if (uk[i] == 0L) {
      cl <- .wp_stray_clitic(ctx, u0)
      if (is.null(cl)) {
        i <- i + 1L
        next
      }
      r <- .wp_resolve(pack, cp[(u0 + 1L):(u0 + length(cl))])
      emit(r$key, r$row, u0, u0 + length(cl))
      while (i <= nu && us[i] < u0 + length(cl)) i <- i + 1L
      next
    }
    if (special) {
      hit <- .wp_vocab_match(ctx, us, ue, i)
      if (!is.null(hit)) {
        emit(hit$key, hit$row, u0, ue[hit$j - 1L])
        i <- hit$j
        next
      }
    }
    if (uk[i] == 2L) {
      i <- i + 1L
      next
    }
    hit <- if (masking) .wp_masked_match(ctx, u0) else NULL
    if (!is.null(hit)) {
      emit(hit$key, hit$row, u0, hit$end)
      end <- hit$end
    } else {
      end <- .wp_word_end(ctx, i)
      .wp_word_tokens(ctx, u0, end)
    }
    while (i <= nu && us[i] < end) i <- i + 1L
  }
  list(w$keys, w$rows, w$s, w$e)
}

# Code points [a, b) of the chunk, 0-based.
.wp_slice <- function(ctx, a, b) if (b > a) ctx$cp[(a + 1L):b] else integer(0)

# The clitic starting with ' at `pos`, if the chunk ends or a non-letter follows it.
# Catches the `'s` left after a vocabulary match, masked word or digit run
# (`covid-19's`, `90's`), which would otherwise leave a stray word `s`.
.wp_stray_clitic <- function(ctx, pos) {
  chars <- ctx$pack$chars
  # Only a code point whose lowercase starts with ' can begin a stray clitic.
  if (!(ctx$cp[pos + 1L] %in% chars$apostrophe)) return(NULL)
  for (cl in ctx$pack$clitics) {
    end <- pos + length(cl)
    if (cl[1] != 39L) next
    seg <- .wp_map(.wp_slice(ctx, pos, min(end, ctx$n)), chars$lower)
    if (identical(seg, cl) && !(end < ctx$n && ctx$kind[end + 1L] == 1L)) return(cl)
  }
  NULL
}

.wp_case_forms <- function(pack, form) {
  chars <- pack$chars
  lower <- .wp_map(form, chars$lower)
  if (!pack$meta$case_sensitive) return(list(lower))
  forms <- list(form)
  title <- c(.wp_map(form[1], chars$upper), .wp_map(form[-1], chars$lower))
  for (f in list(lower, title)) {
    if (!any(vapply(forms, identical, TRUE, f))) forms[[length(forms) + 1L]] <- f
  }
  forms
}

# Key and row of a form (code points) under the case rule; the lowercased
# form and -1 when absent.
.wp_resolve <- function(pack, form) {
  for (f in .wp_case_forms(pack, form)) {
    key <- intToUtf8(f)
    row <- .wp_row(pack$vocab, key)
    if (row >= 0L) return(list(key = key, row = row))
  }
  list(key = intToUtf8(.wp_map(form, pack$chars$lower)), row = -1L)
}

# (index after the last unit, key, row) of the longest qualifying key from unit i, or NULL.
.wp_vocab_match <- function(ctx, us, ue, i) {
  pack <- ctx$pack
  start <- us[i]
  best <- NULL
  for (j in i:length(us)) {
    cand <- .wp_slice(ctx, start, ue[j])
    prefixed <- FALSE
    for (f in .wp_case_forms(pack, cand)) {
      if (.wp_has_prefix(pack$vocab, intToUtf8(f))) {
        prefixed <- TRUE
        break
      }
    }
    if (!prefixed) break
    if (.wp_qualifies(pack$chars, cand)) {
      r <- .wp_resolve(pack, cand)
      if (r$row >= 0L) best <- list(j = j + 1L, key = r$key, row = r$row)
    }
  }
  best
}

.wp_qualifies <- function(chars, cand) {
  k <- .wp_kind(chars, cand)
  if (!any(k == 1L) || !any(k != 1L & cand != 39L & cand != 45L)) return(FALSE)
  last <- length(cand)
  if (k[last] == 1L || k[last] == 2L) return(TRUE)
  # A final period only with another before it: keeps `u.s.`, rejects junk keys like `it.`.
  cand[last] == 46L && any(cand[-last] == 46L)
}

.wp_masked_match <- function(ctx, u0) {
  # The masked run: letters and *, starting at this unit's letter.
  longest <- max(as.integer(names(ctx$pack$by_len)))
  i <- u0 + 1L
  while (i < ctx$n && i - u0 <= longest && (ctx$kind[i + 1L] == 1L || ctx$cp[i + 1L] == 42L)) i <- i + 1L
  run <- .wp_slice(ctx, u0, i)
  if (!any(run == 42L) || length(run) > longest) return(NULL)
  low <- .wp_map(run, ctx$pack$chars$lower)
  if (length(low) != length(run)) return(NULL)
  words <- ctx$pack$by_len[[as.character(length(run))]]
  hits <- Filter(function(wd) all(low == 42L | low == wd), words)
  if (length(hits) != 1L) return(NULL)
  r <- .wp_resolve(ctx$pack, hits[[1]])
  list(key = r$key, row = r$row, end = u0 + length(run))
}

# End of the longest match of L+(['-]L+)* starting at letter unit i. Letter
# units are maximal runs and other units single characters, so the word grows
# by a ' or - unit followed by a letter unit.
.wp_word_end <- function(ctx, i) {
  uk <- ctx$uk
  us <- ctx$us
  nu <- length(uk)
  while (i + 2L <= nu && uk[i + 1L] == 0L && ctx$cp[us[i + 1L] + 1L] %in% c(39L, 45L) && uk[i + 2L] == 1L) {
    i <- i + 2L
  }
  ctx$ue[i]
}

.wp_word_tokens <- function(ctx, ws, we) {
  word <- .wp_slice(ctx, ws, we)
  nw <- length(word)
  # Lowercased once when every code point maps to one code point; then the
  # lowercased tail is the tail of the lowercased word.
  low <- ctx$pack$chars$lower$single[word + 1L]
  whole <- !anyNA(low)
  cut <- we
  for (cl in ctx$pack$clitics) {
    k <- length(cl)
    if (nw <= k) next
    tail <- if (whole) low[(nw - k + 1L):nw] else .wp_map(word[(nw - k + 1L):nw], ctx$pack$chars$lower)
    if (length(tail) == k && all(tail == cl)) {
      cut <- we - k
      break
    }
  }
  .wp_stem_tokens(ctx, ws, cut)
  if (cut < we) {
    r <- .wp_resolve(ctx$pack, .wp_slice(ctx, cut, we))
    ctx$emit(r$key, r$row, cut, we)
  }
}

.wp_stem_tokens <- function(ctx, ss, se) {
  stem <- .wp_slice(ctx, ss, se)
  r <- .wp_resolve(ctx$pack, stem)
  if (r$row >= 0L || !any(stem == 45L)) {
    ctx$emit(r$key, r$row, ss, se)
    return(invisible())
  }
  # Each part between hyphens is a word of its own; a part has no '-', so the
  # recursion ends at the next level. An empty part gives no token.
  cuts <- c(ss - 1L, ss + which(stem == 45L) - 1L, se)
  for (p in seq_len(length(cuts) - 1L)) {
    a <- cuts[p] + 1L
    b <- cuts[p + 1L]
    if (b > a) .wp_word_tokens(ctx, a, b)
  }
}

# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------

.wp_l2_rows <- function(x) {
  norms <- sqrt(rowSums(x * x))
  norms[norms == 0] <- 1
  x / norms
}

# Int8 rows to unit rows with ABTT `m` applied. `abtt` row 1 is the mean,
# rows 2..11 the principal components. Computed in double precision; the
# Python reader works in float32, and the self-test allows for the difference.
.wp_dequantize <- function(q, abtt, m) {
  limit <- nrow(abtt) - 1L
  if (m < 0L || m > limit) stop("abtt must be between 0 and ", limit, ", got ", m)
  x <- .wp_l2_rows(q * 1)
  if (m > 0L) {
    x <- sweep(x, 2L, abtt[1, ])
    comps <- abtt[2:(m + 1L), , drop = FALSE]
    x <- x - (x %*% t(comps)) %*% comps
    x <- .wp_l2_rows(x)
  }
  x
}

#' Read embedding rows
#'
#' @param pack A `wordpack` from [wordpack_open()].
#' @param rows 0-based rows, in the order wanted.
#' @param abtt Number of ABTT components to remove (0 to 10).
#' @return A numeric matrix, one row per entry of `rows`.
#' @export
wordpack_read_rows <- function(pack, rows, abtt = 1L) {
  idx <- as.numeric(rows)
  if (!length(idx)) return(matrix(numeric(0), 0L, pack$dims))
  if (anyNA(idx) || min(idx) < 0 || max(idx) >= pack$n) stop("rows must be between 0 and ", pack$n - 1)
  uniq <- sort(unique(idx))
  inverse <- match(idx, uniq)
  q <- matrix(0L, length(uniq), pack$dims)
  # Each contiguous run of sorted unique rows is one seek + read.
  cuts <- which(diff(uniq) != 1)
  starts <- c(1L, cuts + 1L)
  stops <- c(cuts, length(uniq))
  con <- file(pack$path, "rb")
  on.exit(close(con))
  for (r in seq_along(starts)) {
    a <- starts[r]
    b <- stops[r]
    nbytes <- (b - a + 1) * pack$dims
    seek(con, pack$matrix_offset + uniq[a] * pack$dims)
    data <- readBin(con, "integer", nbytes, size = 1L, signed = TRUE)
    if (length(data) != nbytes) .wp_error(pack$path, ": file is truncated (matrix)")
    q[a:b, ] <- matrix(data, ncol = pack$dims, byrow = TRUE)
  }
  .wp_dequantize(q, pack$abtt, as.integer(abtt))[inverse, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

#' Run the pack's self-test
#'
#' Tokenizes the pack's stored texts and reads its stored rows, and stops at
#' the first difference from the stored output.
#'
#' @param pack A `wordpack` from [wordpack_open()].
#' @return `pack`, invisibly.
#' @export
wordpack_selftest <- function(pack) {
  label <- paste0(pack$meta$name, "-", pack$meta$version, " (", pack$path, ")")
  prefix <- paste0(label, ": self-test failed (reader format version ", .wp_format_version, ")")
  st <- pack$selftest_data
  texts <- vapply(st$texts, function(e) e$text, "")
  got <- wordpack_tokenize(pack, texts)
  pairs <- function(x) if (length(x)) do.call(rbind, lapply(x, function(p) as.integer(unlist(p)))) else matrix(integer(0), 0L, 2L)
  for (i in seq_along(st$texts)) {
    e <- st$texts[[i]]
    s <- got$sents[got$sents$doc == i, ]
    t <- got$tokens[got$tokens$doc == i, ]
    fields <- list(
      sents = list(cbind(s$start - 1L, s$end), pairs(e$sents)),
      keys = list(t$key, as.character(unlist(e$keys))),
      spans = list(cbind(t$start - 1L, t$end), pairs(e$spans)),
      sent = list(t$sent - 1L, as.integer(unlist(e$sent))),
      rows = list(t$row, as.integer(unlist(e$rows)))
    )
    for (f in names(fields)) {
      a <- fields[[f]][[1]]
      b <- fields[[f]][[2]]
      same <- length(a) == length(b) && all(dim(as.matrix(a)) == dim(as.matrix(b))) && all(a == b)
      if (!same) {
        .wp_error(prefix, ": text ", i - 1L, " ", deparse(e$text), ", field '", f, "': expected ",
                  paste(format(b), collapse = " "), ", got ", paste(format(a), collapse = " "))
      }
    }
  }
  for (v in st$vectors) {
    got_v <- wordpack_read_rows(pack, v$row, v$abtt)[1, ]
    diff <- max(abs(got_v - as.numeric(unlist(v$values))), 0)
    if (!isTRUE(diff <= 1e-5)) {
      .wp_error(prefix, ": row ", v$row, " at abtt=", v$abtt, " differs by ", format(diff, digits = 3))
    }
  }
  invisible(pack)
}
