# Port of ssdiff/embeddings.py (the wordpack part). Change the two together.

# Rows of a pack kept in memory: the most frequent rows, and the default
# vocabulary limit of neighbour search and clustering. Python's _RAM_TOP_N;
# change them together.
.RAM_TOP_N <- 50000L

# "<name>-<version>" -> the loaded ssd_embeddings at ABTT 1. Results keep a
# reference to their pack and find its rows here.
.registry <- new.env(parent = emptyenv())

#' Load a wordpack
#'
#' Loads word embeddings from a `.wordpack` file, by tag, by file or by the
#' pack's full name. `x` is read in this order:
#'
#' 1. A tag (a string with no dot): `"default"`, `"en"`, `"en-dolma"` or
#'    `"en-wiki"`. The tag's message is shown and its pack is loaded.
#' 2. A path where a file exists. It must end in `.wordpack`.
#' 3. A full name such as `"en-glove2024-dolma-1"`: the pack is loaded from
#'    the cache ([wordpack_cache()]) and downloaded there first if it is
#'    missing. A cached pack loads with no network request.
#'
#' The pack must pass its self-test. It loads at ABTT 1 with unit rows. The
#' `min(n, 50000)` most frequent rows are kept in memory (about 120 MB for a
#' 300-dimensional pack); other rows are read from the file when a text needs
#' them, so the file must stay in place. Each loaded pack stays in a session
#' registry, from which saved results find it again, so its memory is freed
#' only when the R session ends, not by `rm()`.
#'
#' @param x A tag, a path to a `.wordpack` file, or a pack's full name.
#' @return An object of class `ssd_embeddings`, a list with `pack` (the open
#'   [wordpack][wordpack_open]), `name`, `version`, `language`, `abtt`,
#'   `vectors` (the in-memory rows, unit length), `stopwords` and
#'   `neighbor_filter` (a regular expression for words hidden from neighbour
#'   lists).
#' @seealso [set_abtt()]
#' @examples
#' \dontrun{
#' emb <- load_wordpack()          # the default English pack
#' emb <- load_wordpack("en-wiki")
#' emb <- load_wordpack("path/to/pack.wordpack")
#' }
#' @export
load_wordpack <- function(x = "default") {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop("x must be a wordpack tag, a path to a .wordpack file or a pack's full name", call. = FALSE)
  }
  is_name <- !grepl(".", x, fixed = TRUE)
  if (is_name && x %in% names(.WORDPACK_TAGS)) {
    tag <- .WORDPACK_TAGS[[x]]
    message(tag$message)
    return(.load_pack_file(.fetch_resolve(sub("\\.wordpack$", "", tag$file))))
  }
  if (file.exists(x) && !dir.exists(x)) {
    if (!endsWith(tolower(x), ".wordpack")) {
      stop(sprintf("'%s': R ssdiff reads .wordpack files only.", x), call. = FALSE)
    }
    return(.load_pack_file(x))
  }
  if (is_name && .is_full_name(x)) return(.load_pack_file(.fetch_resolve(x)))
  stop(sprintf("No wordpack tag or file named '%s'. Tags: %s.", x, .tag_list()), call. = FALSE)
}

# Opens a pack file, runs its self-test, reads the top rows at ABTT 1 and
# adds the result to the registry.
.load_pack_file <- function(path) {
  pack <- wordpack_open(path)
  wordpack_selftest(pack)
  top <- min(pack$n, .RAM_TOP_N)
  meta <- pack$meta
  nf <- meta$neighbor_filter
  emb <- structure(list(
    pack = pack,
    name = meta$name,
    version = meta$version,
    language = meta$language,
    abtt = 1L,
    vectors = wordpack_read_rows(pack, seq_len(top) - 1L, abtt = 1L),
    stopwords = pack$stopwords,
    neighbor_filter = if (is.null(nf) || !nzchar(nf)) NULL else nf
  ), class = "ssd_embeddings")
  assign(.full_name(emb), emb, envir = .registry)
  emb
}

.full_name <- function(emb) paste0(emb$name, "-", emb$version)

.check_embeddings <- function(emb, arg = "emb") {
  if (!inherits(emb, "ssd_embeddings")) {
    stop(arg, " must be embeddings from load_wordpack()", call. = FALSE)
  }
}

#' Set the ABTT of loaded embeddings
#'
#' ABTT ("all-but-the-top") removes the pack's mean and its top `m` principal
#' components from every row, then scales each row to unit length again. A
#' pack stores its rows before ABTT, so any `m` from 0 to 10 can be set, up or
#' down; the rows are read again from the file.
#'
#' @param emb Embeddings from [load_wordpack()].
#' @param m Number of components to remove, a whole number from 0 to 10.
#' @return New embeddings at ABTT `m`; `emb` is unchanged.
#' @examples
#' \dontrun{
#' emb <- load_wordpack()
#' emb3 <- set_abtt(emb, 3)
#' }
#' @export
set_abtt <- function(emb, m) {
  .check_embeddings(emb)
  limit <- nrow(emb$pack$abtt) - 1L
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m != round(m) || m < 0 || m > limit) {
    stop(sprintf("m must be a whole number between 0 and %d, got %s", limit, format(m)), call. = FALSE)
  }
  m <- as.integer(m)
  if (m == emb$abtt) return(emb)
  emb$vectors <- wordpack_read_rows(emb$pack, seq_len(nrow(emb$vectors)) - 1L, abtt = m)
  emb$abtt <- m
  emb
}

#' Print loaded embeddings
#'
#' Shows the pack's full name, its number of words and dimensions, the ABTT
#' and the file path.
#'
#' @param x Embeddings from [load_wordpack()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.ssd_embeddings <- function(x, ...) {
  cat("<ssd_embeddings ", .full_name(x), ": ", formatC(x$pack$n, format = "d", big.mark = ","), " words x ",
      x$pack$dims, " dims, ABTT ", x$abtt, ">\n", sep = "")
  cat("  file: ", x$pack$path, "\n", sep = "")
  invisible(x)
}

# What a result keeps of its pack.
.pack_ref <- function(emb) {
  list(name = emb$name, version = emb$version, abtt = emb$abtt, path = emb$pack$path)
}

# The embeddings a reference names, at its ABTT: from the registry, else
# loaded from the cache by full name.
.embeddings_from_ref <- function(ref) {
  full <- paste0(ref$name, "-", ref$version)
  emb <- .registry[[full]]
  if (is.null(emb)) {
    path <- file.path(.wordpack_cache_dir(), paste0(full, ".wordpack"))
    if (!file.exists(path)) {
      stop(sprintf("This result was fitted on wordpack %s. Load it with load_wordpack(\"%s\") or load_wordpack(\"%s\").",
                   full, full, ref$path), call. = FALSE)
    }
    emb <- .load_pack_file(path)
  }
  set_abtt(emb, ref$abtt)
}

# Unit vectors of `keys` at the embeddings' ABTT: a matrix with one row per
# distinct key of the pack, named by key. Rows outside the in-memory block are
# read from the file in one call; a token that is not a key gets no row.
.word_vectors <- function(emb, keys) {
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys)]
  rows <- wordpack_lookup(emb$pack, keys)
  keys <- keys[rows >= 0L]
  rows <- rows[rows >= 0L]
  out <- matrix(0, length(rows), emb$pack$dims, dimnames = list(keys, NULL))
  held <- rows < nrow(emb$vectors)
  out[held, ] <- emb$vectors[rows[held] + 1L, , drop = FALSE]
  if (!all(held)) out[!held, ] <- wordpack_read_rows(emb$pack, rows[!held], abtt = emb$abtt)
  out
}
