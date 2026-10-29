# Port of ssdiff/wordpack/fetch.py and the WORDPACK_TAGS table of ssdiff/lang_config.py. Change them together.

# Tag -> pack file on the host, and the message shown when the tag is loaded.
.WORDPACK_TAGS <- list(
  default = list(
    file = "en-glove2024-dolma-1.wordpack",
    message = paste0(
      "No embedding named: loading the default English wordpack, en-glove2024-dolma-1 ",
      "(web pages, Reddit, books). For formal text such as news, essays and reports, ",
      "use load_wordpack(\"en-wiki\").")
  ),
  en = list(
    file = "en-glove2024-dolma-1.wordpack",
    message = paste0(
      "wordpack en: en-glove2024-dolma-1 (web pages, Reddit, books). ",
      "For formal text such as news, essays and reports, use load_wordpack(\"en-wiki\").")
  ),
  `en-dolma` = list(
    file = "en-glove2024-dolma-1.wordpack",
    message = "wordpack en-dolma: en-glove2024-dolma-1 (web pages, Reddit, books)."
  ),
  `en-wiki` = list(
    file = "en-glove2024-wikigiga-1.wordpack",
    message = paste0(
      "wordpack en-wiki: en-glove2024-wikigiga-1 (Wikipedia 2024 and newswire 1994-2010). ",
      "For informal or recent text, use load_wordpack(\"en\").")
  )
)

# Base URLs, tried in order. Every host holds the same files under the same names.
.WORDPACK_HOSTS <- "https://resources.s-on.pl/wordpacks"

# A request that fails is sent to the same host this many more times before the next host is tried.
.WORDPACK_RETRIES <- 2L

# A full name is `<name>-<version>` with a whole-number version. Pack names
# never end in a hyphen followed only by digits.
.WORDPACK_FULL_NAME <- "^([a-z0-9-]+)-([0-9]+)$"

.is_full_name <- function(x) grepl(.WORDPACK_FULL_NAME, x)

.tag_list <- function() paste(names(.WORDPACK_TAGS), collapse = ", ")

# The directory pack files are kept in (01 "The pack cache").
.wordpack_cache_dir <- function() {
  dir <- Sys.getenv("WORDPACK_CACHE")
  if (nzchar(dir)) path.expand(dir) else tools::R_user_dir("wordpack", "cache")
}

#' The wordpack cache
#'
#' Lists the packs in the directory that downloaded packs are kept in. The
#' directory is the environment variable `WORDPACK_CACHE` when it is set, else
#' `tools::R_user_dir("wordpack", "cache")`. Python ssdiff keeps its packs in
#' `~/.wordpacks`; set `WORDPACK_CACHE=~/.wordpacks` (for example in
#' `~/.Renviron`) to share one copy of each pack.
#'
#' A pack is about 380 MB. Nothing is deleted automatically: remove a pack
#' you no longer need by deleting its file.
#'
#' @return The cache directory, invisibly. The packs in it and their sizes
#'   are printed.
#' @examples
#' wordpack_cache()
#' @export
wordpack_cache <- function() {
  dir <- .wordpack_cache_dir()
  files <- sort(list.files(dir, pattern = "\\.wordpack$", full.names = TRUE), method = "radix")
  cat("Wordpack cache: ", dir, "\n", sep = "")
  if (!length(files)) {
    cat("  no packs\n")
  } else {
    mb <- formatC(round(file.size(files) / 1e6), format = "d", big.mark = ",")
    cat(sprintf("  %s  %s MB\n", format(sub("\\.wordpack$", "", basename(files))), format(mb, justify = "right")),
        sep = "")
  }
  invisible(dir)
}

# Runs `expr` (a request) and returns NULL, or the error message when it fails.
# Warnings from download.file() carry the reason (an HTTP status, say) and go
# into that message.
.fetch_try <- function(expr) {
  warned <- character(0)
  err <- tryCatch(
    withCallingHandlers({
      expr
      NULL
    }, warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e)
  )
  if (is.null(err)) NULL else paste(c(err, warned), collapse = "; ")
}

# Downloads `url` to `dest`. Cloudflare answers 403 to some default User-Agents.
.fetch_get <- function(url, dest) {
  agent <- paste0("ssdiff/", utils::packageVersion("ssdiff"))
  status <- utils::download.file(url, dest, mode = "wb", quiet = TRUE, headers = c("User-Agent" = agent))
  if (status != 0L) stop("download failed with status ", status, call. = FALSE)
}

# Says that `host` is asked again; after its last try, returns `failures` with this one added.
.fetch_failed <- function(full, host, error, attempt, failures) {
  if (attempt <= .WORDPACK_RETRIES) {
    message(sprintf("wordpack %s: %s failed (%s); trying again (%d of %d)",
                    full, host, error, attempt, .WORDPACK_RETRIES))
    failures
  } else {
    c(failures, sprintf("  %s: %s", host, error))
  }
}

# index.json from the first host that gives it in the right shape.
.fetch_index <- function(full) {
  failures <- character(0)
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp))
  for (host in .WORDPACK_HOSTS) {
    for (attempt in seq_len(.WORDPACK_RETRIES + 1L)) {
      index <- NULL
      err <- .fetch_try({
        .fetch_get(paste0(host, "/index.json"), tmp)
        index <- jsonlite::read_json(tmp)
        packs <- if (is.list(index) && !is.null(names(index))) index[["packs"]]
        if (!is.list(packs) || (length(packs) && is.null(names(packs)))) {
          stop("index.json does not have the expected contents", call. = FALSE)
        }
      })
      if (is.null(err)) return(index)
      failures <- .fetch_failed(full, host, err, attempt, failures)
    }
  }
  stop(sprintf("Wordpack '%s' is not in the cache (%s) and cannot be downloaded: ", full, .wordpack_cache_dir()),
       sprintf("no host answered for index.json in %d tries each:\n", .WORDPACK_RETRIES + 1L),
       paste(failures, collapse = "\n"), call. = FALSE)
}

# The entry of index.json for a full name: file, size, sha256, format_version.
.fetch_entry <- function(index, full) {
  parts <- regmatches(full, regexec(.WORDPACK_FULL_NAME, full))[[1]]
  name <- parts[2]
  version <- parts[3]
  pack <- index[["packs"]][[name]]
  if (is.null(pack)) {
    stop(sprintf("No wordpack named '%s'. Tags: %s.", full, .tag_list()), call. = FALSE)
  }
  versions <- pack[["versions"]]
  found <- versions[[version]]
  if (is.null(found)) {
    have <- names(versions)
    have <- have[order(as.numeric(have))]
    stop(sprintf("No wordpack named '%s'. Versions of %s: %s.", full, name, paste(have, collapse = ", ")),
         call. = FALSE)
  }
  fv <- found[["format_version"]]
  if (!(is.numeric(fv) && length(fv) == 1L && fv == .wp_format_version)) {
    stop(sprintf("Wordpack %s has format version %s; this ssdiff reads format version %d. Upgrade ssdiff.",
                 full, format(fv), .wp_format_version), call. = FALSE)
  }
  found
}

# Downloads the file `entry` describes to `dest`. It appears there only after
# its size and SHA-256 are checked.
.fetch_download <- function(full, entry, dest) {
  dir <- dirname(dest)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  size <- as.numeric(entry[["size"]])
  message(sprintf("downloading wordpack %s (%s MB) to %s", full, formatC(round(size / 1e6), format = "d", big.mark = ","), dest))
  message("Set WORDPACK_CACHE=~/.wordpacks to share packs with Python ssdiff.")
  old <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old))
  failures <- character(0)
  for (host in .WORDPACK_HOSTS) {
    for (attempt in seq_len(.WORDPACK_RETRIES + 1L)) {
      part <- tempfile(paste0(full, "-"), tmpdir = dir, fileext = ".part")
      # The .part file goes on any exit, an interrupt included; after a rename it is already gone.
      err <- tryCatch(.fetch_try({
        .fetch_get(paste0(host, "/", entry[["file"]]), part)
        got_size <- file.size(part)
        got_sha <- unname(tools::sha256sum(part))
        want_sha <- as.character(entry[["sha256"]] %||% "")
        if (!isTRUE(got_size == size) || !identical(got_sha, want_sha)) {
          stop(sprintf("got %.0f bytes with SHA-256 %s; the index lists %.0f bytes with SHA-256 %s",
                       got_size, got_sha, size, want_sha), call. = FALSE)
        }
        # Another session may have finished the same download first.
        if (!file.rename(part, dest) && !file.exists(dest)) {
          stop("cannot move the download to ", dest, call. = FALSE)
        }
      }), finally = unlink(part))
      if (is.null(err)) {
        message(sprintf("downloaded wordpack %s", full))
        return(invisible(dest))
      }
      failures <- .fetch_failed(full, host, err, attempt, failures)
    }
  }
  stop(sprintf("Cannot fetch wordpack '%s': no host gave its file in %d tries each:\n", full, .WORDPACK_RETRIES + 1L),
       paste(failures, collapse = "\n"), call. = FALSE)
}

# The path of the cached pack for a full name, downloading it when it is
# missing. A full name in the cache makes no network request.
.fetch_resolve <- function(full) {
  path <- file.path(.wordpack_cache_dir(), paste0(full, ".wordpack"))
  if (file.exists(path)) return(path)
  .fetch_download(full, .fetch_entry(.fetch_index(full), full), path)
  path
}
