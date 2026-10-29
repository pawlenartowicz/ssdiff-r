# Fixture packs and a temporary pack cache for the wordpack tests.

test_pack <- function() test_path("fixtures", "wordpack", "test.wordpack")

# Runs `code` with WORDPACK_CACHE and R_USER_CACHE_DIR pointing at new
# temporary directories and with an empty session registry; restores all
# three afterwards.
with_empty_cache <- function(code) {
  old <- Sys.getenv(c("WORDPACK_CACHE", "R_USER_CACHE_DIR"), unset = NA)
  saved <- as.list(.registry, all.names = TRUE)
  on.exit({
    for (v in names(old)) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old[v]))
    rm(list = ls(.registry, all.names = TRUE), envir = .registry)
    list2env(saved, envir = .registry)
  })
  Sys.setenv(WORDPACK_CACHE = tempfile("wordpack-cache-"), R_USER_CACHE_DIR = tempfile("r-cache-"))
  rm(list = ls(.registry, all.names = TRUE), envir = .registry)
  force(code)
}

# A host directory in the layout of the real one: two copies of the fixture
# pack listed as test-1 and test-2 in index.json. `changes`
# overrides fields of every index entry. Returns the host's file:// URL.
make_host <- function(changes = list()) {
  dir <- tempfile("wordpack-host-")
  dir.create(dir)
  src <- test_pack()
  versions <- list()
  for (version in c("1", "2")) {
    file <- paste0("test-", version, ".wordpack")
    file.copy(src, file.path(dir, file))
    entry <- list(file = file, size = file.size(src), sha256 = unname(tools::sha256sum(src)), format_version = 1L)
    entry[names(changes)] <- changes
    versions[[version]] <- entry
  }
  index <- list(index_version = 1L, packs = list(`test` = list(language = "en", versions = versions)))
  jsonlite::write_json(index, file.path(dir, "index.json"), auto_unbox = TRUE, digits = NA)
  file_url(dir)
}

file_url <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  paste0(ifelse(startsWith(path, "/"), "file://", "file:///"), path)
}

# Files in the cache directory, .part files included.
cache_files <- function() list.files(Sys.getenv("WORDPACK_CACHE"), all.files = TRUE, no.. = TRUE)
