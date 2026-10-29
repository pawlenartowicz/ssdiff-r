test_that("the cache directory is WORDPACK_CACHE when set, else R's user cache", with_empty_cache({
  expect_equal(.wordpack_cache_dir(), Sys.getenv("WORDPACK_CACHE"))
  Sys.unsetenv("WORDPACK_CACHE")
  expect_equal(.wordpack_cache_dir(), tools::R_user_dir("wordpack", "cache"))
  expect_true(startsWith(.wordpack_cache_dir(), Sys.getenv("R_USER_CACHE_DIR")))
}))

test_that("wordpack_cache() returns the directory and lists its packs with sizes", with_empty_cache({
  dir <- Sys.getenv("WORDPACK_CACHE")
  expect_output(got <- wordpack_cache(), "no packs")
  expect_equal(got, dir)
  dir.create(dir)
  writeBin(raw(3e6), file.path(dir, "test-1.wordpack"))
  expect_output(wordpack_cache(), "test-1  3 MB", fixed = TRUE)
}))

test_that("a full name not in the cache is downloaded and checked against the index", with_empty_cache({
  local_mocked_bindings(.WORDPACK_HOSTS = make_host())
  msgs <- capture_messages(path <- .fetch_resolve("test-2"))
  expect_true(any(grepl("Set WORDPACK_CACHE=~/.wordpacks to share packs with Python ssdiff.", msgs, fixed = TRUE)))
  expect_equal(path, file.path(Sys.getenv("WORDPACK_CACHE"), "test-2.wordpack"))
  expect_identical(unname(tools::sha256sum(path)), unname(tools::sha256sum(test_pack())))
  expect_equal(cache_files(), "test-2.wordpack")
}))

test_that("a full name in the cache makes no request", with_empty_cache({
  local_mocked_bindings(.fetch_get = function(url, dest) stop("requested ", url))
  dir <- Sys.getenv("WORDPACK_CACHE")
  dir.create(dir)
  file.copy(test_pack(), file.path(dir, "test-1.wordpack"))
  expect_equal(.fetch_resolve("test-1"), file.path(dir, "test-1.wordpack"))
}))

test_that("a download that fails its size or SHA-256 check leaves no file in the cache", with_empty_cache({
  for (changes in list(list(sha256 = strrep("0", 64)), list(size = 1L))) {
    local_mocked_bindings(.WORDPACK_HOSTS = make_host(changes))
    msgs <- character(0)
    err <- withCallingHandlers(
      tryCatch(.fetch_resolve("test-1"), error = function(e) conditionMessage(e)),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    expect_match(err, "Cannot fetch wordpack 'test-1': no host gave its file in 3 tries each", fixed = TRUE)
    expect_match(err, "SHA-256", fixed = TRUE)
    expect_equal(sum(grepl("trying again", msgs, fixed = TRUE)), 2L)
    expect_length(cache_files(), 0L)
  }
}))

test_that("a download cut off by an error or an interrupt leaves no file in the cache", with_empty_cache({
  host <- make_host()
  local_mocked_bindings(.WORDPACK_HOSTS = host)
  index <- suppressMessages(.fetch_index("test-1"))
  cut <- function(signal) function(url, dest) {
    writeBin(raw(100), dest)
    signal()
  }
  entry <- .fetch_entry(index, "test-1")
  dest <- file.path(Sys.getenv("WORDPACK_CACHE"), "test-1.wordpack")
  local_mocked_bindings(.fetch_get = cut(function() stop("connection reset")))
  expect_error(suppressMessages(.fetch_download("test-1", entry, dest)), "connection reset", fixed = TRUE)
  expect_length(cache_files(), 0L)

  interrupt <- structure(class = c("interrupt", "condition"), list(message = "", call = NULL))
  local_mocked_bindings(.fetch_get = cut(function() signalCondition(interrupt)))
  got <- tryCatch(suppressMessages(.fetch_download("test-1", entry, dest)), interrupt = function(e) "interrupted")
  expect_identical(got, "interrupted")
  expect_length(cache_files(), 0L)
}))

test_that("a format version this package cannot read asks for an upgrade", with_empty_cache({
  local_mocked_bindings(.WORDPACK_HOSTS = make_host(list(format_version = 2L)))
  expect_error(.fetch_resolve("test-1"), "Upgrade ssdiff", fixed = TRUE)
  expect_length(cache_files(), 0L)
}))

test_that("an unknown full name lists the versions or the tags", with_empty_cache({
  local_mocked_bindings(.WORDPACK_HOSTS = make_host())
  expect_error(.fetch_resolve("test-7"), "No wordpack named 'test-7'. Versions of test: 1, 2.",
               fixed = TRUE)
  expect_error(.fetch_resolve("pl-nothing-1"), "No wordpack named 'pl-nothing-1'. Tags: default, en, en-dolma, en-wiki.",
               fixed = TRUE)
}))

test_that("with no host and an empty cache the error names every host", with_empty_cache({
  hosts <- file_url(file.path(tempdir(), c("gone-a", "gone-b")))
  local_mocked_bindings(.WORDPACK_HOSTS = hosts)
  err <- tryCatch(suppressMessages(.fetch_resolve("test-1")), error = conditionMessage)
  expect_match(err, "Wordpack 'test-1' is not in the cache (", fixed = TRUE)
  expect_match(err, "no host answered for index.json in 3 tries each", fixed = TRUE)
  for (h in hosts) expect_match(err, h, fixed = TRUE)
}))

test_that("a host that answers index.json with something else is skipped", with_empty_cache({
  portal <- tempfile("portal-")
  dir.create(portal)
  good <- make_host()
  for (answer in c("<html>Sign in to this network</html>", "5", "{\"status\": \"blocked\"}")) {
    writeLines(answer, file.path(portal, "index.json"))
    local_mocked_bindings(.WORDPACK_HOSTS = c(file_url(portal), good))
    expect_equal(suppressMessages(.fetch_resolve("test-1")),
                 file.path(Sys.getenv("WORDPACK_CACHE"), "test-1.wordpack"))
    unlink(file.path(Sys.getenv("WORDPACK_CACHE"), "test-1.wordpack"))
  }
}))

test_that("every tag names a pack file under a full name", {
  for (tag in names(.WORDPACK_TAGS)) {
    expect_false(grepl(".", tag, fixed = TRUE))
    file <- .WORDPACK_TAGS[[tag]]$file
    expect_true(endsWith(file, ".wordpack"))
    expect_true(.is_full_name(sub("\\.wordpack$", "", file)))
  }
})
