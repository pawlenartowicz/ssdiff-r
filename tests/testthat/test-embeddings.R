test_that("load_wordpack() on the fixture pack: metadata, ABTT 1, unit top rows", with_empty_cache({
  emb <- load_wordpack(test_pack())
  expect_s3_class(emb, "ssd_embeddings")
  expect_equal(emb$name, "test")
  expect_equal(emb$version, "1")
  expect_equal(emb$language, "en")
  expect_equal(emb$abtt, 1)
  expect_equal(emb$neighbor_filter, ".*\\d")
  pack <- wordpack_open(test_pack())
  expect_equal(emb$stopwords, pack$stopwords)
  expect_true(is.double(emb$vectors))
  expect_equal(dim(emb$vectors), c(min(pack$n, .RAM_TOP_N), pack$dims))
  expect_equal(sqrt(rowSums(emb$vectors^2)), rep(1, pack$n))
  expect_equal(emb$vectors, wordpack_read_rows(pack, seq_len(pack$n) - 1L, abtt = 1L))
}))

test_that("load_wordpack() runs the pack's self-test", with_empty_cache({
  local_mocked_bindings(wordpack_selftest = function(pack) stop("self-test ran"))
  expect_error(load_wordpack(test_pack()), "self-test ran", fixed = TRUE)
}))

test_that("a tag loads its pack file from the cache with no request", with_empty_cache({
  local_mocked_bindings(.load_pack_file = function(path) path,
                        .fetch_get = function(url, dest) stop("requested ", url))
  dir <- Sys.getenv("WORDPACK_CACHE")
  dir.create(dir)
  for (f in c("en-glove2024-dolma-1.wordpack", "en-glove2024-wikigiga-1.wordpack")) file.create(file.path(dir, f))
  expect_message(got <- load_wordpack(), "No embedding named: loading the default English wordpack", fixed = TRUE)
  expect_equal(got, file.path(dir, "en-glove2024-dolma-1.wordpack"))
  expect_message(got <- load_wordpack("en-wiki"), "use load_wordpack(\"en\").", fixed = TRUE)
  expect_equal(got, file.path(dir, "en-glove2024-wikigiga-1.wordpack"))
  expect_equal(suppressMessages(load_wordpack("en-dolma")), file.path(dir, "en-glove2024-dolma-1.wordpack"))
}))

test_that("a full name loads from the cache, downloading it when missing", with_empty_cache({
  local_mocked_bindings(.WORDPACK_HOSTS = make_host())
  emb <- suppressMessages(load_wordpack("test-2"))
  expect_equal(emb$pack$path, normalizePath(file.path(Sys.getenv("WORDPACK_CACHE"), "test-2.wordpack")))
}))

test_that("a file that is not a .wordpack raises the wordpack-only error", {
  txt <- tempfile(fileext = ".txt")
  writeLines("dog 0.1 0.2", txt)
  expect_error(load_wordpack(txt), "R ssdiff reads .wordpack files only", fixed = TRUE)
})

test_that("a name that is no tag, file or full name lists the tags without a request", with_empty_cache({
  local_mocked_bindings(.fetch_get = function(url, dest) stop("requested ", url))
  for (x in c("xx", "EN", "english", "missing.kv", "test-1.wordpack")) {
    expect_error(load_wordpack(x), sprintf("No wordpack tag or file named '%s'. Tags: default, en, en-dolma, en-wiki.", x),
                 fixed = TRUE)
  }
}))

test_that("set_abtt() reads the rows again at m and leaves its input unchanged", with_empty_cache({
  emb <- load_wordpack(test_pack())
  before <- emb$vectors
  for (m in c(0, 3, 10)) {
    got <- set_abtt(emb, m)
    expect_equal(got$abtt, m)
    expect_equal(got$vectors, wordpack_read_rows(emb$pack, seq_len(nrow(before)) - 1L, abtt = m))
  }
  expect_equal(set_abtt(set_abtt(emb, 5), 1)$vectors, before)
  expect_equal(emb$abtt, 1)
  expect_identical(emb$vectors, before)
}))

test_that("set_abtt() takes whole numbers from 0 to 10 only", with_empty_cache({
  emb <- load_wordpack(test_pack())
  expect_error(set_abtt(emb, 11), "m must be a whole number between 0 and 10, got 11", fixed = TRUE)
  expect_error(set_abtt(emb, -1), "between 0 and 10")
  expect_error(set_abtt(emb, 1.5), "between 0 and 10")
  expect_error(set_abtt(list(), 1), "load_wordpack()", fixed = TRUE)
}))

test_that("print() shows the full name, words, dimensions, ABTT and file", with_empty_cache({
  emb <- set_abtt(load_wordpack(test_pack()), 2)
  out <- capture.output(print(emb))
  expect_match(out[1], "test-1: 48 words x 16 dims, ABTT 2", fixed = TRUE)
  expect_match(out[2], emb$pack$path, fixed = TRUE)
}))

test_that("a pack reference resolves from the registry at its ABTT", with_empty_cache({
  emb <- load_wordpack(test_pack())
  ref <- .pack_ref(set_abtt(emb, 3))
  expect_equal(ref[c("name", "version", "abtt")], list(name = "test", version = "1", abtt = 3L))
  local_mocked_bindings(.load_pack_file = function(path) stop("loaded ", path))
  got <- .embeddings_from_ref(ref)
  expect_equal(got$abtt, 3)
  expect_equal(got$vectors, set_abtt(emb, 3)$vectors)
}))

test_that("with an empty registry a reference loads its pack from the cache", with_empty_cache({
  emb <- load_wordpack(test_pack())
  ref <- .pack_ref(emb)
  rm(list = ls(.registry), envir = .registry)
  dir <- Sys.getenv("WORDPACK_CACHE")
  dir.create(dir)
  file.copy(test_pack(), file.path(dir, "test-1.wordpack"))
  got <- .embeddings_from_ref(ref)
  expect_equal(got$pack$path, normalizePath(file.path(dir, "test-1.wordpack")))
  expect_equal(got$vectors, emb$vectors)
  expect_true(exists("test-1", envir = .registry, inherits = FALSE))
}))

test_that("with an empty registry and cache a reference names its pack", with_empty_cache({
  ref <- .pack_ref(load_wordpack(test_pack()))
  rm(list = ls(.registry), envir = .registry)
  expect_error(.embeddings_from_ref(ref),
               "This result was fitted on wordpack test-1. Load it with load_wordpack(\"test-1\")",
               fixed = TRUE)
}))

test_that("word vectors come from memory or the file at the embeddings' ABTT", with_empty_cache({
  emb <- set_abtt(load_wordpack(test_pack()), 2)
  # Keep 5 rows in memory, so "dog" (row 15) and "people" (row 38) are read from the file.
  emb$vectors <- emb$vectors[1:5, ]
  got <- .word_vectors(emb, c("people", "the", "no-such-word", "dog", "the"))
  expect_equal(rownames(got), c("people", "the", "dog"))
  expect_equal(unname(got), wordpack_read_rows(emb$pack, c(38, 0, 15), abtt = 2))
}))
