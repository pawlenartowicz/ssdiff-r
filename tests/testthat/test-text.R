fixture_texts <- function(emb) vapply(emb$pack$selftest_data$texts, function(e) e$text, "")

test_that("tokens of strings are the reader's keys minus the pack's stopwords", with_empty_cache({
  emb <- load_wordpack(test_pack())
  texts <- fixture_texts(emb)
  got <- .tokenize_texts(texts, emb)
  tok <- wordpack_tokenize(emb$pack, texts)$tokens
  expected <- lapply(seq_along(texts), function(d) {
    keys <- tok$key[tok$doc == d]
    keys[!keys %in% emb$stopwords]
  })
  expect_identical(got$docs, expected)
  expect_identical(got$docs, lapply(got$pre_docs, `[[`, "doc_lemmas"))
  expect_identical(got$tokenizer, "wordpack:test-1")
}))

test_that("a chunk or key longer than an R name (10,000 bytes) is tokenized and looked up uncached", with_empty_cache({
  pack <- load_wordpack(test_pack())$pack
  tok <- wordpack_tokenize(pack, paste("hello", strrep("ab", 6000)))$tokens
  expect_identical(nchar(tok$key), c(5L, 12000L))   # as Python's reader
  expect_identical(wordpack_lookup(pack, strrep("a", 20000)), -1L)
  expect_error(wordpack_lookup(pack, NA), "keys has missing values.", fixed = TRUE)
}))

test_that("a text keeps Python's PreprocessedDoc fields with 1-based positions", with_empty_cache({
  emb <- load_wordpack(test_pack())
  p <- .tokenize_texts("I don't like the dog. It's well-known!", emb)$pre_docs[[1]]
  expect_identical(p$raw, "I don't like the dog. It's well-known!")
  expect_identical(p$sents_surface, c("I don't like the dog.", "It's well-known!"))
  expect_identical(p$sents_lemmas, list(c("do", "n't", "like", "dog"), c("it", "'s", "well-known")))
  expect_identical(p$doc_lemmas, c("do", "n't", "like", "dog", "it", "'s", "well-known"))
  expect_equal(unname(p$sent_char_spans), rbind(c(1, 21), c(23, 38)))
  expect_identical(substring(p$raw, p$sent_char_spans[, 1], p$sent_char_spans[, 2]), p$sents_surface)
  expect_equal(p$token_to_sent, c(1, 1, 1, 1, 2, 2, 2))
  # "i" (1st) and "the" (5th) are stopwords; positions count them.
  expect_equal(p$sents_kept_idx, list(c(2, 3, 4, 6), c(1, 2, 3)))
}))

test_that("NA gives an empty document", with_empty_cache({
  emb <- load_wordpack(test_pack())
  got <- .tokenize_texts(c("the dog swims", NA, ""), emb)
  expect_identical(got$docs[[2]], character(0))
  expect_identical(got$docs[[3]], character(0))
  expect_identical(got$pre_docs[[2]]$raw, "")
  expect_length(got$pre_docs[[2]]$sents_surface, 0L)
  expect_identical(got$docs[[1]], c("dog", "swims"))
}))

test_that("token lists are kept as given and have no sentence structure", {
  texts <- list(c("The", "dog", "the", "a"), NULL, character(0), "dog")
  got <- .tokenize_texts(texts, NULL)
  expect_identical(got$docs, list(c("The", "dog", "the", "a"), character(0), character(0), "dog"))
  expect_null(got$pre_docs)
  expect_null(got$tokenizer)
  expect_named(got, c("docs", "pre_docs", "tokenizer"))
})

test_that("any other input raises the two-forms error", {
  msg <- "texts must be a character vector (one text per document) or a list of character vectors"
  for (texts in list(1:3, factor("a"), list("a", 1), list(list("a")), data.frame(text = "a"))) {
    expect_error(.tokenize_texts(texts, NULL), msg, fixed = TRUE)
  }
})
