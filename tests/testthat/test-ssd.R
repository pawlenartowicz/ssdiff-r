test_that("strings and the same tokens as token lists give the same vectors", with_empty_cache({
  emb <- load_wordpack(test_pack())
  d <- golden_texts()
  texts <- d$text[1:40]
  from_strings <- ssd(texts, emb, lexicon = d$lexicon)
  from_tokens <- ssd(from_strings$docs, emb, lexicon = d$lexicon)
  expect_identical(from_strings$tokenizer, "wordpack:test-1")
  expect_length(from_strings$pre_docs, 40L)
  expect_null(from_tokens$tokenizer)
  expect_null(from_tokens$pre_docs)
  expect_identical(from_tokens$docs, from_strings$docs)
  expect_identical(from_tokens$x, from_strings$x)
  expect_identical(from_tokens$keep, from_strings$keep)
}))

test_that("use_full_doc defaults to TRUE without a lexicon and FALSE with one", with_empty_cache({
  emb <- load_wordpack(test_pack())
  docs <- list(c("dog", "cat", "swim"), c("blue", "green"))
  expect_true(ssd(docs, emb)$use_full_doc)
  expect_equal(ssd(docs, emb)$keep, c(TRUE, TRUE))
  s <- ssd(docs, emb, lexicon = "dog")
  expect_false(s$use_full_doc)
  expect_equal(s$keep, c(TRUE, FALSE))
  # A lexicon with full-document vectors stays possible.
  expect_equal(ssd(docs, emb, lexicon = "dog", use_full_doc = TRUE)$x, ssd(docs, emb)$x)
}))

test_that("embeddings = NULL loads the default pack", with_empty_cache({
  emb <- load_wordpack(test_pack())
  local_mocked_bindings(load_wordpack = function(x = "default") {
    expect_equal(x, "default")
    emb
  })
  s <- ssd(list(c("dog", "cat")))
  expect_equal(s$pack, .pack_ref(emb))
}))

test_that("ssd() checks texts and lexicon before it loads the default pack", with_empty_cache({
  local_mocked_bindings(load_wordpack = function(x = "default") stop("the pack was loaded"))
  expect_error(ssd(character(0)), "texts is empty", fixed = TRUE)
  expect_error(ssd(1:3), "texts must be a character vector", fixed = TRUE)
  expect_error(ssd("a dog", lexicon = 1), "lexicon must be a character vector of words", fixed = TRUE)
}))

test_that("a factor of texts is read as its labels", with_empty_cache({
  emb <- load_wordpack(test_pack())
  texts <- c("the dog and the cat", "blue", "the dog")
  expect_identical(ssd(factor(texts), emb, lexicon = "dog"), ssd(texts, emb, lexicon = "dog"))
}))

test_that("a text without a lexicon holds tokens only, and every fit raises", with_empty_cache({
  emb <- load_wordpack(test_pack())
  docs <- list(c("dog", "cat"), c("cat", "duck"), c("dog", "swim"))
  s <- ssd(docs, emb, use_full_doc = FALSE)
  expect_null(s$x)
  expect_identical(s$docs, docs)
  msg <- "This text has no lexicon. Use with_lexicon(s, lexicon) or use_full_doc = TRUE."
  for (kind in c("pls", "regression", "multipls", "groups", "groups_median")) {
    expect_error(.outcome_rows(s, c(1, 2, 3), NULL, kind), msg, fixed = TRUE)
  }
  # The lexicon tools need tokens only.
  lex <- suggest_lexicon(s, c(1, 2, 3), min_docs = 1)
  expect_setequal(lexicon_tokens(lex), c("dog", "cat", "duck", "swim"))
  expect_equal(evaluate_lexicon(s, NULL, "dog")$summary$docs_any, 2L)
}))

test_that("with_lexicon() equals ssd() with that lexicon and uses ssd()'s defaults", with_empty_cache({
  emb <- load_wordpack(test_pack())
  d <- golden_texts()
  texts <- d$text[1:40]
  base <- ssd(texts, emb, window = 5, sif_a = NULL)
  expect_identical(with_lexicon(base, d$lexicon), ssd(texts, emb, lexicon = d$lexicon))
  expect_identical(with_lexicon(base, d$lexicon, window = 2, sif_a = NULL, use_full_doc = TRUE),
                   ssd(texts, emb, lexicon = d$lexicon, window = 2, sif_a = NULL, use_full_doc = TRUE))
}))

test_that("the object keeps a reference to its pack, not the pack's rows", with_empty_cache({
  emb <- set_abtt(load_wordpack(test_pack()), 3)
  s <- ssd(list(c("dog", "cat")), emb)
  expect_equal(s$pack, list(name = "test", version = "1", abtt = 3L, path = emb$pack$path))
  expect_false(any(vapply(s, function(f) inherits(f, c("ssd_embeddings", "wordpack")), TRUE)))
  # with_lexicon() finds the rows again at the text's ABTT.
  expect_identical(with_lexicon(s, "dog", use_full_doc = TRUE)$x, s$x)
}))

test_that("window and sif_a are checked as in Python", with_empty_cache({
  emb <- load_wordpack(test_pack())
  docs <- list(c("dog", "cat"))
  expect_error(ssd(docs, emb, window = 0), "window must be >= 1, got 0", fixed = TRUE)
  expect_error(ssd(docs, emb, window = 1.5), "window must be >= 1", fixed = TRUE)
  expect_error(ssd(docs, emb, sif_a = 0), "sif_a must be > 0, got 0", fixed = TRUE)
  expect_error(ssd(docs, emb, sif_a = -1), "sif_a must be > 0, got -1", fixed = TRUE)
  s <- ssd(docs, emb)
  expect_error(with_lexicon(s, "dog", window = 0), "window must be >= 1", fixed = TRUE)
  expect_error(with_lexicon(s, "dog", sif_a = 0), "sif_a must be > 0", fixed = TRUE)
  expect_error(ssd(list(), emb), "texts is empty", fixed = TRUE)
}))

test_that("print() shows documents, kept documents, D, |L| and the tokenizer", with_empty_cache({
  emb <- load_wordpack(test_pack())
  out <- capture.output(print(ssd(c("the dog and the cat", "blue", NA), emb, lexicon = c("dog", "x"))))
  expect_match(out[1], "3 documents, 1 with a vector, D = 16, |L| = 2", fixed = TRUE)
  expect_match(out[2], "wordpack:test-1", fixed = TRUE)
  expect_match(out[3], "fit_pls()  fit_multipls()  fit_regression()  fit_groups()", fixed = TRUE)
  out <- capture.output(print(ssd(list("dog"), emb, use_full_doc = FALSE)))
  expect_match(out[1], "1 documents, no vectors (tokens only), |L| = 0", fixed = TRUE)
}))
