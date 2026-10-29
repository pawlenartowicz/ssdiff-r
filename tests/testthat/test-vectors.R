# Three documents of the fixture pack's words; "zzz" is not a key. Counts over
# all documents: dog 3, cat 3, swim 1, duck 1, blue 1, zzz 1 (total 10).
hand_docs <- list(c("dog", "cat", "swim"), c("cat", "dog", "cat", "dog", "duck"), c("blue", "zzz"))

unit <- function(v) v / sqrt(sum(v^2))

test_that("document vectors of a hand-built text match values computed by hand", with_empty_cache({
  emb <- load_wordpack(test_pack())
  v <- function(w) wordpack_read_rows(emb$pack, wordpack_lookup(emb$pack, w), abtt = 1L)[1, ]
  dog <- v("dog"); cat <- v("cat"); swim <- v("swim"); duck <- v("duck"); blue <- v("blue")

  for (sif in list(1e-3, NULL)) {
    # Weights of a word seen 3 and 1 times.
    a3 <- if (is.null(sif)) 1 else sif / (sif + 3 / 10)
    a1 <- if (is.null(sif)) 1 else sif / (sif + 1 / 10)

    # Seed mode, lexicon "dog", window 2. Document 3 has no hit.
    s <- ssd(hand_docs, emb, lexicon = "dog", window = 2, sif_a = sif)
    d1 <- (a3 * cat + a1 * swim) / (a3 + a1)
    occ1 <- (2 * a3 * cat + a3 * dog) / (3 * a3)                  # hit 2: cat, cat, dog
    occ2 <- (a3 * dog + a3 * cat + a1 * duck) / (2 * a3 + a1)     # hit 4: dog, cat, duck
    expect_equal(s$keep, c(TRUE, TRUE, FALSE))
    expect_equal(s$x, rbind(unit(d1), unit((occ1 + occ2) / 2)), tolerance = 1e-12)

    # Full-document mode; "zzz" has no vector, so document 3 is "blue".
    s <- ssd(hand_docs, emb, sif_a = sif)
    f1 <- (a3 * dog + a3 * cat + a1 * swim) / (2 * a3 + a1)
    f2 <- (2 * a3 * cat + 2 * a3 * dog + a1 * duck) / (4 * a3 + a1)
    expect_equal(s$keep, c(TRUE, TRUE, TRUE))
    expect_equal(s$x, rbind(unit(f1), unit(f2), unit(blue)), tolerance = 1e-12)
  }
}))

test_that("a document with no word that has a vector is not kept", with_empty_cache({
  emb <- load_wordpack(test_pack())
  docs <- list(c("dog", "zzz"), character(0), c("zzz", "qqq"), c("cat", "dog", "zzz"))
  full <- ssd(docs, emb)
  expect_equal(full$keep, c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(full$n_kept, 2L)
  expect_equal(nrow(full$x), 2L)
  # Seed mode: "dog" in document 1 has only "zzz" in its window.
  seed <- ssd(docs, emb, lexicon = "dog", window = 1)
  expect_equal(seed$keep, c(FALSE, FALSE, FALSE, TRUE))
}))
