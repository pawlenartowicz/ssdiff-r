test_that("the pack's neighbor_filter hides words, matched at their start as Python's re.match", with_empty_cache({
  emb <- load_wordpack(test_pack())
  expect_equal(emb$neighbor_filter, ".*\\d")
  expect_equal(.skip_neighbor(emb, c("covid-19", "4chan", "dog")), c(TRUE, TRUE, FALSE))
  emb$neighbor_filter <- "\\d"
  expect_equal(.skip_neighbor(emb, c("a1", "1a")), c(FALSE, TRUE))
  # Without a filter: a digit anywhere, or an uppercase first letter.
  emb$neighbor_filter <- NULL
  expect_equal(.skip_neighbor(emb, c("Dog", "a1", "dog", "Émile", "café")),
               c(TRUE, TRUE, FALSE, TRUE, FALSE))
}))

test_that("the English packs' neighbor_filter hides the words Python's re.match hides", {
  # The filter of en-glove2024-dolma-1 and -wikigiga-1 (wordpack-compiler
  # languages/en); the flags are Python's. On every key of both packs R and
  # Python agree (2026-10-27).
  emb <- list(neighbor_filter = paste0(
    ".*\\d|^[A-Z\u00c0\u00c1\u00c2\u00c3\u00c4\u00c5\u00c6\u00c7\u00c8\u00c9\u00ca\u00cb\u00cc\u00cd\u00ce",
    "\u00cf\u00d0\u00d1\u00d2\u00d3\u00d4\u00d5\u00d6\u00d8\u00d9\u00da\u00db\u00dc\u00dd\u00de\u0104\u0106",
    "\u0118\u0141\u0143\u00d3\u015a\u0179\u017b\u0108\u011c\u0124\u0134\u015c\u016c\u0160\u017d\u0110]"))
  words <- c("dog", "Dog", "\u00e9lan", "\u00c9lan", "\u0141\u00f3d\u017a", "\u0142\u00f3d\u017a", "abc1", "x\u0663",
             "x\uff13", "x\u00b2", "\u03a9mega", "\u01c5x", "b2b", "\u00df", "\u015ci", "\u015di")
  expect_identical(.skip_neighbor(emb, words),
                   c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE,
                     FALSE))
})

test_that(".filtered_neighbors() leaves out filtered words and sorts by cosine", with_empty_cache({
  emb <- load_wordpack(test_pack())
  rows <- wordpack_lookup(emb$pack, c("covid-19", "4chan")) + 1L
  nb <- .filtered_neighbors(emb, colSums(emb$vectors[rows, ]), topn = 10)
  expect_equal(nrow(nb), 10L)
  expect_false(any(c("covid-19", "4chan") %in% nb$word))
  expect_false(is.unsorted(rev(nb$cos)))
  expect_equal(nb$cos, drop(emb$vectors[nb$row, ] %*% .unit_vector(colSums(emb$vectors[rows, ]))))
  expect_equal(nrow(.filtered_neighbors(emb, emb$vectors[1, ], topn = 100)), emb$pack$n - 2L)
}))

test_that("neighbour search and clustering are limited to the rows in RAM", with_empty_cache({
  emb <- load_wordpack(test_pack())
  emb$vectors <- emb$vectors[1:30, ]           # as if only 30 rows were held
  beta <- colSums(load_wordpack(test_pack())$vectors[39:41, ])
  expect_warning(nb <- .filtered_neighbors(emb, beta, topn = 20), "clamped to 30", fixed = TRUE)
  expect_true(all(nb$row <= 30))
  expect_warning(cl <- .cluster_top_neighbors(emb, beta, topn = 20, k = 3), "restrict_vocab=50000 clamped to 30")
  words <- unlist(lapply(cl, function(c) c$words$word))
  expect_true(all(wordpack_lookup(emb$pack, words) < 30))
}))

test_that(".cluster_top_neighbors(): sizes, order, the neg side and the shared generator", with_empty_cache({
  emb <- load_wordpack(test_pack())
  beta <- colSums(emb$vectors[c(16, 18, 40), ])
  pos <- .cluster_top_neighbors(emb, beta, topn = 30, k = 4)
  expect_lte(length(pos), 4L)
  expect_true(all(vapply(pos, `[[`, 0L, "size") >= 2L))
  expect_equal(sum(vapply(pos, function(c) nrow(c$words), 0L)), sum(vapply(pos, `[[`, 0L, "size")))
  expect_false(is.unsorted(-vapply(pos, `[[`, 0, "centroid_cos_beta")))
  for (c in pos) expect_false(is.unsorted(-c$words$cos_centroid))

  neg <- .cluster_top_neighbors(emb, beta, topn = 30, k = 4, side = "neg")
  flipped <- .cluster_top_neighbors(emb, -beta, topn = 30, k = 4)
  # The same query: the same clusters, cosines with beta negated, so the same order.
  expect_equal(lapply(neg, function(c) c$words$word), lapply(flipped, function(c) c$words$word))
  expect_equal(vapply(neg, `[[`, 0, "centroid_cos_beta"), -vapply(flipped, `[[`, 0, "centroid_cos_beta"))
  expect_false(is.unsorted(vapply(neg, `[[`, 0, "centroid_cos_beta")))

  set.seed(1)
  auto <- .cluster_top_neighbors(emb, beta, topn = 30)
  set.seed(2)
  expect_identical(.cluster_top_neighbors(emb, beta, topn = 30), auto)
  expect_error(.cluster_top_neighbors(emb, beta, topn = 1), "Not enough neighbors to cluster.", fixed = TRUE)
  expect_length(.cluster_top_neighbors(emb, beta, topn = 3, k = 10, min_cluster_size = 1), 3L)
}))
