# Twelve points in three clear groups; the k-means references are Python's
# ssdiff.utils.math on the same matrix.
km_points <- function() {
  matrix(c(0, 0, 1, 0.1, 0, 1, 0, 0.2, 1, 5, 5, 0, 5.1, 5, 0, 5, 5.3, 0, -4, 2, 2, -4.2, 2, 2,
           -4, 2.1, 2.2, 1, 1, 1, 0.5, 0.2, 0.9, 4.8, 5.2, 0.1), ncol = 3, byrow = TRUE)
}

test_that(".standardize() uses the population SD and scales constant columns by 1", {
  x <- cbind(c(1, 2, 3, 4), 5, c(1, 1, 1, 1 + 1e-13))
  z <- .standardize(x)
  expect_equal(z$mean, c(2.5, 5, 1 + 2.5e-14))
  expect_equal(z$scale, c(sqrt(1.25), 1, 1))
  expect_equal(z$x[, 1], (1:4 - 2.5) / sqrt(1.25))
  expect_equal(z$x[, 2], rep(0, 4))
})

test_that(".pca_fit_transform() reconstructs the centred data; the component signs do not matter", {
  x <- km_points()
  full <- .pca_fit_transform(x, 3)
  expect_equal(full$z %*% full$components, sweep(x, 2, colMeans(x)))
  expect_equal(tcrossprod(full$components), diag(3))
  expect_equal(sum(full$evr), 1)
  # The projection onto the first two components equals the eigen-decomposition's,
  # whatever sign either method gives a component.
  two <- .pca_fit_transform(x, 2)
  e <- eigen(stats::cov(x), symmetric = TRUE)
  v <- e$vectors[, 1:2]
  expect_equal(crossprod(two$components), tcrossprod(v))
  expect_equal(two$evr, e$values[1:2] / sum(e$values))
})

test_that(".l2_normalize_rows() and .unit_vector() give unit length and leave zeros at zero", {
  m <- .l2_normalize_rows(rbind(c(3, 4), c(0, 0)))
  expect_equal(m, rbind(c(0.6, 0.8), c(0, 0)))
  expect_equal(.unit_vector(c(3, 4)), c(0.6, 0.8))
  expect_equal(.unit_vector(c(1e-13, 0)), c(0, 0))
})

test_that("the k-means++ draws follow the shared generator", {
  x <- km_points()
  # First centre: floor(u0 * 12) with seed 2137's first uniform 0.9498900428996606, so 11, row 12.
  expect_equal(.kmeans_plus_plus(x, 1, .rng(2137)), 12)
  # Python's _kmeans_plus_plus(x, 4, Rng(5)) picks rows 4, 2, 6, 9 (0-based).
  expect_equal(.kmeans_plus_plus(x, 4, .rng(5)), c(5, 3, 7, 10))
})

test_that(".kmeans() gives Python's clusters and is deterministic", {
  x <- km_points()
  km <- .kmeans(x, 3)
  expect_equal(km$labels - 1L, c(1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 0))
  expect_equal(km$inertia, 1.6264999999999985, tolerance = 1e-12)
  expect_equal(km$centers, rbind(c(4.975, 5.125, 0.025), c(0.32, 0.28, 0.98), c(-12.2, 6.1, 6.2) / 3))
  python <- list(`2` = c(1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0), `4` = c(1, 1, 1, 0, 0, 0, 2, 2, 2, 3, 1, 0),
                 `5` = c(1, 1, 1, 0, 0, 4, 2, 2, 2, 3, 1, 4))
  for (k in names(python)) expect_equal(.kmeans(x, as.integer(k), random_state = 11)$labels - 1L, python[[k]])

  set.seed(1)
  before <- .Random.seed
  expect_identical(.kmeans(x, 4, random_state = 11), .kmeans(x, 4, random_state = 11))
  expect_identical(.Random.seed, before)
  expect_error(.kmeans(x, 13), "Cannot request k=13 clusters from n=12 samples.", fixed = TRUE)
})

test_that(".kmeans() handles identical points (zero distances, empty clusters)", {
  km <- .kmeans(matrix(1, 5, 2), 3)
  expect_equal(km$labels, rep(1L, 5))
  expect_equal(km$inertia, 0)
})

test_that(".silhouette() on a hand-computed case; .kmeans_auto_k() picks Python's k", {
  d <- .pairwise_euclidean(matrix(c(0, 1, 10, 11)))
  expect_equal(d[1, ], c(0, 1, 10, 11))
  expect_equal(.silhouette(d, c(1, 1, 2, 2)), mean(c(9.5 / 10.5, 8.5 / 9.5)))
  expect_equal(.silhouette(d, c(1, 1, 1, 1)), 0)

  auto <- .kmeans_auto_k(km_points(), 2, 5)
  expect_equal(auto$k, 3L)
  expect_equal(auto$labels - 1L, c(1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1, 0))
  # k_min above k_max: Python's range() is empty, so only k_min is fitted.
  expect_identical(.kmeans_auto_k(km_points(), 4, 3)[c("labels", "k")], c(.kmeans(km_points(), 4)["labels"], k = 4L))
})
