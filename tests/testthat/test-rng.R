# Same numbers as ssdiff/tests/utils/test_rng.py in the Python package.

# Exact decimal string of hi * 2^32 + lo, built in base-1e6 limbs.
u64_dec <- function(hi, lo) {
  carry_up <- function(l) {
    carry <- 0
    for (i in seq_along(l)) {
      t <- l[i] + carry
      carry <- floor(t / 1e6)
      l[i] <- t - carry * 1e6
    }
    while (carry > 0) {
      l <- c(l, carry - floor(carry / 1e6) * 1e6)
      carry <- floor(carry / 1e6)
    }
    l
  }
  vapply(seq_along(hi), function(k) {
    l <- carry_up(hi[k])
    for (j in 1:2) l <- carry_up(l * 65536)
    l <- carry_up(c(l[1] + lo[k], l[-1]))
    while (length(l) > 1L && l[length(l)] == 0) l <- l[-length(l)]
    l <- rev(l)
    paste0(sprintf("%.0f", l[1]), paste(sprintf("%06.0f", l[-1]), collapse = ""))
  }, "")
}

raw_dec <- function(seed, start, n) {
  z <- .rng_raw_at(seed, start, n)
  u64_dec(z$hi, z$lo)
}

test_that("raw values match the hardcoded splitmix64 vectors", {
  expect_identical(raw_dec(0, 0, 5), c(
    "16294208416658607535", "7960286522194355700", "487617019471545679",
    "17909611376780542444", "1961750202426094747"))
  expect_identical(raw_dec(2137, 0, 5), c(
    "17522378519535027745", "364756094709023243", "8236747292942766741",
    "8881960951598843736", "8369509784638573016"))
  expect_identical(raw_dec(1234567, 0, 5), c(
    "6457827717110365317", "3203168211198807973", "9817491932198370423",
    "4593380528125082431", "16408922859458223821"))
})

test_that("uniforms match the hardcoded vector", {
  expect_identical(
    .rng_uniform_at(2137, 0, 5),
    c(0.9498900428996606, 0.019773467515542542, 0.44651496546113223,
      0.48149206798274413, 0.45371203455719133))
  # A 32-bit half of exactly 2^31 reaches the XOR in the first value.
  expect_identical(
    .rng_uniform_at(11895787, 0, 3),
    c(0.0008695978769349688, 0.8180201584014879, 0.6516966238145004))
})

test_that("uniforms equal Python's bit for bit at random seeds and counters", {
  cases <- jsonlite::read_json(test_path("fixtures", "rng", "python_uniforms.json"))
  expect_length(cases, 200L)
  # Python writes 13 hex digits ("0x1.8000000000000p-1"), C's %a drops trailing zeros.
  canon <- function(h) {
    h <- sub("^0x0\\.0+p\\+0$", "0x0p+0", h)
    sub("\\.?0*p", "p", h)
  }
  for (case in cases) {
    expected <- canon(unlist(case$values))
    got <- canon(sprintf("%a", .rng_uniform_at(case$seed, case$start, length(expected))))
    expect_identical(got, expected, label = sprintf("seed %.0f, start %.0f", case$seed, case$start))
  }
})

test_that("draws in order from one generator match Python's", {
  r <- .rng(2137)
  expect_identical(.rng_permutation(r, 12) - 1L, c(1L, 11L, 5L, 2L, 4L, 3L, 9L, 8L, 6L, 7L, 10L, 0L))
  expect_equal(.rng_choice(r, c(0.5, 0, 2, 1.25, 0, 3), 20),
               c(3, 5, 2, 5, 2, 3, 2, 2, 2, 2, 5, 5, 5, 2, 5, 5, 5, 5, 5, 2))
  expect_identical(.rng_integers(r, 1000), 25)
  expect_identical(.rng_integers(r, 7), 1)

  expect_identical(.rng_permutation(.rng(2137), 10) - 1L, c(1L, 5L, 2L, 4L, 3L, 9L, 8L, 6L, 7L, 0L))
  r <- .rng(2137)
  expect_identical(vapply(1:5, function(i) .rng_integers(r, 10), 0), c(9, 0, 4, 4, 4))
})

test_that("chunked draws equal one draw", {
  whole <- .rng_uniform_at(99, 0, 1000)
  parts <- unlist(lapply(seq(0, 750, by = 250), function(s) .rng_uniform_at(99, s, 250)))
  expect_identical(parts, whole)
  expect_true(min(whole) >= 0 && max(whole) < 1)

  a <- .rng(7)
  expect_identical(c(.rng_uniform(a, 3), .rng_uniform(a, 3)), .rng_uniform(.rng(7), 6))
})

test_that("the seed must be a whole number in [0, 2^53)", {
  for (bad in list(-1, 2^53, 1.5, NA, NA_real_, NULL, c(1, 2), "1", TRUE)) {
    expect_error(.rng(bad), "seed must be a whole number", label = deparse(bad))
  }
  expect_silent(.rng(0))
  expect_silent(.rng(2^53 - 1))
  expect_silent(.rng(2137L))
})

test_that("choice never returns an index whose weight is 0", {
  idx <- .rng_choice(.rng(3), c(0, 2, 0, 1, 0), 5000)
  expect_setequal(unique(idx), c(1, 3))
})
