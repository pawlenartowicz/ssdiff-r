# Port of ssdiff/utils/rng.py. Change the two together.

# Counter-based generator: value i (i = 0, 1, ...) is splitmix64's
# mix(seed + (i + 1) * 0x9E3779B97F4A7C15 mod 2^64). Mirrors
# ssdiff/utils/rng.py in the Python package; change the two together.
#
# R has no 64-bit integer, so a uint64 is a pair of doubles (hi, lo), each a
# whole number in [0, 2^32), and products are taken over 16-bit limbs. Every
# intermediate value stays below 2^53, so all arithmetic is exact. Quotients
# use floor(x / 2^k), which is exact and much faster than %/% and %% on doubles.

# A 64-bit constant written in hex, as four 16-bit limbs, lowest first.
.u64_limbs <- function(hex) {
  rev(strtoi(substring(hex, c(1, 5, 9, 13), c(4, 8, 12, 16)), 16L))
}
.rng_gamma <- .u64_limbs("9E3779B97F4A7C15")
.rng_m1 <- .u64_limbs("BF58476D1CE4E5B9")
.rng_m2 <- .u64_limbs("94D049BB133111EB")

# (hi, lo) * c mod 2^64, c as four limbs.
.u64_mul <- function(x, c) {
  lo1 <- floor(x$lo / 65536)
  hi1 <- floor(x$hi / 65536)
  a <- list(x$lo - lo1 * 65536, lo1, x$hi - hi1 * 65536, hi1)
  out <- vector("list", 4L)
  carry <- 0
  for (k in 1:4) {
    t <- carry
    for (i in 1:k) t <- t + a[[i]] * c[k - i + 1L]
    carry <- floor(t / 65536)
    out[[k]] <- t - carry * 65536
  }
  list(hi = out[[3]] + out[[4]] * 65536, lo = out[[1]] + out[[2]] * 65536)
}

# bitwXor takes 32-bit signed integers, and R has no -2^31 (it is NA_integer_):
# XOR the low 31 bits and put bit 31 back as a double.
.u32_xor <- function(a, b) {
  ha <- a >= 2^31
  hb <- b >= 2^31
  as.double(bitwXor(as.integer(a - ha * 2^31), as.integer(b - hb * 2^31))) + (ha != hb) * 2^31
}

# x ^ (x >> s), for 0 < s < 32.
.u64_xorshift <- function(x, s) {
  q <- floor(x$hi / 2^s)
  list(hi = .u32_xor(x$hi, q),
       lo = .u32_xor(x$lo, floor(x$lo / 2^s) + (x$hi - q * 2^s) * 2^(32 - s)))
}

# Raw outputs at counter positions start .. start + n - 1, as (hi, lo).
.rng_raw_at <- function(seed, start, n) {
  i <- start + seq_len(n)                       # i + 1 for i = start, ...
  ihi <- floor(i / 2^32)
  z <- .u64_mul(list(hi = ihi, lo = i - ihi * 2^32), .rng_gamma)
  shi <- floor(seed / 2^32)
  lo <- z$lo + (seed - shi * 2^32)
  carry <- lo >= 2^32
  hi <- z$hi + shi + carry
  z <- list(hi = hi - (hi >= 2^32) * 2^32, lo = lo - carry * 2^32)
  z <- .u64_mul(.u64_xorshift(z, 30), .rng_m1)
  z <- .u64_mul(.u64_xorshift(z, 27), .rng_m2)
  .u64_xorshift(z, 31)
}

# Uniform doubles in [0, 1): (raw >> 11) * 2^-53.
.rng_uniform_at <- function(seed, start, n) {
  z <- .rng_raw_at(seed, start, n)
  (z$hi * 2^21 + floor(z$lo / 2^11)) * 2^-53
}

# A generator: the seed and a counter that every draw advances.
.rng <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      seed != floor(seed) || seed < 0 || seed >= 2^53) {
    stop("seed must be a whole number with 0 <= seed < 2^53, got ", format(seed), call. = FALSE)
  }
  r <- new.env(parent = emptyenv())
  r$seed <- as.double(seed)
  r$counter <- 0
  r
}

.rng_uniform <- function(r, n) {
  u <- .rng_uniform_at(r$seed, r$counter, n)
  r$counter <- r$counter + n
  u
}

# One integer in [0, n): floor(u * n). 0-based, as in Python.
.rng_integers <- function(r, n) floor(.rng_uniform(r, 1) * n)

# A permutation of 1..n: the stable order of n uniforms (numpy's stable argsort, plus 1).
.rng_permutation <- function(r, n) order(.rng_uniform(r, n), method = "radix")

# size weighted draws with replacement, 0-based indices into w.
# findInterval(x, c) equals numpy's searchsorted(c, x, side = "right").
.rng_choice <- function(r, w, size) {
  c <- cumsum(w)
  u <- .rng_uniform(r, size)
  pmin(findInterval(u * c[length(c)], c), length(w) - 1L)
}
