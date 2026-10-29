# fit_regression() with several texts, and texts built from columns of data
# (spec 03, "Texts from columns", "Rows and texts").

# The golden text as `essay`, and its token lists with every 20th document from
# the 8th emptied as `diary` (12 documents without a vector).
two_texts <- function() {
  r <- regression_text()
  docs <- r$s$docs
  gaps <- seq(8L, length(docs), by = 20L)
  docs[gaps] <- list(character(0))
  r$diary <- ssd(docs, r$emb)
  r$gaps <- gaps
  r
}

test_that("several texts: rows without a vector in any text drop, each text has its own K", with_empty_cache({
  r <- two_texts()
  res <- fit_regression(list(essay = r$s, diary = r$diary), y ~ cov + ssd(essay) + group * ssd(diary),
                        data = r$data, fixed_k = 4)
  expect_identical(res$stats$n_kept, 228L)
  expect_identical(res$texts[[1]]$keep, !seq_len(240) %in% r$gaps)
  expect_identical(vapply(res$texts, `[[`, 0L, "n_no_vector"), c(0L, 12L))
  expect_identical(vapply(res$texts, `[[`, "", "name"), c("essay", "diary"))
  expect_identical(vapply(res$texts, `[[`, 0L, "k"), c(4L, 4L))
  expect_identical(names(res$units), c("ssd(essay)", "ssd(diary)", "grouplow:ssd(diary)", "groupmid:ssd(diary)"))
  expect_identical(vapply(res$units, `[[`, 0L, "text"), c(1L, 2L, 2L, 2L), ignore_attr = TRUE)
  expect_identical(unname(vapply(res$units, `[[`, "", "test")), c(NA, "ssd(diary)", "group:ssd(diary)", "group:ssd(diary)"))
  expect_true(all(c(".ssd_essay_pca_1", ".ssd_diary_pca_4") %in% names(coef(res$model))))
  expect_identical(res$test$label, c(NA, "essay", "diary", "ssd(diary)", "group:ssd(diary)"))

  # Each text's test drops that text's terms only.
  d <- res$model$model
  essay <- paste0(".ssd_essay_pca_", 1:4)
  diary <- paste0(".ssd_diary_pca_", 1:4)
  m0 <- lm(reformulate(c("cov", "group", diary, paste0("group:", diary)), ".ssd_y"), data = d)
  expect_equal(res$test$pvalue[2], anova(m0, res$model)[2L, "Pr(>F)"], tolerance = 1e-10)
  m0 <- lm(reformulate(c("cov", "group", essay), ".ssd_y"), data = d)
  expect_equal(res$test$pvalue[3], anova(m0, res$model)[2L, "Pr(>F)"], tolerance = 1e-10)
  expect_equal(res$test$df[1:3], c(16, 4, 12))
}))

test_that("several texts with a sweep: one K per text", with_empty_cache({
  r <- two_texts()
  res <- fit_regression(list(essay = r$s, diary = r$diary), y ~ ssd(essay) + ssd(diary), data = r$data,
                        k_min = 2, k_max = 8, k_step = 2)
  for (t in res$texts) {
    expect_identical(t$sweep$k, c(2L, 4L, 6L, 8L))
    expect_identical(t$k, t$sweep$k[which.max(t$sweep$joint_score)])
  }
  expect_identical(res$units[["ssd(diary)"]]$fit_info$n_components, res$texts[[2]]$k)
}))

test_that("texts from columns: built with every default, with one pack load for all", with_empty_cache({
  r <- two_texts()
  d <- golden_texts()
  data <- cbind(r$data, essay = d$text, diary = rev(d$text), stringsAsFactors = FALSE)
  loads <- 0L
  local_mocked_bindings(load_wordpack = function(...) {
    loads <<- loads + 1L
    r$emb
  })
  res <- fit_regression(y ~ ssd(essay) + ssd(diary), data = data, fixed_k = 3)
  expect_identical(loads, 1L)
  expect_identical(res$texts[[1]]$use_full_doc, TRUE)
  expect_identical(res$texts[[1]]$x, ssd(d$text, r$emb)$x)

  # The formula first is the empty list first.
  same <- fit_regression(list(), y ~ ssd(essay) + ssd(diary), data = data, fixed_k = 3)
  expect_identical(same$test, res$test)

  # An element of s wins over the column of the same name.
  res <- fit_regression(list(essay = r$s), y ~ ssd(essay) + ssd(diary), data = data, fixed_k = 3)
  expect_identical(res$texts[[1]]$lexicon, r$s$lexicon)
  expect_identical(res$texts[[2]]$use_full_doc, TRUE)

  # One text from a column has no text level: its components are .ssd_pca_*.
  res <- fit_regression(y ~ ssd(essay), data = data, fixed_k = 3)
  expect_identical(res$terms$term, "ssd(essay)")
  expect_true(".ssd_pca_1" %in% names(coef(res$model)))
}))

test_that("texts that are neither elements nor columns, or not named, raise", with_empty_cache({
  r <- two_texts()
  expect_error(fit_regression(y ~ ssd(essay), data = r$data),
               "the formula names the texts 'essay', which are neither columns of data nor elements of s", fixed = TRUE)
  expect_error(fit_regression(list(essay = r$s, diary = r$diary), y ~ ssd(essay), data = r$data),
               "the texts 'diary' were passed but the formula names only 'essay'", fixed = TRUE)
  expect_error(fit_regression(list(essay = r$s), y ~ ssd(), data = r$data), "bare ssd() needs one text", fixed = TRUE)
  expect_error(fit_regression(list(essay = r$s), y ~ cov, data = r$data), "has no ssd(name) term", fixed = TRUE)
  expect_error(fit_regression(list(essay = r$s), y ~ ssd(essay)), "fit_regression() needs data = df.", fixed = TRUE)
  expect_error(fit_regression(list(essay = r$s, diary = with_lexicon(r$s, "life")), y ~ ssd(essay) + ssd(diary),
                              data = r$data[1:10, ]), "nrow(data)=10", fixed = TRUE)
  short <- ssd(r$s$docs[1:10], r$emb)
  expect_error(fit_regression(list(essay = r$s, diary = short), y ~ ssd(essay) + ssd(diary), data = r$data),
               "the texts differ in length: essay = 240, diary = 10.", fixed = TRUE)
}))

test_that("a dotted text name works on the fastglmm path", with_empty_cache({
  skip_if_not_installed("fastglmm")
  r <- two_texts()
  res <- suppressWarnings(fit_regression(list(my.essay = r$s, diary = r$diary), y ~ ssd(my.essay) + ssd(diary) + (1 | g),
                                         data = r$data, fixed_k = 3))
  expect_identical(res$terms$term, c("ssd(my.essay)", "ssd(diary)"))
  expect_true(".ssd_my.essay_pca_1" %in% names(.reg_coef(res$model)))
  expect_true(all(startsWith(setdiff(names(fastglmm::fixef(res$model)), "(Intercept)"), "ssd__")))
  expect_identical(names(coef(res)), "(Intercept)")
}))
