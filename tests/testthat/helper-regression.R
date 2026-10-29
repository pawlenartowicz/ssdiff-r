# The golden text (240 documents) with its data columns, for the
# fit_regression tests: y, cov, group (three labels), a binary yb (y2 > 0), a count
# yc and a grouping factor g of 24 levels. Call inside with_empty_cache().
regression_text <- function() {
  d <- golden_texts()
  emb <- load_wordpack(test_pack())
  list(emb = emb, s = ssd(d$text, emb, lexicon = d$lexicon),
       data = data.frame(y = d$y, cov = d$cov, group = d$group, yb = as.numeric(d$y2 > 0),
                         yc = pmax(0, round(2 + d$y + d$cov)), g = factor(rep(seq_len(24), each = 10)),
                         stringsAsFactors = FALSE))
}
