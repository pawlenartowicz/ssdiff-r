# Port of ssdiff/results/contrast_result.py. Change the two together.

#' The difference of two gradients
#'
#' `contrast()` reads the direction in which two gradients differ:
#' `unit(g1) - unit(g2)`, scaled to unit length, so that neither the
#' outcomes' scales nor the backends enter the difference. The gradients can
#' come from any fits on embeddings of the same dimension, for example
#' `r1$gradient` and `r2$gradient`, or two units,
#' `res$units[["sexM:ssd(essay)"]]$gradient`.
#'
#' The result has words and clusters ([words()], [clusters()],
#' [cluster_words()]) in `emb`, and no documents: [snippets()],
#' [cluster_snippets()], [docs()] and [doc_text()] raise. There is no test;
#' `cos_g1_g2` is the only statement about how the two gradients relate.
#'
#' @param g1,g2 Two numeric vectors of the embeddings' dimension.
#' @param emb Embeddings from [load_wordpack()], read for the words.
#' @param labels The names of the two gradients.
#' @return An object of class `c("ssd_contrast", "ssd_result")`, a list with
#'   `beta` and `gradient` (both the unit difference), `stats` (a one-row
#'   data frame: `cos_g1_g2`, the cosine of the two gradients, and
#'   `diff_norm`, the length of `unit(g1) - unit(g2)` before scaling),
#'   `labels`, and `texts` (one entry holding the pack of `emb`).
#' @examples
#' \dontrun{
#' r1 <- fit_pls(s, "score", data = df)
#' r2 <- fit_pls(s, "age", data = df)
#' cc <- contrast(r1$gradient, r2$gradient, emb, labels = c("score", "age"))
#' words(cc)
#' }
#' @export
contrast <- function(g1, g2, emb, labels = c("g1", "g2")) {
  .check_embeddings(emb)
  if (!is.numeric(g1) || !is.numeric(g2)) stop("g1 and g2 must be numeric vectors", call. = FALSE)
  g1 <- as.double(g1)
  g2 <- as.double(g2)
  if (length(g1) != length(g2)) {
    stop(sprintf("g1 and g2 must have the same length, got %d and %d.", length(g1), length(g2)), call. = FALSE)
  }
  if (length(g1) != emb$pack$dims) {
    stop(sprintf("g1 and g2 have length %d, but the embeddings have %d dimensions.", length(g1), emb$pack$dims),
         call. = FALSE)
  }
  if (length(labels) != 2L) stop("labels must hold two names, one per gradient.", call. = FALSE)
  u1 <- .unit_vector(g1)
  u2 <- .unit_vector(g2)
  diff <- u1 - u2
  diff_norm <- sqrt(sum(diff^2))
  if (!is.finite(diff_norm) || diff_norm < 1e-12) {
    stop("g1 and g2 are the same direction (or zero): the contrast is empty.", call. = FALSE)
  }
  beta <- diff / diff_norm
  structure(list(
    beta = beta,
    gradient = .unit_vector(beta),
    stats = data.frame(cos_g1_g2 = sum(u1 * u2), diff_norm = diff_norm),
    labels = as.character(labels),
    texts = list(list(name = NULL, docs = NULL, pre_docs = NULL, tokenizer = NULL, pack = .pack_ref(emb),
                      lexicon = NULL, window = NULL, sif_a = NULL, use_full_doc = NULL, x = NULL, keep = NULL,
                      n_no_vector = NULL)),
    fit_args = NULL
  ), class = c("ssd_contrast", "ssd_result"))
}
