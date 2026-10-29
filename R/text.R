# Port of the wordpack route of ssdiff/utils/text.py and the text handling of SSD.__init__ in ssdiff/ssd.py. Change them together.

# Tokens of a text in one of the two input forms (02 "Texts"): `docs` (tokens
# per document), `pre_docs` (sentence structure for snippets, NULL for token
# lists) and `tokenizer` (the record Python keeps, NULL for token lists).
.tokenize_texts <- function(texts, emb) {
  if (.texts_form(texts) == "strings") {
    .check_embeddings(emb, "embeddings")
    texts <- unname(texts)
    texts[is.na(texts)] <- ""
    pre_docs <- .preprocess_texts_wordpack(texts, emb$pack, emb$stopwords)
    return(list(docs = lapply(pre_docs, `[[`, "doc_lemmas"), pre_docs = pre_docs,
                tokenizer = paste0("wordpack:", .full_name(emb))))
  }
  docs <- lapply(unname(texts), function(d) if (is.null(d)) character(0) else unname(d))
  list(docs = docs, pre_docs = NULL, tokenizer = NULL)
}

# "strings" or "tokens", the form of `texts`; stops on any other input.
.texts_form <- function(texts) {
  if (is.character(texts)) return("strings")
  if (is.list(texts) && !is.object(texts) &&
      all(vapply(texts, function(d) is.null(d) || is.character(d), TRUE))) {
    return("tokens")
  }
  stop("texts must be a character vector (one text per document) or a list of character vectors ",
       "(one token list per document)", call. = FALSE)
}

# Per text, Python's PreprocessedDoc fields. A token is dropped when its key,
# lowercased by the pack's rule, is a stopword. Sentence and token positions
# count from 1; sent_char_spans are inclusive and index `raw`.
.preprocess_texts_wordpack <- function(texts, pack, stopwords) {
  tok <- wordpack_tokenize(pack, texts)
  n_docs <- length(texts)
  sents <- tok$sents
  t <- tok$tokens

  keys <- unique(t$key)
  lowered <- vapply(keys, function(k) intToUtf8(.wp_map(utf8ToInt(k), pack$chars$lower)), "", USE.NAMES = FALSE)
  kept <- !(lowered %in% stopwords)[match(t$key, keys)]

  # Each sentence gets one id across all texts; sents are ordered by text, then sentence.
  n_sents <- tabulate(sents$doc, n_docs)
  first_sent <- cumsum(c(0L, n_sents))[seq_len(n_docs)]
  sent_id <- first_sent[t$doc] + t$sent
  # Position of each token within its sentence, stopwords counted.
  pos <- seq_along(sent_id) - match(sent_id, sent_id) + 1L
  by_sent <- factor(sent_id[kept], levels = seq_len(nrow(sents)))
  sent_keys <- unname(split(t$key[kept], by_sent))
  sent_pos <- unname(split(pos[kept], by_sent))
  doc_of <- factor(t$doc[kept], levels = seq_len(n_docs))
  doc_keys <- unname(split(t$key[kept], doc_of))
  doc_sent <- unname(split(t$sent[kept], doc_of))

  lapply(seq_len(n_docs), function(d) {
    s <- first_sent[d] + seq_len(n_sents[d])
    raw <- tok$text[d]
    list(
      raw = raw,
      sents_surface = if (length(s)) substring(raw, sents$start[s], sents$end[s]) else character(0),
      sents_lemmas = sent_keys[s],
      doc_lemmas = doc_keys[[d]],
      sent_char_spans = cbind(start = sents$start[s], end = sents$end[s]),
      token_to_sent = doc_sent[[d]],
      sents_kept_idx = sent_pos[s]
    )
  })
}
