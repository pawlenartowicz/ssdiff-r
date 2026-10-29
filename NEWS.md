# ssdiff 0.1.0

* First release: an R port of the Python package ssdiff 4.0.1, tested
  against reference results written by Python.
* Wordpack embeddings by tag, name or file, with a shared download cache
  (`load_wordpack()`, `wordpack_cache()`, `set_abtt()`).
* Texts from strings or token lists (`ssd()`, `with_lexicon()`) and lexicon
  tools (`suggest_lexicon()`, `evaluate_lexicon()`).
* Four fits: `fit_regression()` (PCA regression with covariates, `ssd()`
  terms, several texts and, through fastglmm, random effects), `fit_pls()`
  and `fit_multipls()` (through plskit), and `fit_groups()`; `contrast()` of
  two gradients.
* Accessors for words, clusters, snippets and documents, `print()`,
  `summary()`, `predict()`, R's model generics, and the offline results page
  (`save_page()`).
