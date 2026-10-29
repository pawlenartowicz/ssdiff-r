# ssdiff

Supervised Semantic Differential (SSD) for R. SSD finds the semantic
dimension in a set of texts that goes with a numeric outcome or separates
groups. It builds a vector per document from the words around a few seed
words (or from the whole document), fits a gradient in word-embedding space
by PCA regression, PLS or group contrasts, and explains it through the nearest
words of each pole, clusters of those words, text snippets and per-document
scores. This package is an R port of the Python package
[ssdiff](https://github.com/hplisiecki/Supervised-Semantic-Differential),
tested to give its numbers on the same input.

## Installation

```r
install.packages("ssdiff")

# or the development version from r-universe:
install.packages("ssdiff", repos = c("https://pawlenartowicz.r-universe.dev", getOption("repos")))
```

Two backends are on r-universe, not CRAN. Install them once with:

```r
ssdiff::install_backends()
```

| Function | Needs |
|---|---|
| `fit_pls()`, `fit_multipls()`, `fit_regression(ssd_gradients > 1)` | plskit |
| `fit_regression()` with random effects (`(1 \| g)`) or a family base R does not fit the same way | fastglmm |
| everything else | nothing beyond CRAN |

A function that needs a backend says so, with the install command, when it
is called without it.

## Example

```r
library(ssdiff)

emb <- load_wordpack("en")      # downloaded once, then read from the cache
lex <- c("people", "life", "work", "time", "family")
essay <- ssd(df$essay, embeddings = emb, lexicon = lex)

res <- fit_regression(essay, score ~ age + sex * ssd(), data = df)
summary(res)
words(res, unit = "ssd()")      # nearest words of each pole
clusters(res, unit = "ssd()")
coef(res); confint(res); AIC(res)

r1 <- fit_pls(essay, "score", data = df)   # needs plskit
r2 <- fit_pls(essay, "age", data = df)
cc <- contrast(r1$gradient, r2$gradient, emb)

save_page(res, "results.html")  # the offline results page, as in Python
```

`fit_groups(essay, "group", data = df)` compares groups, and
`fit_regression(score ~ age + ssd(essay), data = df)` builds the text from
the column `essay` with every default.

## The pack cache

A wordpack is about 380 MB. It is downloaded once into
`tools::R_user_dir("wordpack", "cache")`, or into the folder in the
environment variable `WORDPACK_CACHE` when it is set. Set
`WORDPACK_CACHE=~/.wordpacks` to share packs with Python ssdiff.
`wordpack_cache()` lists the packs there and their sizes; delete old ones by
hand.

## What R does not do

R loads wordpacks only (no `.txt`, `.vec`, `.bin` or `.kv` embedding files)
and has no profile mode, view classes, table exporters or `report()`; tables
are data frames and `save_page()` writes the results page.

## Citation

Plisiecki, Lenartowicz, Pokropek, Małyska and Flakus (2025).
<https://doi.org/10.31234/osf.io/gvrsb_v1>
