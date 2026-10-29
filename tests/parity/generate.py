"""Write the files the R tests compare against: Python ssdiff's results.

``run.py`` runs this in a virtual environment holding the ssdiff release
pinned in ``parity.toml``; CI and R CMD check never run it. Both parts, or
the one named:
    python tests/parity/generate.py [golden|reference]

``golden``: committed files, on the test pack (the committed
``tests/testthat/fixtures/wordpack/test.wordpack``). After a run, a changed
file in ``git status`` means Python's results moved.
  - ``fixtures/golden/``: ``texts.json``, the golden cases (``GOLDEN``,
    ``two_texts_case``, ``contrast_case``, ``lexicon_cases``), the document
    vectors of their texts, and ``versions.json``.
  - ``fixtures/page/<kind>.json``: ``page_data()`` of one fit per result kind.
  - ``fixtures/rng/python_uniforms.json``: uniforms of ``ssdiff/utils/rng.py``.
  - ``fixtures/lexicon/python_small.json``: lexicon results on a 12-document text.

``reference``: every case (``case_list``) on the real pack (``pack`` in
parity.toml, fetched by name into Python's pack cache), written into the
git-ignored ``out/`` with its ``texts.json`` and ``pack.json`` (the pack's name
and the file Python read, which R loads too).

Case files: ``<case>.json``, one per fit, ``contrast()`` or lexicon call; each
names its texts under ``texts``. ``vectors_<text>.json``: every document vector
of one text configuration (``TEXT_CONFIGS``), with Python's 0-based ``doc_id``;
floats are written so they read back to the same double.

``texts.json``: synthetic texts drawn from the pack's vocabulary with a planted
signal, the columns ``text``, ``y``, ``y2`` (``Y`` is both), ``group`` (three
labels) and ``cov``, and the ``lexicon``. Two latent outcomes z1, z2 ~ N(0, 1)
per document. Each context word ``w`` is drawn with probability proportional to
``exp(TILT * (z1 * p1(w) + z2 * p2(w)))``, where p1 and p2 are the word's
standardized projections on the pack's two pole directions (``*_TEXTS``).
``y = z1``, ``y2 = z2``.
"""

from __future__ import annotations

import json
import re
import sys
import tomllib
from importlib.metadata import version
from pathlib import Path

import numpy as np

import ssdiff
from ssdiff import SSD, Embeddings
from ssdiff.results.continuous_result import ClustersView, ClusterWordsView
from ssdiff.results.page import _clean, _view_table, page_data
from ssdiff.utils.rng import uniform_at
from ssdiff.wordpack.fetch import resolve

HERE = Path(__file__).parent
FIXTURES = HERE.parent / "testthat" / "fixtures"
GOLDEN_DIR = FIXTURES / "golden"
OUT = HERE / "out"
TEST_PACK = FIXTURES / "wordpack" / "test.wordpack"
PINS = tomllib.loads((HERE / "parity.toml").read_text(encoding="utf-8"))

# --- Texts -----------------------------------------------------------------
SEED = 2137
N_DOCS = 240
LEXICON = ["people", "like"]  # words of both packs, a stopword in neither
# Per pack: the frequency ranks the context words come from, and the poles of
# the two planted directions.
REFERENCE_TEXTS = {"ranks": (100, 4000), "poles": (("good", "bad"), ("money", "love"))}
GOLDEN_TEXTS = {"ranks": (0, None), "poles": (("dog", "cat"), ("blue", "green"))}
# Weak on purpose: a strong signal puts every permutation p-value at its floor
# 1 / (n_perm + 1), which any permutation sequence reproduces, so R's generator
# would go untested.
TILT = 0.1
N_TOP_WORDS = 20
N_VECTOR_ROWS = 5
N_TOKEN_TEXTS = 3

# Clustering arguments of the cluster tables, written into each case as
# ``cluster_args``; Python's defaults otherwise. Rounding decides the order of
# words whose cosines tie, which R cannot reproduce: a two-word cluster's words
# are equally close to their centroid, so every case leaves such clusters out.
CLUSTER_ARGS = {"min_cluster_size": 3}

# How each text is built from the ``text`` column of texts.json: as strings, or
# as token lists (``edit_tokens``: the tokenizer's lists as they are, ``raw`` or
# ``gaps``); with the lexicon, or from full documents. Each has its own vectors file.
TEXT_CONFIGS = {
    "text": {"token_lists": False, "edit": None, "lexicon": LEXICON},
    "tokens": {"token_lists": True, "edit": None, "lexicon": LEXICON},
    "tokens_raw": {"token_lists": True, "edit": "raw", "lexicon": LEXICON},
    "tokens_gaps": {"token_lists": True, "edit": "gaps", "lexicon": None},
}
# ``raw``: stopwords that are pack words, which the string path drops.
RAW_STOPWORDS = ("the", "and", "of")
# ``gaps``: documents 7, 27, 47, ... are emptied, so they have no vector.
GAP_START, GAP_EVERY = 7, 20


def sweep(emb) -> dict:
    """The K sweep of fit_regression: Python's default k_max (120), stopped below
    the rank of the document vectors. The pack's ABTT takes one direction out of
    every word vector, so they span dims - 1 dimensions, and at K >= dims - 1 the
    fit's last component is rounding noise in both languages."""
    return {"k_max": min(120, emb.vector_size - 2)}


def case_list(sweep: dict) -> list:
    """Each case: (name, text configuration, ssd y, call, positional args, keyword
    args). Every fit kind runs with default arguments plus one non-default case."""
    return [
        ("pls_default", "text", "y", "fit_pls", [], {}),
        ("pls_k2", "text", "y", "fit_pls", [], {"k": 2}),
        ("pls_gradients3", "text", "y", "fit_pls", [], {"ssd_gradients": 3}),
        ("pls_formula", "text", "y", "fit_pls", ["~ cov"], {}),
        ("multipls_default", "text", ["y", "y2"], "fit_multipls", [], {}),
        ("multipls_rotate_outcomes", "text", ["y", "y2"], "fit_multipls", [], {"rotate": "outcomes"}),
        ("regression_default", "text", "y", "fit_regression", [], sweep),
        ("regression_fixed_k10", "text", "y", "fit_regression", [], {"fixed_k": 10}),
        ("regression_formula", "text", "y", "fit_regression", ["~ cov"], sweep),
        ("groups_default", "text", "group", "fit_groups", [], {}),
        ("groups_median_split", "text", "y", "fit_groups", [], {"median_split": True}),
        # The same fit as pls_default on the token lists of the same texts.
        ("tokens_pls_default", "tokens", "y", "fit_pls", [], {}),
        # Fixed K: at a large swept K the three text terms leave few residual df.
        ("regression_interaction", "text", "y", "fit_regression", ["~ group*ssd(text)"], {"fixed_k": 10}),
        ("regression_gradients3", "text", "y", "fit_regression", [], {"ssd_gradients": 3, **sweep}),
        # Token lists used as given: upper-case forms and stopwords stay in.
        ("tokens_raw_regression_fixed_k10", "tokens_raw", "y", "fit_regression", [], {"fixed_k": 10}),
    ]


# The cases of ``case_list`` on the test pack: one per fit kind.
GOLDEN = ("pls_default", "multipls_default", "regression_default", "groups_default")
# ssdiff.fit_regression with two texts: name in the formula -> text configuration.
# The diary has documents without a vector, so the two texts keep different rows.
TWO_TEXTS = {"essay": "text", "diary": "tokens_gaps"}
TWO_TEXTS_FORMULA = "y ~ cov + ssd(essay) + ssd(diary)"

# --- Small lexicon text ----------------------------------------------------
# 12 token-list documents on the test pack: quantile bins with ties, 3 and 2
# groups (Yates' correction), seeds with and without an outcome.
SMALL_DOCS = [
    ["dog", "cat", "swim"], ["dog", "dogs", "blue"], ["cat", "duck", "green"], ["dog", "people", "like"],
    ["cat", "swim", "never"], ["duck", "blue", "happens"], ["dog", "cat", "green"], ["people", "like", "they"],
    ["dogs", "duck", "swim"], ["cat", "blue", "never"], ["dog", "green", "with"], ["duck", "people", "are"],
]
SMALL_Y = [0.5, 1.2, -0.3, 2.0, 0.1, -1.0, 0.7, 1.5, -0.2, 0.0, 0.9, -0.8]
SMALL_G3 = ["a", "b", "c"] * 4
SMALL_G2 = ["x", "y"] * 6

# --- Uniforms --------------------------------------------------------------
N_SEEDS, N_VALUES = 200, 50


def make_texts(emb, ranks, poles) -> dict:
    rng = np.random.default_rng(SEED)
    stop = set(emb.stopwords)
    lo, hi = ranks
    vocab = [w for w in emb.index_to_key[lo:hi]
             if re.fullmatch(r"[a-z]+", w) and w not in stop and w not in LEXICON]
    V = np.vstack([emb[w] for w in vocab]).astype(np.float64)

    def projection(a, b):
        p = V @ (emb[a] - emb[b])
        return (p - p.mean()) / p.std()

    p1, p2 = (projection(a, b) for a, b in poles)
    z = rng.standard_normal((N_DOCS, 2))
    texts = []
    for z1, z2 in z:
        logits = TILT * (z1 * p1 + z2 * p2)
        prob = np.exp(logits - logits.max())
        prob /= prob.sum()
        sentences = []
        for _ in range(rng.integers(3, 7)):
            sent = list(rng.choice(vocab, size=rng.integers(6, 11), p=prob))
            if rng.random() < 0.8:
                sent.insert(int(rng.integers(0, len(sent) + 1)), str(rng.choice(LEXICON)))
            sentences.append(" ".join(sent).capitalize() + ".")
        texts.append(" ".join(sentences))
    y, y2 = z[:, 0], z[:, 1]
    cut = np.quantile(y + 0.5 * rng.standard_normal(N_DOCS), [1 / 3, 2 / 3])
    group = np.where(y < cut[0], "low", np.where(y < cut[1], "mid", "high"))
    cov = 0.3 * y + rng.standard_normal(N_DOCS)
    return {"lexicon": LEXICON, "text": texts, "y": y.tolist(), "y2": y2.tolist(),
            "group": group.tolist(), "cov": cov.tolist()}


def _columns(data: dict) -> dict:
    """texts.json as the ``data`` of a fit."""
    return {k: np.asarray(v) for k, v in data.items() if k != "lexicon"}


def text_column(data: dict, emb, lexicon) -> SSD:
    """The ``text`` column as a text, named after the column as Python names a DataFrame's."""
    ssd = SSD(data["text"], embeddings=emb, lexicon=lexicon)
    ssd.name = "text"
    return ssd


def edit_tokens(docs: list[list[str]], edit: str | None, emb) -> list[list[str]]:
    """Token lists made from the tokenizer's *docs*.

    ``None``: as they are. ``raw``: forms the string path never yields, which a
    token-list text keeps: the first token upper-cased (even documents), the
    first lexicon word capitalised (every third document; that seed no longer
    matches), and a ``RAW_STOPWORDS`` word inserted at position 2 (odd
    documents). ``gaps``: every ``GAP_EVERY``-th document from ``GAP_START``
    emptied.
    """
    out = [list(d) for d in docs]
    if edit == "gaps":
        for i in range(GAP_START, len(out), GAP_EVERY):
            out[i] = []
    elif edit == "raw":
        stop = set(emb.stopwords)
        assert all(w in stop and w in emb for w in RAW_STOPWORDS), "RAW_STOPWORDS must be pack stopwords"
        for i, d in enumerate(out):
            if i % 2 == 0 and d:
                d[0] = d[0].upper()
            if i % 3 == 0:
                j = next((j for j, t in enumerate(d) if t in LEXICON), None)
                if j is not None:
                    d[j] = d[j].capitalize()
            if i % 2 == 1:
                d.insert(min(2, len(d)), RAW_STOPWORDS[i % 3])
        assert not any(t in emb for d in out for t in d if t != t.lower()), "an upper-case form is a pack key"
    elif edit is not None:
        raise ValueError(f"unknown token edit {edit!r}")
    return out


def _table_or_none(view):
    return None if view is None else _view_table(view)


def _unit(leaf, cluster_args: dict, clusters: bool = False) -> dict:
    pos, neg = leaf.clusters.pos(**cluster_args), leaf.clusters.neg(**cluster_args)
    unit = {
        "beta": _clean(np.asarray(leaf.beta)),
        "words": {"pos": _view_table(leaf.words.pos(N_TOP_WORDS)),
                  "neg": _view_table(leaf.words.neg(N_TOP_WORDS))},
        "cluster_words": _view_table(ClusterWordsView(list(pos._words_rows) + list(neg._words_rows))),
    }
    if clusters:
        rows = list(pos._rows) + list(neg._rows)
        unit["clusters"] = _view_table(ClustersView(leaf, _rows_override=rows), drop=("top_snippet",))
    return unit


def _key(key) -> str:
    return key if isinstance(key, str) else "_".join(map(str, key))


def _units(result, cluster_args: dict) -> dict:
    """The main fit (when ``beta`` is one vector) and every leaf: gradients, dimensions, pairs,
    terms or texts. A text with several terms nests its terms."""
    units = {} if isinstance(result.beta, dict) else {"main": _unit(result, cluster_args)}
    for key in result.keys() if hasattr(result, "keys") else []:
        leaf = result[key]
        units[_key(key)] = _units(leaf, cluster_args) if isinstance(leaf.beta, dict) else _unit(leaf, cluster_args)
    return units


def _results(result) -> dict:
    """The result's tables, keys and units. With several texts, ``fit_info`` and
    ``sweep`` are keyed by text name."""
    out = {}
    for table in ("stats", "test", "fit_info", "pairs", "pls_info", "sweep"):
        view = getattr(result, table, None)
        if isinstance(view, dict):
            out[table] = {k: _table_or_none(v) for k, v in view.items()}
        elif view is not None:
            out[table] = _view_table(view)
    if getattr(result, "_whole_test", None) is not None:
        # fit_regression: one Wald row per result key (labelled by term, as
        # ``pvalues``), the whole-model test, the covariates, and the rows kept.
        cols = ("label", "name", "pvalue", "statistic", "df", "df_resid")
        out["terms"] = {"columns": list(cols),
                        "rows": [[_clean(getattr(t, c)) for c in cols] for t in result._term_tests]}
        out["model_test"] = {k: _clean(getattr(result, k)) for k in ("fvalue", "f_pvalue", "llr", "llr_pvalue")}
        out["covariates"] = _view_table(result.covariates)
        out["n_no_vector"] = _clean(result.n_no_vector)
        out["kept_doc_id"] = _clean(np.asarray(result._obs_index))
    out["keys"] = [_key(k) for k in result.keys()] if hasattr(result, "keys") else []
    out["cluster_args"] = CLUSTER_ARGS
    out["units"] = _units(result, CLUSTER_ARGS)
    return out


def _text_info(config: str, ssd) -> dict:
    """How the text is built and which file holds its document vectors."""
    return {
        "config": config,
        "column": "text",
        "token_lists": ssd.tokenizer is None,
        "lexicon": sorted(ssd.lexicon) or None,
        "window": ssd.window,
        "sif_a": ssd.sif_a,
        "use_full_doc": ssd.use_full_doc,
        "vectors": f"vectors_{config}.json",
    }


def _first_vectors(ssd) -> dict:
    x, has_vector, _ = ssd._vectors()
    return {"doc_id": np.flatnonzero(has_vector)[:N_VECTOR_ROWS].tolist(), "x": _clean(x[:N_VECTOR_ROWS])}


def fit_case(texts, df, name, config, y, call, args, kwargs) -> dict:
    ssd = texts[config]
    outcome = f"{y} {args[0]}" if args else y
    result = getattr(ssd, call)(outcome, data=df, **kwargs)
    out = {
        "case": name,
        "ssd": {"y": y, "lexicon": TEXT_CONFIGS[config]["lexicon"]},
        "texts": {config: _text_info(config, ssd)},
        "call": call,
        "args": list(args),
        "kwargs": kwargs,
        "tokens": ssd.docs[:N_TOKEN_TEXTS],
        "n_kept": ssd.n_kept,
        "doc_vectors": _first_vectors(ssd),
    }
    if TEXT_CONFIGS[config]["token_lists"]:
        out["token_lists"] = ssd.docs
    return out | _results(result)


def two_texts_case(texts, df, sweep: dict) -> dict:
    """``ssdiff.fit_regression(formula, data, essay=..., diary=...)``: one model, two texts."""
    passed = {n: texts[c] for n, c in TWO_TEXTS.items()}
    result = ssdiff.fit_regression(TWO_TEXTS_FORMULA, df, **sweep, **passed)
    out = {
        "case": "regression_two_texts",
        "texts": {n: _text_info(c, passed[n]) for n, c in TWO_TEXTS.items()},
        "call": "ssdiff.fit_regression",
        "args": [TWO_TEXTS_FORMULA],
        "kwargs": sweep,
        "tokens": {n: t.docs[:N_TOKEN_TEXTS] for n, t in passed.items()},
        "n_kept": {n: t.n_kept for n, t in passed.items()},
        "doc_vectors": {n: _first_vectors(t) for n, t in passed.items()},
        "token_lists": {n: t.docs for n, t in passed.items() if t.tokenizer is None},
    }
    return out | _results(result)


def contrast_case(texts, df) -> dict:
    """``ssdiff.contrast`` of the gradients of ``fit_pls`` on ``y`` and on ``y2``."""
    ssd = texts["text"]
    fits = [{"call": "fit_pls", "args": [y], "kwargs": {}} for y in ("y", "y2")]
    g1, g2 = (ssd.fit_pls(f["args"][0], data=df).gradient for f in fits)
    labels = ("y", "y2")
    c = ssdiff.contrast(g1, g2, ssd.embeddings, labels=labels)
    return {
        "case": "contrast",
        "texts": {"text": _text_info("text", ssd)},
        "call": "ssdiff.contrast",
        "args": ["g1", "g2", "embeddings"],
        "kwargs": {"labels": list(labels)},
        "g1": {**fits[0], "gradient": _clean(g1)},
        "g2": {**fits[1], "gradient": _clean(g2)},
        "stats": _view_table(c.stats),
        "keys": [],
        "cluster_args": CLUSTER_ARGS,
        "units": {"main": _unit(c, CLUSTER_ARGS, clusters=True)},
    }


def _lexicon_tables(result) -> dict:
    return {"stats": _view_table(result.stats), "suggestions": _view_table(result.suggestions),
            "summary": _table_or_none(result.summary)}


def lexicon_cases(ssd, df) -> list[dict]:
    y = np.asarray(df["y"])
    return [
        {"case": "suggest_lexicon", "call": "suggest_lexicon", "args": ["y"], "kwargs": {},
         **_lexicon_tables(ssd.suggest_lexicon(y))},
        {"case": "evaluate_lexicon", "call": "evaluate_lexicon", "args": ["y", LEXICON], "kwargs": {},
         **_lexicon_tables(ssd.evaluate_lexicon(y, LEXICON))},
    ]


def _write(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, allow_nan=False, indent=1) + "\n", encoding="utf-8")


def write_vectors(out: Path, config: str, ssd) -> None:
    """Every document vector of one text, one row per line. ``doc_id`` is
    Python's 0-based document position; ``x`` row i belongs to ``doc_id[i]``."""
    x, has_vector, _ = ssd._vectors()
    head = {
        "text": _text_info(config, ssd),
        "doc_id_base": 0,
        "n_docs": len(ssd),
        "n_kept": int(has_vector.sum()),
        "dims": int(x.shape[1]),
        "doc_id": np.flatnonzero(has_vector).tolist(),
    }
    # repr-exact: json writes each float64 as the shortest string that reads back to it.
    rows = ",\n  ".join(json.dumps(r, allow_nan=False) for r in x.tolist())
    meta = json.dumps(head, ensure_ascii=False, indent=1)  # ends in "\n}"
    (out / f"vectors_{config}.json").write_text(meta[:-2] + ',\n "x": [\n  ' + rows + "\n ]\n}\n", encoding="utf-8")


def write_cases(out: Path, emb, data: dict, names=None) -> None:
    """texts.json, the cases of ``case_list`` named in *names* (all when None),
    the two-text, contrast and lexicon cases, and the vectors of the texts they use."""
    df = _columns(data)
    sw = sweep(emb)
    cases = [c for c in case_list(sw) if names is None or c[0] in names]
    used = {c[1] for c in cases} | set(TWO_TEXTS.values())
    tokens = SSD(data["text"], embeddings=emb).docs
    texts = {
        config: SSD(edit_tokens(tokens, c["edit"], emb), embeddings=emb, lexicon=c["lexicon"]) if c["token_lists"]
        else text_column(data, emb, c["lexicon"])
        for config, c in TEXT_CONFIGS.items() if config in used
    }
    out.mkdir(exist_ok=True)
    for old in out.glob("*.json"):  # a case dropped from the list leaves no file behind
        old.unlink()
    _write(out / "texts.json", data)
    for config, text in texts.items():
        write_vectors(out, config, text)
    for name, config, y, call, args, kwargs in cases:
        _write(out / f"{name}.json", fit_case(texts, df, name, config, y, call, args, kwargs))
    for case in (two_texts_case(texts, df, sw), contrast_case(texts, df), *lexicon_cases(texts["text"], df)):
        _write(out / f"{case['case']}.json", case)
    print(f"{out}: {len(list(out.glob('*.json')))} files")


def write_pages(emb, data: dict) -> None:
    """``page_data()`` of one fit per result kind (and the two-text regression and
    the contrast), with ``created`` fixed so reruns only differ where the output changed."""
    df = _columns(data)
    ssd = text_column(data, emb, LEXICON)
    diary = SSD(ssd.docs[::-1], embeddings=emb, lexicon=LEXICON)  # token lists: no snippets
    fits = {
        "pls": ssd.fit_pls(df["y"], k=2, n_splits=20, random_state=42),
        "pca_ols": ssd.fit_regression(df["y"], fixed_k=None, k_min=2, k_max=6, k_step=1),
        "groups": ssd.fit_groups(df["group"], n_perm=50, random_state=42),
        "regression": ssdiff.fit_regression(
            "y ~ ssd(essay) + group*ssd(diary)", df, essay=ssd, diary=diary,
            k_min=2, k_max=6, k_step=1),
        "multipls": ssd.fit_multipls(["y", "y2"], data=df, ssd_gradients=2, n_perm=19, n_splits=5,
                                     random_state=0),
    }
    fits["contrast"] = ssdiff.contrast(fits["pls"].gradient, fits["pca_ols"].gradient, emb,
                                       labels=("pls", "pca_ols"))
    for kind, result in fits.items():
        # Small limits keep the committed files small; the shape is the same at any limit.
        page = page_data(result, title=f"Example ({kind})", snippets_per_side=10,
                         snippets_per_cluster=5, doc_texts=3)
        page["created"] = "2026-01-01T00:00:00Z"
        _write(FIXTURES / "page" / f"{kind}.json", page)
    print(f"{FIXTURES / 'page'}: {len(fits)} files")


def write_uniforms() -> None:
    """Uniforms at random seeds and counter starts, as float.hex() strings."""
    meta = np.random.default_rng(20261027)
    seeds = [int(s) for s in meta.integers(0, 2**53, N_SEEDS)]
    starts = [int(s) for s in meta.integers(0, 2**40, N_SEEDS)]
    # Edge seeds, both at counter 0.
    seeds[:2] = [0, 2**53 - 1]
    starts[:2] = [0, 0]
    cases = [
        {"seed": s, "start": t, "values": [float(u).hex() for u in uniform_at(s, t, N_VALUES)]}
        for s, t in zip(seeds, starts)
    ]
    path = FIXTURES / "rng" / "python_uniforms.json"
    path.write_text("[\n" + ",\n".join(json.dumps(c, separators=(",", ":")) for c in cases) + "\n]\n",
                    encoding="utf-8")
    print(path)


def write_lexicon_small(emb) -> None:
    s = SSD(SMALL_DOCS, embeddings=emb, use_full_doc=False)
    y, g3, g2 = np.asarray(SMALL_Y), np.asarray(SMALL_G3), np.asarray(SMALL_G2)
    results = {
        "continuous": s.suggest_lexicon(y, min_docs=2, top_k=50, n_bins=3),
        "categorical3": s.suggest_lexicon(g3, min_docs=2, var_type="categorical"),
        "seeds": s.suggest_lexicon(seeds=["dog", "cat"], top_k=3, min_docs=1),
        "seeds_y": s.suggest_lexicon(y, seeds=["cat", "dog"], top_k=2, min_docs=1, n_bins=2),
        "evaluate2": s.evaluate_lexicon(g2, ["dog", "cat", "swim", "dog"], var_type="categorical"),
        "evaluate_cont": s.evaluate_lexicon(y, ["dog", "duck"], n_bins=3),
    }
    path = FIXTURES / "lexicon" / "python_small.json"
    _write(path, {"docs": SMALL_DOCS, "y": SMALL_Y, "g3": SMALL_G3, "g2": SMALL_G2,
                  **{k: _lexicon_tables(r) for k, r in results.items()}})
    print(path)


def main(parts: list[str]) -> None:
    if "golden" in parts:
        emb = Embeddings.load(str(TEST_PACK))
        data = make_texts(emb, **GOLDEN_TEXTS)
        write_cases(GOLDEN_DIR, emb, data, GOLDEN)
        write_pages(emb, data)
        write_uniforms()
        write_lexicon_small(emb)
        _write(GOLDEN_DIR / "versions.json", {**{d: version(d) for d in ("ssdiff", "plskit", "glmm", "numpy")},
                                              "pack": PINS["pack"], "testpack": PINS["testpack"]})
    if "reference" in parts:
        path = resolve(PINS["pack"])
        emb = Embeddings.load(path)
        write_cases(OUT, emb, make_texts(emb, **REFERENCE_TEXTS))
        _write(OUT / "pack.json", {"pack": PINS["pack"], "path": path})


if __name__ == "__main__":
    main(sys.argv[1:] or ["golden", "reference"])
