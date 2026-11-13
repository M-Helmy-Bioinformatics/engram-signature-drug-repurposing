"""Core functions: engram signature -> drug connectivity ranking.

Pipeline
--------
1. build_signature():   differential-expression table -> up / down gene lists
2. mouse_to_human():    map mouse gene symbols to human symbols
3. connectivity_scores(): score every drug profile against the signature
4. benchmark_auc():     do known memory-enhancing compounds rank above others?

Input conventions
-----------------
- DE table: CSV from scanpy's rank_genes_groups_df with columns
  `names`, `logfoldchanges`, `pvals_adj` (the Severance repo writes this as
  results/de_engram_vs_other.csv).
- Drug matrix: pandas DataFrame, rows = compounds, columns = HUMAN gene
  symbols, values = differential-expression z-scores after treatment
  (for example LINCS L1000 level-5 signatures).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# Genes used to LABEL engram-like cells. They are excluded from the signature
# by default, otherwise the signature is partly circular (cells were selected
# for high IEG expression, so IEGs are trivially "differential").
LABEL_GENES = {"Fos", "Npas4", "Arc", "Egr1"}


def build_signature(de: pd.DataFrame, n: int = 100, padj: float = 0.05,
                    exclude=LABEL_GENES) -> tuple[list[str], list[str]]:
    """Return (up_genes, down_genes), each up to n genes, ranked by log fold change."""
    need = {"names", "logfoldchanges", "pvals_adj"}
    if not need.issubset(de.columns):
        raise ValueError(f"DE table needs columns {sorted(need)}")
    d = de.dropna(subset=["names", "logfoldchanges", "pvals_adj"])
    d = d[(d["pvals_adj"] < padj) & (~d["names"].isin(exclude))]
    up = d[d["logfoldchanges"] > 0].sort_values("logfoldchanges", ascending=False)
    down = d[d["logfoldchanges"] < 0].sort_values("logfoldchanges")
    return up["names"].head(n).tolist(), down["names"].head(n).tolist()


def mouse_to_human(symbols, mapping: dict[str, str] | None = None) -> list[str]:
    """Map mouse -> human gene symbols.

    If `mapping` (mouse symbol -> human symbol, e.g. exported from Ensembl
    BioMart) is given, it is used and unmapped genes are dropped. Otherwise a
    rough heuristic is used (upper-case the symbol), which works for many but
    NOT all genes. Use a real ortholog table for any serious result.
    """
    if mapping is not None:
        return [mapping[s] for s in symbols if s in mapping]
    return [s.upper() for s in symbols]


def connectivity_scores(drugs: pd.DataFrame, up: list[str], down: list[str],
                        n_perm: int = 500, seed: int = 0) -> pd.DataFrame:
    """Score each compound against the signature.

    raw score  = mean z(up genes) - mean z(down genes)
    A high positive score means the drug moves genes in the same direction as
    the engram signature (candidate "mimic"); a strongly negative score means
    the opposite direction (candidate "reverser").

    The raw score is converted to a permutation z-score using random gene sets
    of the same sizes, so scores are comparable across compounds.
    """
    up = [g for g in up if g in drugs.columns]
    down = [g for g in down if g in drugs.columns]
    if len(up) < 5 or len(down) < 5:
        raise ValueError(f"Too few signature genes found in drug matrix "
                         f"(up={len(up)}, down={len(down)}). Check gene symbols.")
    Z = drugs.to_numpy(dtype=float)
    cols = {g: i for i, g in enumerate(drugs.columns)}
    iu = np.array([cols[g] for g in up])
    idn = np.array([cols[g] for g in down])
    raw = Z[:, iu].mean(axis=1) - Z[:, idn].mean(axis=1)

    rng = np.random.default_rng(seed)
    null = np.empty((n_perm, Z.shape[0]))
    n_genes = Z.shape[1]
    for k in range(n_perm):
        a = rng.choice(n_genes, size=len(iu), replace=False)
        b = rng.choice(n_genes, size=len(idn), replace=False)
        null[k] = Z[:, a].mean(axis=1) - Z[:, b].mean(axis=1)
    mu, sd = null.mean(axis=0), null.std(axis=0)
    sd[sd == 0] = np.nan
    out = pd.DataFrame({"raw_score": raw, "perm_z": (raw - mu) / sd},
                       index=drugs.index)
    out["rank"] = out["perm_z"].rank(ascending=False, method="min")
    out.attrs["n_up_used"], out.attrs["n_down_used"] = len(up), len(down)
    return out.sort_values("perm_z", ascending=False)


def benchmark_auc(scores: pd.DataFrame, positives: list[str]) -> dict:
    """AUC: probability that a known positive outranks a random other compound.

    0.5 = no better than chance; 1.0 = all positives above all others.
    """
    pos = [p for p in positives if p in scores.index]
    if len(pos) < 3:
        return {"n_positives_found": len(pos), "auc": float("nan"),
                "p_value": float("nan")}
    neg = scores.index.difference(pos)
    a = scores.loc[pos, "perm_z"].dropna()
    b = scores.loc[neg, "perm_z"].dropna()
    u, p = mannwhitneyu(a, b, alternative="greater")
    return {"n_positives_found": len(pos), "auc": float(u / (len(a) * len(b))),
            "p_value": float(p)}
