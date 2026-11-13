"""End-to-end test on SYNTHETIC data (no downloads needed).

This checks that the code works: planted "mimic" drugs should rank at the top
and planted "reverser" drugs at the bottom. It says nothing about real biology.

    python run_demo.py
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent / "src"))
from engram_cmap import build_signature, mouse_to_human, connectivity_scores, benchmark_auc

rng = np.random.default_rng(1)
n_genes, n_drugs = 2000, 800
genes = [f"Gene{i}" for i in range(n_genes)]           # "mouse" style symbols
human = [g.upper() for g in genes]

# Fake DE table: first 60 genes up, next 60 down, rest noise
lfc = rng.normal(0, 0.2, n_genes)
lfc[:60] += 2.0
lfc[60:120] -= 2.0
de = pd.DataFrame({"names": genes, "logfoldchanges": lfc,
                   "pvals_adj": np.where(np.abs(lfc) > 1, 1e-6, 0.5)})

up_m, down_m = build_signature(de, n=50)
up, down = mouse_to_human(up_m), mouse_to_human(down_m)

# Fake drug matrix: noise, plus 10 planted mimics and 10 planted reversers
Z = rng.normal(0, 1, (n_drugs, n_genes))
direction = np.zeros(n_genes)
direction[:60], direction[60:120] = 1, -1
mimics = [f"mimic_{i}" for i in range(10)]
reversers = [f"reverser_{i}" for i in range(10)]
names = mimics + reversers + [f"drug_{i}" for i in range(n_drugs - 20)]
Z[:10] += direction * 1.0
Z[10:20] -= direction * 1.0
drugs = pd.DataFrame(Z, index=names, columns=human)

scores = connectivity_scores(drugs, up, down, n_perm=300)
res = benchmark_auc(scores, mimics)

print(f"Signature genes used: up={scores.attrs['n_up_used']}, down={scores.attrs['n_down_used']}")
print("Top 10 compounds:\n", scores.head(10)[["raw_score", "perm_z"]].round(2))
print("Benchmark (planted mimics):", res)
print("Bottom 3:", scores.tail(3).index.tolist())

Path("results").mkdir(exist_ok=True)
scores.to_csv("results/demo_ranking.csv")
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(scores["perm_z"].dropna(), bins=50, color="#999999")
for name in mimics:
    ax.axvline(scores.loc[name, "perm_z"], color="#2a9d8f", alpha=0.6)
ax.set_xlabel("Connectivity z-score")
ax.set_ylabel("Compounds")
ax.set_title("SYNTHETIC demo: planted mimics (green) rank at the top")
plt.tight_layout()
fig.savefig("results/demo_ranking.png", dpi=200)
