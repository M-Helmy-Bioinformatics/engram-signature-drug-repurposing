"""Run the pipeline on REAL inputs.

    python run_pipeline.py \
        --de data/de_engram_vs_other.csv \
        --drugs data/lincs_signatures.csv \
        --benchmark data/benchmark_compounds.csv \
        --ortholog-map data/mouse_human_orthologs.csv

Inputs (see data/README.md):
  --de            DE table (columns: names, logfoldchanges, pvals_adj)
  --drugs         CSV, rows = compounds, columns = human gene symbols,
                  values = z-scores (first column = compound name)
  --benchmark     CSV with a `compound` column (known memory-enhancing
                  compounds, each with a literature citation you verified)
  --ortholog-map  optional CSV with columns `mouse`, `human`
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent / "src"))
from engram_cmap import build_signature, mouse_to_human, connectivity_scores, benchmark_auc

ap = argparse.ArgumentParser()
ap.add_argument("--de", required=True)
ap.add_argument("--drugs", required=True)
ap.add_argument("--benchmark")
ap.add_argument("--ortholog-map")
ap.add_argument("--n", type=int, default=100, help="genes per direction")
ap.add_argument("--out", default="results/ranking.csv")
args = ap.parse_args()

de = pd.read_csv(args.de)
up_m, down_m = build_signature(de, n=args.n)
mapping = None
if args.ortholog_map:
    om = pd.read_csv(args.ortholog_map)
    mapping = dict(zip(om["mouse"], om["human"]))
up, down = mouse_to_human(up_m, mapping), mouse_to_human(down_m, mapping)
print(f"Signature: {len(up)} up, {len(down)} down (human symbols)")

drugs = pd.read_csv(args.drugs, index_col=0)
scores = connectivity_scores(drugs, up, down)
Path(args.out).parent.mkdir(exist_ok=True)
scores.to_csv(args.out)
print(scores.head(20)[["raw_score", "perm_z"]].round(2))

if args.benchmark:
    pos = pd.read_csv(args.benchmark)["compound"].dropna().tolist()
    print("Benchmark:", benchmark_auc(scores, pos))
print(f"Saved {args.out}")
