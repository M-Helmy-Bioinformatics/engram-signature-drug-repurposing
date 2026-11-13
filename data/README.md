# Data

Nothing in this folder is committed except this file and the benchmark template (see `.gitignore`).

## 1. Engram DE table (`de_engram_vs_other.csv`)
Produced by the companion repo `severance-engram-analysis` (`python notebooks/analysis.py` on a real engram scRNA-seq dataset, output `results/de_engram_vs_other.csv`). Needs columns `names`, `logfoldchanges`, `pvals_adj`.

Record the dataset here:

| Field | Value |
|-------|-------|
| GEO accession | _fill in_ |
| Paper | _fill in_ |
| Region / species | _fill in_ |

## 2. Drug signatures (`lincs_signatures.csv`)
Source: LINCS L1000 (GEO series GSE92742 for phase I and GSE70138 for phase II; also available via the CLUE platform, clue.io). The files are large and use the `.gctx` format; the Python package `cmapPy` can read them.

You need a CSV with **rows = compounds, columns = human gene symbols, values = z-scores**, first column = compound name. Start with a subset (for example, one cell line and a limited set of compounds) to keep it manageable. Note that gene identifiers in LINCS files may be Entrez IDs and need converting to symbols.

Record exactly which release, cell lines, doses and times you used.

## 3. Ortholog table (`mouse_human_orthologs.csv`)
Columns `mouse`, `human`. Export from Ensembl BioMart (mouse genes with human orthologs; keep one-to-one orthologs).

## 4. Benchmark compounds (`benchmark_compounds.csv`)
Columns: `compound`, `evidence` (citation). Add only compounds whose memory-enhancing evidence you have **read and verified yourself**. Compound names must match the names used in your drug signatures file.
