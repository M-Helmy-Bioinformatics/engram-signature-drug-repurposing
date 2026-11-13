# Engram Signature Drug Repurposing

**Status: work in progress.** The code runs end to end on synthetic data. Results on real data are not yet available (see [Roadmap](#roadmap)).

## Question

Neurons that encode a memory ("engram cells") switch on a characteristic set of genes. Can public drug-response data suggest existing compounds that push cells toward that same gene-expression pattern, as *hypotheses* for memory modulators?

## Approach

1. **Engram signature.** From a public hippocampal single-cell RNA-seq dataset, label engram-like cells by immediate-early-gene (IEG) score (*Fos*, *Npas4*, *Arc*, *Egr1*), run differential expression (engram-like vs other), and keep the top up- and down-regulated genes. The IEGs used for labeling are excluded from the signature to avoid circularity.
2. **Ortholog mapping.** Convert mouse gene symbols to human (an ortholog table such as an Ensembl BioMart export is recommended; the built-in upper-casing is only a rough fallback).
3. **Connectivity scoring.** For each compound, `score = mean z(up genes) - mean z(down genes)`, converted to a permutation z-score using random gene sets of equal size. High positive = same direction as the signature ("mimic"); strongly negative = opposite ("reverser").
4. **Benchmark.** Test whether compounds with published memory-enhancing evidence rank above other compounds (AUC, Mann-Whitney test). The benchmark list must be curated from the literature, with every entry's citation verified (`data/benchmark_compounds.csv`).

## Repository layout

```
├── src/engram_cmap.py        # signature, ortholog mapping, scoring, benchmark
├── run_demo.py               # synthetic end-to-end test (no downloads)
├── run_pipeline.py           # run on real inputs
├── data/                     # inputs (not committed; see data/README.md)
├── results/                  # outputs
├── requirements.txt
└── LICENSE
```

## Quick start

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate      Mac/Linux: source .venv/bin/activate
pip install -r requirements.txt
python run_demo.py
```

The demo plants 10 "mimic" compounds in random data and checks they rank at the top. **It validates the code, not the biology.**

## Running on real data

See `data/README.md` for where to get the inputs. Then:

```bash
python run_pipeline.py --de data/de_engram_vs_other.csv \
    --drugs data/lincs_signatures.csv \
    --benchmark data/benchmark_compounds.csv \
    --ortholog-map data/mouse_human_orthologs.csv
```

The DE table comes from the companion repository (`severance-engram-analysis`, `notebooks/analysis.py`).

## Limitations

- IEG-high cells mark recent neuronal activity, not necessarily a specific memory.
- LINCS L1000 profiles come largely from non-neuronal (mostly cancer) cell lines, so effects in neurons may differ.
- Connectivity scoring yields many false positives. Output is a ranked list of **candidate hypotheses**, not findings.
- Nothing here is validated experimentally.
- Mouse-to-human ortholog mapping loses genes and can mis-map paralogs.

## Roadmap

- [ ] Choose and process a real engram scRNA-seq dataset (record GEO accession)
- [ ] Obtain LINCS L1000 signatures and a subset suitable for the analysis
- [ ] Build and verify the benchmark compound list with citations
- [ ] Run, inspect, and report results (figures in `results/`)
- [ ] Sensitivity checks (signature size, scoring method)

## AI-assistance note

The initial code scaffold was developed with assistance from an AI model (Claude). The author is responsible for reviewing, testing and understanding all code and results.

## Data and credits

- LINCS L1000: Subramanian et al., *Cell* 2017 (cite the exact release used).
- Engram scRNA-seq: cite the source study and GEO accession you use.

## License

MIT (code). Data remain under their original terms.
