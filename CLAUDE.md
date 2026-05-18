# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

SCELSE fork of the Alm Lab AdaptML tool (Hunt et al., *Science* 320:1081, 2008). Two parallel implementations of the same algorithm:

1. **Python 3 CLI** — original pipeline, ported from Python 2.5 to modern Python 3 + NumPy/SciPy
2. **`adaptml.html`** — complete JS reimplementation in a single self-contained browser file

## Running the browser app

```bash
python3 -m http.server 8000
# visit http://localhost:8000/adaptml.html
```

A local server is required only for the "Load example" button (`fetch` is blocked on `file://`). Without it, users can paste the tree manually.

## Running the Python CLI

```bash
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt

# End-to-end on the Vibrio example:
python3 wrapper/WrAdaptMLFile.py example/adaptml.file      # Step 1: habitat learning
python3 wrapper/WrapLikelihoodFile.py example/likelihood.file  # Step 2: cluster identification
```

Outputs go to `example/output/`. The two config files (`adaptml.file`, `likelihood.file`) are the canonical entry points; the underlying scripts can also be invoked directly with keyword= arguments.

## Architecture

### Python pipeline — two-stage design

**Stage 1 — Habitat learning** (`habitats/trunk/AdaptML.py`)
- Reads an unrooted Newick tree + optional params
- EM loop: `LearnLiks` (Felsenstein pruning) → `EstimateStates` (marginal posterior per node) → `LearnRates` (update habitat emission matrix + µ)
- Merges similar habitats (sum-of-squared-diffs < `collapse_thresh`)
- Writes `habitat.matrix` and `mu.val`

**Stage 2 — Cluster identification** (`clusters/trunk/JointML.py`)
- Reads outputs from Stage 1
- Joint Viterbi assigns one best habitat to every node
- Leaf-label shuffling (`clusters/getstats/rand_JointML.py`) builds per-node empirical null distributions
- `GetLikelihoods.py` converts randomization outputs to per-node likelihood thresholds
- Significant clusters: subtree ≥ 90% habitat-coherent AND joint likelihood > empirical threshold (`thresh`)
- Writes iTOL files: `full.file`, `cluster.file`, `prune.file`, `bars.file`, `itol.tree`, `strain.names`

**Shared modules** (both stages have their own copy under `habitats/trunk/` and `clusters/trunk/`):
- `ML.py` — Felsenstein pruning, `LearnLiks`, `EstimateStates`, `LearnRates`/`NodeRate`
- `multitree.py` — tree container, `LeafShuffle` for randomization
- `node.py` — node class, all recursive tree operations (`PieCharts`, `FulliTol`, `ClusterTest`, etc.)
- `branch.py` — branch class

`clusters/getstats/` mirrors the cluster stage but with `rand_` prefix (randomization variants that shuffle leaf labels and recompute likelihoods).

### Wrapper scripts

`wrapper/WrAdaptMLFile.py` / `WrapLikelihoodFile.py` — parse the `.file` configs and call the underlying scripts via `subprocess` using `sys.executable`.

### Browser app (`adaptml.html`)

Single file (~1900 lines), all CSS + JS inlined. Structure:
- **Newick parser** — handles bootstrap labels, zero-length branches
- **`rootify()`** — midpoint-splits the outgroup edge to create a rooted tree; assigns `parent`/`children`/`leaf_nodes` on every node
- **Habitat EM** — JS port of `habitats/trunk/ML.py` + `AdaptML.py`; runs in an inline Web Worker (Blob URL) so the UI stays responsive; posts progress events
- **Joint Viterbi** — JS port of `clusters/trunk/ML.py`
- **`runRandomizations()`** — shuffles leaf ecology labels N times, records per-node likelihoods, computes percentile threshold
- **`findClusters()`** — post-order divergence detection + `tryCluster` descent (≥90% coherence + likelihood > threshold); each node guarded with `if (!node.cluster_root)` to prevent double-counting
- **Output writers** — `writeNewick` (embeds `N00001…` internal IDs), `writeTreeColorsDataset` (TREE_COLORS range — population clade bands), `writeSymbolsDataset` (DATASET_SYMBOL — habitat circles), `writeColorstripDataset` (DATASET_COLORSTRIP — per-position ecology rings), `writeMultibarDataset` (DATASET_MULTIBAR), plus JSON/TSV outputs and an inline ZIP writer (STORED/uncompressed, CRC32)

### Key algorithm parameters (publication defaults)

| Parameter | Default | Notes |
|---|---|---|
| `init_hab_num` | 16 | Initial habitat count for EM |
| `collapse_thresh` | 0.10 | Merge habitats with sum-sq-diff < this |
| `converge_thresh` | 0.001 | Stop EM when habitat matrix change < this |
| `thresh` / `thresh_p` | 0.9999 | Empirical p-value (p < 0.0001); use 0.95 for quick exploration |
| `rand_iters` | 10000 | Randomization iterations for empirical null |
| `seed` | 2727 | RNG seed for reproducibility |
| Coherence cutoff | 90% | Subtree must be ≥90% same-habitat to qualify as a population |

### numpy 2.x gotcha

`from numpy import *` shadows builtin `sum`. Any file that sums `dict.values()` must use `_builtin_sum`:
```python
from builtins import sum as _builtin_sum
scale = _builtin_sum(habitat_matrix[habitat].values())
```
This is already applied in `habitats/trunk/AdaptML.py` and `clusters/trunk/JointML.py`.

### iTOL output files

| File | iTOL format | Visualises |
|---|---|---|
| `tree.nwk` | Newick (internal nodes named `N00001…`) | Tree topology |
| `dataset_population_strips.txt` | `TREE_COLORS range` | Population clade wedges (alternating blue/gray) |
| `dataset_habitat_symbols.txt` | `DATASET_SYMBOL` | Habitat circles at cluster roots |
| `dataset_ecology_pos*.txt` | `DATASET_COLORSTRIP` | Per-position ecology rings (one file per position) |
| `habitats.json` | JSON | Emission probability matrix |
| `populations.tsv` | TSV | Cluster members |

In iTOL, set **"Colored ranges" → Cover → Clade** (or Tree) to render `dataset_population_strips.txt` as filled wedge sectors rather than thin arcs.
