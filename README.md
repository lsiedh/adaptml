AdaptML
=======

Automatically partition a gene phylogeny by using genetic and ecological similarity.

This is a SCELSE fork of Lawrence David's [almlab/adaptml](https://github.com/almlab/adaptml).
It ships **two** ways to run the same analysis:

1. **Python 3 CLI** — the original Python pipeline, ported to modern Python.
2. **`adaptml.html`** — a single self-contained browser app (no install required) that runs the entire AdaptML algorithm in JavaScript and emits files in the **current** iTOL upload format.

Both reproduce Figure 1A–B of [Hunt, David, Gevers, Preheim, Alm & Polz, *Science* 320:1081 (2008)](http://www.ncbi.nlm.nih.gov/pubmed/18497299).

---

## Getting the files from GitHub

**GitHub** is a website that hosts code. You do not need a GitHub account to download this tool — just follow one of the two options below.

### Option A — Download a ZIP (easiest, no Git required)

1. Open **<https://github.com/lsiedh/adaptml>** in your browser.
2. Click the green **"< > Code"** button near the top-right of the page.
3. Choose **"Download ZIP"**.
4. Unzip the downloaded file. You will get a folder called `adaptml-main` (or similar). Rename it to `adaptml` if you like, and note where it is — you will `cd` into it for every command below.

### Option B — Clone with Git (easier to update later)

If you already have Git installed (see prerequisites below), open a terminal and run:

```bash
git clone https://github.com/lsiedh/adaptml.git
cd adaptml
```

To update the code later: `git pull`.

---

## 1. Python 3 CLI

### Prerequisites — install Python 3 and Git

You need **Python 3.9 or newer** and **Git**. Check with `python3 --version` and `git --version`; if either is missing, follow the steps for your OS below.

**macOS**

```bash
# Install Homebrew if you don't have it (https://brew.sh)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Python and Git
brew install python git
```

**Linux (Debian / Ubuntu)**

```bash
sudo apt update
sudo apt install -y python3 python3-venv python3-pip git
```

**Linux (Fedora / RHEL)**

```bash
sudo dnf install -y python3 python3-pip git
```

**Windows**

1. Install Python 3 from [python.org/downloads](https://www.python.org/downloads/) — during setup, tick **"Add python.exe to PATH"**.
2. Install Git from [git-scm.com/download/win](https://git-scm.com/download/win).
3. Run the commands below from **PowerShell** or **Git Bash**, substituting `python` for `python3` and `.venv\Scripts\activate` for `source .venv/bin/activate`.

### Get the code

Download the ZIP or clone with Git as described in the **"Getting the files from GitHub"** section at the top of this document, then `cd` into the folder before running the commands below.

### Install Python dependencies

```bash
python3 -m venv .venv
source .venv/bin/activate           # Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
```

### Run the included Vibrionaceae example

```bash
python3 wrapper/WrAdaptMLFile.py     example/adaptml.file
python3 wrapper/WrapLikelihoodFile.py example/likelihood.file
```

Outputs land under `example/output/`:

| File | Meaning |
|---|---|
| `habitat.matrix` | Inferred emission probability matrix (habitat → environment) |
| `mu.val` | Inferred transition rate µ |
| `stats.file` | EM convergence log |
| `emp_trees/` | Per-iteration randomization likelihoods used for empirical significance |
| `thresh.file` | Per-node empirical likelihood threshold |
| `itol.tree` | Rooted Newick tree for iTOL |
| `full.file` | Full iTOL dataset (legacy format) |
| `cluster.file` | Significant clusters / habitat circles |
| `prune.file`, `prune.tree`, `bars.file` | Pruned-tree summary (Fig 1B equivalent) |
| `habitat.file` | Per-leaf habitat assignment |
| `lik.file` | Joint likelihood + AIC |
| `strain.names` | Member strains of every significant cluster |

### Run on your own data

Make a copy of `example/adaptml.file` and `example/likelihood.file`, edit:

- `tree=…` — Newick tree where leaf names follow `EcologyID_SequenceID` (see `readme.pdf` §2.1.1.1)
- `outgroup=…` — name of the outgroup leaf
- `init_hab_num=16`, `collapse_thresh=0.1`, `converge_thresh=0.001`, `rateopt=avg` — defaults match the published paper
- `color=…` — single-space-delimited color spec; see `example/color.file`
- `thresh=0.9999` — empirical p-value threshold (use 0.9999 for the published p < 0.0001; 0.95 for quicker exploration)

Or call the two component scripts directly (see `readme.pdf` §2 for the full syntax).

---

## 2. `adaptml.html` (single-file browser app)

`adaptml.html` is a self-contained HTML/JS/CSS file. Open it in any modern browser:

```bash
python3 -m http.server 8000
# then visit http://localhost:8000/adaptml.html

# To stop the server:
pkill -f "http.server 8000"
```

(A local web server is needed only so the "Load example" button can `fetch` the bundled `example/vibrio.hsp60.tree`; without it you can still paste a tree by hand.)

The app re-implements the AdaptML algorithm in JavaScript with publication-grade defaults (see below), runs in a Web Worker so the page stays responsive, and produces a ZIP of modern iTOL dataset files ready to drag-and-drop onto your tree at [itol.embl.de](https://itol.embl.de/).

### Publication-grade defaults

| Parameter | Default | Notes |
|---|---|---|
| Initial habitats | 16 | Matches SOM |
| Collapse threshold | 0.10 | Merges similar habitats |
| Convergence threshold | 0.001 | EM stop criterion |
| Empirical significance percentile | **0.9999** | p < 0.0001; use 0.95 for quick exploration (~1 min) |
| Randomization iterations | **10000** | ~10 min on the vibrio dataset |
| Random seed | 2727 | For reproducibility |

### Output files and iTOL upload

Download the ZIP from the results panel and upload to [itol.embl.de](https://itol.embl.de/):

| File | iTOL format | Visualises |
|---|---|---|
| `tree.nwk` | Newick (internal nodes named `N00001…`) | Tree topology — upload this first |
| `dataset_population_strips.txt` | `TREE_COLORS range` | Ecological population clade wedges (alternating blue/gray, numbered) |
| `dataset_habitat_symbols.txt` | `DATASET_SYMBOL` | Habitat circles at cluster root nodes |
| `dataset_ecology_pos*.txt` | `DATASET_COLORSTRIP` | Per-position ecology rings (one file per ecology position) |
| `habitats.json` | JSON | Learned emission probability matrix |
| `populations.tsv` | TSV | Per-population: habitat, leaf count, member names |
| `stats.tsv` | TSV | EM convergence log |

**iTOL tip for population wedges:** after uploading `dataset_population_strips.txt`, click the **"Colored ranges"** dataset in the iTOL panel and set **Cover → Clade** (or **Tree**) so each population fills a solid wedge sector rather than a thin outer arc.

---

## References

- Hunt DE, David LA, Gevers D, Preheim SP, Alm EJ, Polz MF. *Resource Partitioning and Sympatric Differentiation Among Closely Related Bacterioplankton*. Science 320(5879):1081–1085 (2008). doi:10.1126/science.1157890
- Original tool guide: see `readme.pdf` in this repo.
