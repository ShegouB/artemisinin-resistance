# Artemisinin Resistance Mutation Scanner

An open-source, lightweight, and offline-deployable Python pipeline to screen *Plasmodium falciparum* *kelch13* (*K13*) sequence isolates for WHO-validated artemisinin resistance mutations.

By utilizing a resilient **6-frame translation algorithm powered by peptide match scoring** and **alignment-aware coordinate mapping**, this tool processes partial Sanger sequence amplicons and supports reverse-strand inputs without manual curation.

## Key Features

- **Offline-First Design:** Operates entirely locally. Ideal for regional malaria reference laboratories with limited internet connectivity.
- **Lightweight GUI Launcher:** Features a native, dependency-minimal Tkinter graphical interface designed for non-bioinformaticians, with real-time log output and asynchronous multi-threaded execution.
- **6-Frame Translation Heuristic:** Detects the correct coding frame and strand (direct or reverse-complement) of Sanger amplicons using peptide matching against the 3D7 reference.
- **Alignment-Aware Coordinate Mapping:** Maps WHO resistance positions through BLASTP alignment HSPs, handling gaps, insertions, and truncations.
- **Batch GenBank Fetching:** Programmatic batch download and filtering to build enriched clinical surveillance cohorts directly from NCBI GenBank.
- **Geographical Metadata Filtering:** Extracts source-country metadata from GenBank records and retains West African isolates.
- **Multi-Format Reporting:** Produces both machine-readable CSV reports and a polished interactive HTML dashboard.

## Directory Structure

```text
artemisinin-resistance/
├── data/
│   ├── reference.fasta              # Curated PF3D7_1343700 protein reference
│   ├── field_isolates.fasta         # Simulated panel of 30 West African isolates
│   ├── ground_truth.tsv             # Simulated genotype ground truth mapping
│   └── real_field_isolates.fasta    # GenBank clinical cohort (n=66)
├── results/
│   ├── blast_db/                    # Local protein BLAST database files
│   ├── blast_hits.tsv               # Tabular BLASTP alignment output
│   ├── resistance_report.csv        # Tabular diagnostic report
│   └── resistance_report.html       # Interactive HTML dashboard
├── scripts/
│   ├── build_field_isolates.py      # Script to simulate 30 mock sequences
│   ├── build_real_isolates.py       # Programmatic GenBank cohort builder
│   ├── scan_resistance.py           # Core BLAST/scanner pipeline
│   ├── validate_scanner.py          # Validation script against ground truth
│   └── gui_launcher.py              # Native Tkinter GUI launcher
└── README.md                        # Documentation
```

## Installation & Prerequisites

This pipeline is built on Python 3 and standard scientific computing dependencies. It requires a local installation of the **NCBI BLAST+** suite.

### 1. Install NCBI BLAST+

Ensure `makeblastdb` and `blastp` are installed and available in your system path:

- **Ubuntu/Debian:** `sudo apt-get install ncbi-blast+`
- **macOS (via Homebrew):** `brew install blast`
- **Windows/Other:** Download installers from [NCBI BLAST+ FTP](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

### 2. Configure Python Environment

It is recommended to manage dependencies within a dedicated Conda environment.

```bash
conda create -n bioarchitect python=3.10 -y
conda activate bioarchitect
conda install -c conda-forge biopython -y
```

If you are not using Conda, you can also install Biopython with pip:

```bash
pip install biopython
```

If you plan to use the GUI on Ubuntu/Debian Linux, install Tkinter:

```bash
sudo apt-get install python3-tk -y
```

## Workflow & Quick Start

You can run the complete analytical pipeline either through the **Command Line Interface (CLI)** or the **Graphical User Interface (GUI)**.

### Option A: Using the Graphical User Interface (Recommended for Clinicians)

To launch the interactive GUI:

```bash
python3 scripts/gui_launcher.py
```

From the interface, you can:
1. Download the GenBank cohort.
2. Select a local FASTA file.
3. Run the resistance scan.
4. Open the HTML dashboard.

### Option B: Using the Command Line Interface (CLI)

#### Step 1: Clone the Repository

```bash
git clone https://github.com/ShegouB/artemisinin-resistance.git
cd artemisinin-resistance
```

#### Step 2: Build the Real Clinical Cohort from GenBank

```bash
python3 scripts/build_real_isolates.py
```

*Note: The script uses an e-mail address (`djagouboris@gmail.com`) for NCBI E-utilities; edit it if needed.*

#### Step 3: Run the Mutation Scanner

```bash
python3 scripts/scan_resistance.py
```

- **Output CSV:** `results/resistance_report.csv`
- **Output HTML:** `results/resistance_report.html`

#### Step 4: Validate the Diagnostic Results

```bash
python3 scripts/validate_scanner.py
```
