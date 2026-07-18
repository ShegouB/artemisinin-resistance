# Artemisinin Resistance Mutation Scanner

A lightweight Python pipeline for screening *Plasmodium falciparum* *kelch13* (*K13*) sequence isolates for WHO-validated artemisinin resistance mutations.

This project is designed for offline use and handles partial Sanger amplicons, reverse-strand inputs, and automatic frame detection.

## What it does

- Detects the correct coding frame and strand using 6-frame translation and 5-mer peptide scoring
- Aligns translated sequences to the PF3D7 K13 reference with BLASTP
- Maps WHO resistance positions through alignment-aware coordinate walking
- Generates both a machine-readable CSV report and an interactive HTML dashboard
- Provides a native Tkinter GUI for non-bioinformaticians
- Supports batch GenBank downloading and West African metadata filtering

## Project layout

```text
artemisinin-resistance/
├── article/
│   ├── article_kelch13_resistance_pipeline.pdf
│   └── kelch13_resistance_pipeline_manuscript.pdf
├── data/
│   ├── reference.fasta              # 3D7 K13 protein reference
│   ├── field_isolates.fasta         # Simulated West African isolates
│   ├── ground_truth.tsv             # Benchmark genotype truth set
│   └── real_field_isolates.fasta    # GenBank clinical cohort
├── results/
│   ├── blast_db/                    # Local BLAST database files
│   ├── blast_hits.tsv               # BLASTP alignment output
│   ├── resistance_report.csv        # Tabular report
│   └── resistance_report.html       # Interactive dashboard
├── scripts/
│   ├── build_field_isolates.py      # Simulate test sequences
│   ├── build_real_isolates.py       # Download and translate GenBank isolates
│   ├── scan_resistance.py           # Core scanner and report generator
│   ├── validate_scanner.py          # Validate outputs against ground truth
│   └── gui_launcher.py              # Tkinter GUI launcher
└── README.md                        # This document
```

## Installation

### Prerequisites

- Python 3.10+
- NCBI BLAST+ (`makeblastdb`, `blastp`)
- Biopython
- `python3-tk` for the GUI on Linux

### Install BLAST+

- Ubuntu/Debian: `sudo apt-get install ncbi-blast+`
- macOS: `brew install blast`
- Windows/Other: download from NCBI

### Create Python environment

```bash
conda create -n bioarchitect python=3.10 -y
conda activate bioarchitect
conda install -c conda-forge biopython -y
```

### Linux GUI dependency

```bash
sudo apt-get install python3-tk -y
```

## Usage

### Option A: Graphical interface

```bash
python3 scripts/gui_launcher.py
```

Steps:

1. Click **Download GenBank Cohort** to fetch the West African isolate set.
2. Choose a local FASTA file.
3. Click **Run Resistance Scan**.
4. Click **Open HTML Dashboard** to view results.

### Option B: Command line

```bash
git clone https://github.com/ShegouB/artemisinin-resistance.git
cd artemisinin-resistance
```

Download and translate the clinical cohort:

```bash
python3 scripts/build_real_isolates.py
```

Run the scanner:

```bash
python3 scripts/scan_resistance.py
```

Outputs:

- `results/resistance_report.csv`
- `results/resistance_report.html`

Validate results:

```bash
python3 scripts/validate_scanner.py
```

## Core algorithms

### 6-frame translation

The pipeline translates each query in six frames (three forward and three reverse-complement) and selects the best match against the 3D7 K13 reference using 5-mer peptide scoring.

### Alignment-aware mapping

The scanner maps WHO resistance positions through BLASTP alignments, handling gaps, insertions, deletions, and truncated amplicons to preserve coordinate accuracy.

## Validation summary

- **Simulated benchmark (n=30):** 100% accuracy on the test set
- **Clinical GenBank cohort (n=66):** 100% concordance with published findings
- Example: isolate `OQ102665` contained the WHO-validated Y493H mutation

## References

1. Ariey, F. et al. (2014). A molecular marker of artemisinin-resistant *Plasmodium falciparum* malaria. *Nature*, 505(7481), 50–55.
2. Ashley, E. A. et al. (2014). Spread of artemisinin resistance in *Plasmodium falciparum* malaria. *NEJM*, 371(5), 411–423.
3. Uwimana, A. et al. (2020). Emergence and clonal expansion of in vitro artemisinin-resistant *Plasmodium falciparum* *kelch13* R561H mutant parasites in Rwanda. *Nature Medicine*, 26(10), 1602–1608.
4. MalariaGEN, Ahouidi, A. et al. (2021). An open dataset of *Plasmodium falciparum* genome variation in 7,000 worldwide samples. *Wellcome Open Research*, 6, 42.
5. Dieng, C. C. et al. (2023). Distribution of *Plasmodium falciparum* *K13* gene polymorphisms across transmission settings in Ghana. *BMC Infectious Diseases*, 23(1), 801.
6. Ajogbasile, F. V. et al. (2022). Molecular profiling of the artemisinin resistance Kelch 13 gene in *Plasmodium falciparum* from Nigeria. *PLOS One*, 17(2), e0264548.
7. Arzika, I. I. et al. (2023). *Plasmodium falciparum* *kelch13* polymorphisms identified after treatment failure with artemisinin-based combination therapy in Niger. *Malaria Journal*, 22, 142.
8. Schmedes, S. et al. (2021). *Plasmodium falciparum* *kelch13* Mutations, 9 Countries in Africa, 2014–2018. *Emerging Infectious Diseases*, 27(7), 1902–1908.
