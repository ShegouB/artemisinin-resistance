# Artemisinin Resistance Mutation Scanner

An open-source, lightweight, and offline-deployable Python pipeline to screen *Plasmodium falciparum* *kelch13* (*K13*) sequence isolates for WHO-validated artemisinin resistance mutations. 

By utilizing a resilient **6-frame translation algorithm powered by a 5-mer peptide match-scoring heuristic** and **alignment-aware coordinate mapping**, this tool robustly processes partial Sanger sequence amplicons and handles reverse-strand sequencing inputs without manual curation.

## Key Features

- **Offline-First Design:** Operates entirely locally. Ideal for regional malaria reference laboratories with limited internet connectivity.
- **Lightweight GUI Launcher:** Features a native, dependency-minimal Tkinter graphical interface designed for non-bioinformaticians, with real-time log capturing and asynchronous multi-threaded executions.
- **6-Frame Translation Heuristic:** Automatically detects the correct coding frame and strand (direct or reverse-complement) of Sanger amplicons using 5-mer match-scoring against the 3D7 reference.
- **Alignment-Aware Coordinate Mapping:** Decouples mutation screening from static coordinate indexing by walking synchronously along BLASTP alignment HSP columns (handling gaps, insertions, and truncations).
- **Batch GenBank Fetching:** Programmatic batch download and filtering module to build geographically enriched clinical surveillance cohorts directly from NCBI GenBank.
- **Geographical Metadata Filtering:** Automatically extracts source countries from GenBank record features or descriptions, strictly retaining West African isolates.
- **Multi-Format Reporting:** Outputs diagnostic results in both machine-readable CSV format and an interactive, highly polished HTML dashboard.

## Directory Structure

```text
artemisinin-resistance/
├── data/
│   ├── reference.fasta              # Curated PF3D7_1343700 3D7 protein reference
│   ├── field_isolates.fasta         # Simulated panel of 30 West African isolates
│   ├── ground_truth.tsv             # Simulated genotype ground truth mapping
│   └── real_field_isolates.fasta    # GenBank clinical cohort (n=66)
├── results/
│   ├── blast_db/                    # Local protein BLAST database files
│   ├── blast_hits.tsv               # Tabular BLASTP alignment output
│   ├── resistance_report.csv        # Tabular diagnostic report
│   └── resistance_report.html       # Polished interactive HTML dashboard
├── scripts/
│   ├── build_field_isolates.py      # Script to simulate 30 mock sequences
│   ├── build_real_isolates.py       # Programmatic download/6-frame translation of GenBank isolates
│   ├── scan_resistance.py           # Core database constructor, BLAST aligner, and mutation scanner
│   ├── validate_scanner.py          # Validation script comparing scanner outputs to ground truth
│   └── gui_launcher.py              # Native asynchronous Tkinter GUI for non-bioinformaticians
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
It is highly recommended to manage dependencies within a dedicated Conda environment.

```bash
# Create and activate environment
conda create -n bioarchitect python=3.10 -y
conda activate bioarchitect

# Install Biopython
conda install -c conda-forge biopython -y
```

*Note for Linux GUI users:* If you plan to run the graphical interface on Linux (Ubuntu/Debian), ensure `python3-tk` is installed on your system:
```bash
sudo apt-get install python3-tk -y
```

## Workflow & Quick Start

You can run the complete analytical pipeline either through the **Command Line Interface (CLI)** or the **Graphical User Interface (GUI)**.

### Option A: Using the Graphical User Interface (Recommended for Clinicians)
To launch the interactive GUI, run the following command in your terminal:
```bash
python3 scripts/gui_launcher.py
```
From the interface, you can:
1. Click **"Download GenBank Cohort"** to fetch and translate the 66 West African isolates.
2. Select any local FASTA file (including your own sequencing data) using the **"Browse"** utility.
3. Click **"Run Resistance Scan"** to perform local alignments and screen for mutations.
4. Click **"Open HTML Dashboard"** to instantly render the polished results in your browser.


### Option B: Using the Command Line Interface (CLI)

#### Step 1: Clone the Repository
```bash
git clone https://github.com/ShegouB/artemisinin-resistance.git
cd artemisinin-resistance
```

#### Step 2: Build the Real Clinical Cohort from GenBank
Run the batch downloader to programmatically fetch, translate (resolving strand offsets), and filter 66 West African clinical sequences representing Benin, Niger, Nigeria, Ghana, and Senegal:
```bash
python3 scripts/build_real_isolates.py
```
*Note: The script utilizes an e-mail address (`djagouboris@gmail.com`) as required by NCBI E-utilities. You may edit this directly in the script.*

#### Step 3: Run the Mutation Scanner
Execute the core pipeline to construct the local database, align query amplicons, map coordinates dynamically, and generate standard reports:
```bash
python3 scripts/scan_resistance.py
```
- **Standard CSV Report:** `results/resistance_report.csv`
- **Polished HTML Dashboard:** `results/resistance_report.html`

#### Step 4: Validate the Diagnostic Results
Verify scanner outputs against the validation catalogue (including simulated ground-truths and real GenBank genotypes):
```bash
python3 scripts/validate_scanner.py
```

## Core Algorithms

### 1. Multi-Frame Heuristic Translation (6-Frame)
To determine the correct open reading frame without prior sequence annotations, the script translates both the direct nucleotide strand $N_{fwd}$ and its reverse-complement $N_{rev}$ across three forward frames, generating six candidate protein sequences $T_{s,f}$ (where $s \in \{fwd, rev\}$ and $f \in \{0, 1, 2\}$) [5, 6]:

$$T_{s,f} = \text{translate}\left(N_{s}[f : 3 \times \lfloor (|N_{s}| - f) / 3 \rfloor]\right)$$

To rank these candidates, we compute a matching score based on the occurrence of 5-mer peptides within the 3D7 reference sequence $R_{3D7}$:

$$S(T_{s,f}) = \sum_{i=1}^{|T_{s,f}|-4} I\left(T_{s,f}[i : i+5] \subset R_{3D7}\right)$$

The sequence with the highest score $S \ge 20$ is selected as the correct ORF. Sequences failing to reach this threshold (e.g. non-homologous contaminants) are automatically rejected.

### 2. Alignment-Aware Coordinate Walking
Traditional index tracking fails on truncated amplicons or insertions/deletions. Instead of static coordinate lookups, our dynamic algorithm uses BLASTP aligned query ($qseq$) and reference ($sseq$) HSP strings to walk along the alignment. It tracks subject coordinates through gaps (`-`), mapping the target WHO position to the exact residue index of the query isolate.

## Diagnostic Validation Results

The pipeline has been benchmarked using a dual-validation strategy:

### A. Simulated Benchmark (n=30)
- **Design:** Mock sequences carrying randomized validated mutations with a 23.3% resistant prevalence.
- **Performance:** Correctly called all 30 sequences (7 mutants, 23 wild-type), achieving **100.0% logical sensitivity and specificity** (30/30 OK).

### B. Geographically Enriched Clinical Panel (n=66)
- **Design:** Real Sanger sequence amplicons downloaded from GenBank, representing Benin, Niger, Nigeria, Ghana, and Senegal [5, 6, 7, 8].
- **Performance:** Resolved previous reverse-strand alignment anomalies (`N/A` statuses). Successfully identified 65 sequences as wild-type and automatically detected the presence of the WHO-validated **Y493H** resistance mutation in clinical isolate **OQ102665** from northern Ghana, showing **100.0% validation concordance** (66/66 OK) with published clinical findings [5].

## References

1. **Ariey, F. et al.** (2014). A molecular marker of artemisinin-resistant *Plasmodium falciparum* malaria. *Nature*, 505(7481), 50–55. https://doi.org/10.1038/nature12876
2. **Ashley, E. A. et al.** (2014). Spread of artemisinin resistance in *Plasmodium falciparum* malaria. *NEJM*, 371(5), 411–423. https://doi.org/10.1056/NEJMoa1314981
3. **Uwimana, A. et al.** (2020). Emergence and clonal expansion of *kelch13* R561H mutant parasites in Rwanda. *Nature Medicine*, 26(7), 1061–1064. https://doi.org/10.1038/s41591-020-0896-3
4. **MalariaGEN, Ahouidi, A. et al.** (2021). An open dataset of *Plasmodium falciparum* genome variation in 7,000 worldwide samples. *Wellcome Open Research*, 6, 42. https://doi.org/10.12688/wellcomeopenres.16168.1
5. **Dieng, C. C., et al.** (2023). Distribution of *Plasmodium falciparum* *K13* gene polymorphisms across transmission settings in Ghana. *Malaria Journal*, 22, 345. https://doi.org/10.1186/s12936-023-04768-4
6. **Ajogbasile, F. V., et al.** (2022). Molecular profiling of the artemisinin resistance Kelch 13 gene in *Plasmodium falciparum* from Nigeria. *PLOS One*, 17(2), e0264548. https://doi.org/10.1371/journal.pone.0264548
7. **Arzika, I. I. et al.** (2023). *Plasmodium falciparum* *kelch13* polymorphisms identified after treatment failure with artemisinin-based combination therapy in Niger. *Malaria Journal*, 22, 142. https://doi.org/10.1186/s12936-023-04565-x
8. **Schmedes, S. et al.** (2021). *Plasmodium falciparum* *kelch13* Mutations, 9 Countries in Africa, 2014–2018. *Emerging Infectious Diseases*, 27(1), 295–298. https://doi.org/10.3201/eid2701.201502
