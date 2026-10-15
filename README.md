# 🧬 Triplet Repeat Detection and DNA Sequence Analysis Using Python

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![pandas](https://img.shields.io/badge/pandas-data%20analysis-150458?logo=pandas&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-visualisation-11557c)
![Tests](https://img.shields.io/badge/pytest-5%20passing-brightgreen?logo=pytest&logoColor=white)
![Course](https://img.shields.io/badge/Course-Introduction%20to%20Biology-8A2BE2)

A Python project that **detects trinucleotide (triplet) repeats in DNA** and analyses **CAG, CGG and GAA repeat loci** across the human genome using a published open dataset. These three repeat motifs are linked to **Huntington's disease (CAG)**, **Fragile X syndrome (CGG)** and **Friedreich's ataxia (GAA)**.

> Mini project for **Introduction to Biology** · Team 14 · Amrita Vishwa Vidyapeetham, Coimbatore

---

## 📌 Table of contents

- [What this project does](#-what-this-project-does)
- [How triplet detection works](#-how-triplet-detection-works)
- [How GC content is calculated](#-how-gc-content-is-calculated)
- [Key results](#-key-results)
- [Result graphs](#-result-graphs)
- [Project structure](#-project-structure)
- [Getting started](#-getting-started)
- [Scripts at a glance](#-scripts-at-a-glance)
- [Data validation findings](#-data-validation-findings)
- [Dataset and credits](#-dataset-and-credits)
- [Limitations and future work](#-limitations-and-future-work)
- [Author](#-author)

---

## 🔍 What this project does

The project has two parts:

**1. Raw DNA sequence analyser**
Give it any DNA sequence (A, T, G, C) and it will:
- validate and clean the input (case and whitespace are normalised, invalid characters are rejected)
- count each nucleotide and calculate **GC content**
- split the sequence into triplets and count **triplet frequencies**
- find **consecutive identical triplets** (for example `CAG CAG CAG CAG`) with their start and end positions
- save a text summary and a CSV of the detected repeats

**2. Genome-wide repeat dataset analysis**
Using CAG, CGG and GAA repeat data from Dryad, it:
- validates the data (schema, missing values, duplicates, anomalies)
- summarises repeat lengths per motif and plots their distributions
- compares the three motifs side by side
- measures how many loci are **stable vs polymorphic**
- shows **where in the genome** each motif sits (exonic, intronic, 5′UTR and so on)

---

## 🧠 How triplet detection works

1. **Validate**: strip whitespace, convert to upper case, accept only `A`, `T`, `G`, `C`.
2. **Compose**: count bases and compute GC content = (G + C) / length × 100.
3. **Split** the sequence into non-overlapping triplets in reading frame 0. Leftover bases at the end are reported, not silently dropped.
4. **Scan** the triplet list for runs of two or more identical consecutive triplets and record the motif, repeat count and base positions.

**Example**: `ATGCAGCAGCAGCAGTTACGGCGG`

```text
Sequence length: 24
A: 6   T: 3   G: 9   C: 6        GC content: 62.50%
Triplets: ATG CAG CAG CAG CAG TTA CGG CGG

CAG repeats 4 times (bases 4-15)
CGG repeats 2 times (bases 19-24)
```

---

## 🧮 How GC content is calculated

**GC content** is the percentage of bases in a DNA sequence that are **G (guanine)** or **C (cytosine)**.

**Formula**

```text
GC content (%) = (G + C) / (A + T + G + C) × 100
```

`A`, `T`, `G` and `C` are the number of times each base appears in the sequence. Because the sequence is validated to contain only these four letters, `A + T + G + C` is the same as the sequence length.

**Worked example:** `ATGCAGCAGCAGCAGTTACGGCGG`

| Base | A | T | G | C | Length |
| :---: | :---: | :---: | :---: | :---: | :---: |
| Count | 6 | 3 | 9 | 6 | 24 |

```text
GC content = (9 + 6) / 24 × 100 = 62.5 %
```

**How it is written in the code** (`dna_analysis.py`):

```python
counts = Counter(dna)
gc_content = round((counts["G"] + counts["C"]) / len(dna) * 100, 4)
```

`raw_dna_pipeline.py` uses the same formula and prints the result to 2 decimal places.

**Why it matters:** a G–C pair is joined by 3 hydrogen bonds and an A–T pair by 2, so GC-rich DNA is held together more strongly. The three repeat motifs in this project have very different GC content:

| Motif | Bases | G + C | GC content |
| :---: | :---: | :---: | :---: |
| CAG | C, A, G | 2 of 3 | 66.7 % |
| CGG | C, G, G | 3 of 3 | 100 % |
| GAA | G, A, A | 1 of 3 | 33.3 % |

---

## 📊 Key results

**Locus-level summary** (from the `*_Platinum_by_locus.csv` files)

| Motif | Disease link | Loci | Mean repeat units | Median | Max repeat units | Stable | Polymorphic |
| :---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| **CAG** | Huntington's disease | 11,220 | 5.71 | 5.0 | 51 | 10,492 (93.5%) | 728 (6.5%) |
| **CGG** | Fragile X syndrome | 5,963 | 6.83 | 6.0 | 67 | 4,177 (70.1%) | 1,786 (30.0%) |
| **GAA** | Friedreich's ataxia | 16,105 | 7.22 | 5.5 | 151 | 11,086 (68.8%) | 5,019 (31.2%) |

**Observations from this data**

- **GAA** has the most loci, the longest repeat tract observed (151 units) and the highest share of polymorphic loci (31.2%).
- **CAG** repeats are the most stable: only 6.5% of CAG loci are polymorphic.
- **Genomic location differs a lot by motif:**
  - **CGG**: 35.0% of loci are in the 5′UTR and 25.7% are exonic, far more than CAG (1.7% 5′UTR) or GAA (0.2% 5′UTR).
  - **GAA**: 53.0% intergenic and 34.5% intronic, with only 4.2% exonic.
  - **CAG**: 40.7% intergenic, 34.9% intronic and 13.0% exonic.

<details>
<summary><b>Per-allele summary</b> (from the filtered per-sample call files)</summary>

| Motif | Records | Unique loci | Mean allele units | Median | Max |
| :---: | ---: | ---: | ---: | ---: | ---: |
| CAG | 67,133 | 11,220 | 5.70 | 5 | 51 |
| CGG | 35,753 | 5,963 | 6.84 | 6 | 67 |
| GAA | 96,079 | 16,105 | 7.16 | 5 | 151 |

</details>

---

## 📈 Result graphs

All graphs are generated by the scripts and saved in the [`results/`](results/) folder.

### 1. Stable vs polymorphic loci

Number of repeat loci in each motif that are **stable** (same length in every genome) or **polymorphic** (length varies between genomes).

![Stable vs polymorphic loci by motif](results/repeat_status_comparison.png)

GAA has the most polymorphic loci (5,019) and CAG has the fewest (728).

### 2. Where each motif sits in the genome

Percentage of each motif's loci found in each genomic region.

![Genomic region distribution of CAG, CGG and GAA](results/genomic_region_comparison.png)

CGG loci are mostly in the 5′UTR (35.0%) and exons (25.7%). CAG and GAA loci are mostly intergenic and intronic.

### 3. Repeat-length distributions

How many alleles have each repeat length (number of repeat units per allele).

| CAG | CGG | GAA |
| :---: | :---: | :---: |
| <img src="results/CAG_distribution.png" width="300" alt="CAG repeat-unit distribution"> | <img src="results/CGG_distribution.png" width="300" alt="CGG repeat-unit distribution"> | <img src="results/GAA_distribution.png" width="300" alt="GAA repeat-unit distribution"> |

Most alleles are short. The most common length is 4 units for CAG (42,440 alleles) and for GAA (50,469 alleles), and 6 units for CGG (16,152 alleles). About 96.5% of CAG, 90.5% of CGG and 91.6% of GAA alleles have 10 units or fewer. GAA has a long thin tail that reaches 151 units.

### 4. Motif comparison

| Records by motif | Mean repeat units | Maximum repeat units |
| :---: | :---: | :---: |
| <img src="results/comparison_record_counts.png" width="300" alt="Records by repeat motif"> | <img src="results/comparison_mean_repeat_units.png" width="300" alt="Mean allele repeat units by motif"> | <img src="results/comparison_max_repeat_units.png" width="300" alt="Maximum allele repeat units by motif"> |

- **Records:** CAG 67,133 · CGG 35,753 · GAA 96,079
- **Mean allele repeat units:** CAG 5.70 · CGG 6.84 · GAA 7.16
- **Maximum allele repeat units:** CAG 51 · CGG 67 · GAA 151

### 5. Number of loci per motif

![Number of genomic loci by repeat motif](results/locus_count_comparison.png)

CAG 11,220 loci · CGG 5,963 loci · GAA 16,105 loci.

---

## 🗂 Project structure

```text
DNA-Triplet-Repeat-Detection/
├── dna_analysis.py            # Reusable sequence functions (validated by the tests)
├── raw_dna_pipeline.py        # Interactive pipeline for a raw DNA sequence
├── analyze_dataset.py         # Per-motif summaries and distribution plots
├── compare_motifs.py          # CAG vs CGG vs GAA comparison charts
├── analyze_locus_data.py      # Locus-level statistics
├── plot_locus_comparison.py   # Charts from the locus-level summary
├── analyze_repeat_status.py   # Stable vs polymorphic analysis
├── analyze_genomic_regions.py # Genomic region distribution
├── validate_dataset.py        # Dataset schema and quality checks
├── final_data_validation.py   # Final validation report
├── verify_results.py          # Sanity checks on the filtered data
├── genotype_anomaly_report.py # Genotype anomaly report
├── inspect_*.py, check_*.py   # Small inspection and debugging helpers
├── data/                      # Dryad dataset files + original dataset README
├── results/                   # Generated CSV reports and PNG charts
├── tests/
│   └── test_dna_analysis.py   # pytest unit tests
└── requirements.txt
```

---

## 🚀 Getting started

**Requirements:** Python 3.10 or later.

```bash
# 1. Clone the repository
git clone https://github.com/karthikeyayadavdasari-45/DNA-Triplet-Repeat-Detection.git
cd DNA-Triplet-Repeat-Detection

# 2. (Optional) create a virtual environment
python -m venv venv
venv\Scripts\activate          # Windows
# source venv/bin/activate     # macOS / Linux

# 3. Install dependencies
pip install -r requirements.txt
```

### Run the sequence analyser

```bash
python dna_analysis.py          # quick demo with the example sequence
python raw_dna_pipeline.py      # interactive: type a sequence, or press Enter for the example
```

Results are written to `results/raw_dna_summary.txt` and `results/raw_triplet_repeat_results.csv`.

### Run the dataset analysis

Run these from the project root, because some scripts use relative paths to `data/` and `results/`.

```bash
python validate_dataset.py          # check the data first
python analyze_dataset.py           # per-motif summaries + distribution plots
python compare_motifs.py            # motif comparison charts
python analyze_locus_data.py        # locus-level summary
python plot_locus_comparison.py     # needs analyze_locus_data.py to run first
python analyze_repeat_status.py     # stable vs polymorphic
python analyze_genomic_regions.py   # genomic region distribution
```

### Run the tests

```bash
pytest
```

The tests cover input normalisation, rejection of empty or invalid sequences, nucleotide and GC statistics, and repeat detection.

---

## 🧰 Scripts at a glance

| Script | Output (in `results/`) |
| :--- | :--- |
| `raw_dna_pipeline.py` | `raw_dna_summary.txt`, `raw_triplet_repeat_results.csv` |
| `analyze_dataset.py` | `{CAG,CGG,GAA}_summary.csv`, `_repeat_distribution.csv`, `_distribution.png`, `all_motifs_summary.csv` |
| `compare_motifs.py` | `final_motif_comparison.csv`, `comparison_*.png` |
| `analyze_locus_data.py` | `locus_analysis_summary.csv`, `combined_locus_data.csv` |
| `plot_locus_comparison.py` | `locus_count_comparison.png`, `mean_repeat_units_comparison.png`, `maximum_repeat_units_comparison.png` |
| `analyze_repeat_status.py` | `repeat_status_analysis.csv`, `repeat_status_comparison.png` |
| `analyze_genomic_regions.py` | `genomic_region_counts.csv`, `genomic_region_percentages.csv`, `genomic_region_comparison.png` |
| `final_data_validation.py` | `final_data_validation.csv` |
| `genotype_anomaly_report.py` | `genotype_anomaly_report.csv` |

---

## ✅ Data validation findings

The validation scripts check each motif's filtered call file for missing values, duplicate sample-locus rows and negative repeat counts.

- **No** missing values, **no** duplicate sample-locus rows and **no** negative repeat counts in any of the three files.
- **Anomaly flagged:** a small number of rows have a single-value genotype (`1`, or `0` in a few GAA rows) instead of the usual `a/b` notation such as `1/1` or `0/1`, and every one of those rows also has `Allele2_Units = 0`:

  | Motif | Rows | Numeric-GT rows with zero second allele |
  | :---: | ---: | ---: |
  | CAG | 67,133 | 156 |
  | CGG | 35,753 | 9 |
  | GAA | 96,079 | 480 |

  These rows are reported in `results/genotype_anomaly_report.csv`. They are why the allele-level minimum is `0`, while the locus-level minimum is `1`.
- The filtered call files contain calls for **6 sample IDs** (`NA12877`, `NA12878`, `NA12889`–`NA12892`). The dataset documentation describes the full Platinum cohort as 191 genomes, so per-allele numbers above reflect these 6 samples.
- **Tip for reproducing:** read `Call_ID` as text (`dtype={"Call_ID": str}`). It contains values like `1.100104164`, and reading it with mixed types can miscount unique loci.

---

## 📚 Dataset and credits

Data: **CGG, CAG, and GAA: genome-wide comparison of the disease linked trinucleotide short tandem repeat**, Dryad, DOI [10.5061/dryad.5tb2rbpgt](https://doi.org/10.5061/dryad.5tb2rbpgt).

The dataset was derived from the [MSSNG Project](https://research.mss.ng/) and [Illumina Platinum Genomes](https://emea.illumina.com/platinumgenomes.html). The original dataset documentation is kept in [`data/README.md`](data/README.md). Please follow the dataset's own terms of use on the Dryad page when reusing the data.

The files in `data/` come from the Dryad dataset above, and all credit for the data belongs to its original authors. The analysis code, tests and results in this repository are our own work.

---

## 🔭 Limitations and future work

**Current limitations**
- The raw sequence analyser reads **reading frame 0 only** and detects **consecutive identical triplets**. A repeat that starts off-frame can be reported as a different motif or missed.
- It processes **one sequence at a time** from the command line.
- Dataset analysis uses the **precomputed repeat calls** from Dryad. It does not call repeats from raw sequencing reads.

**Ideas for next steps**
- Scan all three reading frames and the reverse complement
- Accept FASTA files and compare several sequences at once
- Add a simple web or desktop interface
- Add disease-threshold flags per motif and gene-level summaries

---

## 👤 Author

- [@karthikeyayadavdasari-45](https://github.com/karthikeyayadavdasari-45)

Amrita Vishwa Vidyapeetham, Coimbatore
