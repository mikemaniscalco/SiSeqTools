# SiSeqTools

Tools for analyzing silicon-related gene and protein sequences.  
This repository is a growing collection of lightweight, command-line utilities for sequence-based analysis.  

Currently, the repo contains **`sliding_window_aa.py`**, a tool for scanning protein sequences with a sliding window to measure amino acid composition.  

---

## 🔬 Current Tool: `sliding_window_aa.py`

### Description
`sliding_window_aa.py` scans protein FASTA sequences with a sliding window to measure amino acid composition (currently **% Lysine** and **% Serine**).  
Windows that meet one of the following thresholds are reported:

- ≥ **10% Lysine**  
- ≥ **18% Serine**

To reduce redundant output, overlapping windows are consolidated:

1. Windows meeting **both thresholds** are prioritized over single-threshold windows.  
2. **Longer windows** are favored over shorter ones.  
3. Overlaps of ≤25% are allowed, but redundant identical windows are collapsed.  

### Usage
```bash
python sliding_window_aa.py -i input.fasta -w 100-2000 -o output.tsv -s 1 -t 2
````

### Arguments

| Argument           | Description                                                                                |
| ------------------ | ------------------------------------------------------------------------------------------ |
| `-i` / `--input`   | Input FASTA file of amino acid sequences                                                   |
| `-o` / `--output`  | Output file (TSV format)                                                                   |
| `-w` / `--window`  | Window size: either a single integer (e.g., `100`) or a range `MIN-MAX` (e.g., `100-2000`) |
| `-s` / `--step`    | Step size for sliding window (**default: 5**)                                              |
| `-t` / `--threads` | Number of CPU threads (**default: 1**)                                                     |

### Example Output (`output.tsv`)

```
seq_id   window   start   end   percent_K   percent_S
seq1     150      0       150   12.0000     19.3333
seq1     300      50      350   11.5000     20.0000
seq2     200      10      210   10.2000     17.8000
```

---

## 📦 Installation

Clone this repository:

```bash
git clone https://github.com/<your-username>/SiSeqTools.git
cd SiSeqTools
```

Install required dependencies (Biopython):

```bash
pip install biopython
```

---

## ⚙️ Requirements

Python 3.8+
(Script is compatible with Python 3.6+, but Biopython and other modern packages now require Python ≥3.8.)

Biopython (tested with version ≥1.83)
---

## 🛠 Planned Additions

* Additional amino acid scanning utilities
* Motif/peptide finder tools
* Visualization helpers for sliding window outputs
* Sequence filtering and annotation utilities

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
