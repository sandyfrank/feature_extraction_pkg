# Protein Feature Extraction Package

A comprehensive Python package for extracting multiple types of features from protein sequences in FASTA format. This package combines 9 different feature extraction methods and produces a single unified CSV file with all features.

## Features Extracted

The package extracts **9 different types of protein features**:

1. **Amino Acid Composition (AAC)** - G1 (20 features)
   - Normalized frequency of each amino acid

2. **Dipeptide Composition (DCP)** - G2 (400 features)
   - Normalized frequency of dipeptide pairs

3. **Physicochemical Properties (PCP)** - G3 (42 features)
   - 21 physicochemical properties + 21 transition properties

4. **Conjoint Triad Composition (CTriad)** - G4 (343 features)
   - 7×7×7 triad matrix flattened

5. **Entropy Features** - G6 (10 features)
   - Shannon entropy calculated over 10 sequence segments

6. **Protein Parameters** - G7 (7 features)
   - Molecular weight, aromaticity, instability index, GRAVY, isoelectric point, molar extinction coefficients

7. **Quasi-Sequence Order (QSO)** - G8 (50 features)
   - Sequence-order coupling numbers with weights

8. **Word2Vec Features** - G5 (100 features)
   - Word2Vec embeddings based on 3-grams
   - Trained on NCBI 2023 protein database (~37,000 genomes, 80+ million proteins)

9. **ESM Embeddings** - G9 (1280 features)
   - Pre-trained ESM language model embeddings from layer 33

---

## Installation

### Requirements

- Python 3.7+
- pandas, numpy, biopython, scipy, gensim, torch, fair-esm

### Install from GitHub

```bash
pip3 install git+https://github.com/sandyfrank/feature_extraction_pkg.git
```

This installs all dependencies and registers the `extract-proteins` command.

---

## Usage

### Command Line

```bash
# Basic usage
extract-proteins proteins.fasta output_features.csv

# Quick test on the included example file (5 proteins, ~1-2 min)
extract-proteins example_proteins.fasta example_output.csv
```

### Python API

```python
from protein_feature_extraction import extract_all_features

# Extract all features from a FASTA file
features_df = extract_all_features("proteins.fasta", "output_features.csv")

print(f"Sequences    : {features_df.shape[0]}")
print(f"Total features: {features_df.shape[1]}")
```

### Advanced Usage

```python
from protein_feature_extraction import extract_all_features

features_df = extract_all_features(
    fasta_file="proteins.fasta",
    output_file="all_features.csv",
    word2vec_model_path="/path/to/word2vec_model.txt",
    esm_model_name="esm2_t33_650M_UR50D",
    esm_tokens_per_batch=4096,
    esm_seq_length=1022,
    esm_repr_layers=[33],
    verbose=True,
)
```

---

## Parameters

### `extract_all_features()`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fasta_file` | str | required | Path to input FASTA file |
| `output_file` | str | `"all_features.csv"` | Path to save output CSV |
| `word2vec_model_path` | str | None | Path to Word2Vec model file |
| `esm_model_name` | str | `"esm2_t33_650M_UR50D"` | ESM model name |
| `esm_output_dir` | str | None | Directory to save ESM embeddings |
| `esm_tokens_per_batch` | int | `4096` | Tokens per batch for ESM |
| `esm_seq_length` | int | `1022` | Maximum sequence length |
| `esm_repr_layers` | list | `[33]` | ESM representation layers |
| `temp_dir` | str | None | Temporary directory for intermediate files |
| `verbose` | bool | `True` | Print progress messages |

---

## Output

The output CSV file contains:
- **Index**: Protein sequence IDs from the FASTA file
- **Columns**: Combined features from all 9 extraction methods

| Prefix | Feature type | # Features |
|---|---|---|
| `G1_*` | Amino Acid Composition | 20 |
| `G2_*` | Dipeptide Composition | 400 |
| `G3_*` | Physicochemical Properties | 42 |
| `G4_*` | Conjoint Triad Composition | 343 |
| `G5_*` | Word2Vec Features | 100 |
| `G6_*` | Entropy Features | 10 |
| `G7_*` | Protein Parameters | 7 |
| `G8_*` | QSO Features | 50 |
| `G9_*` | ESM Embeddings | 1280 |
| | **Total** | **2252** |

---

## Examples

### Example 1: Basic Feature Extraction

```python
from protein_feature_extraction import extract_all_features

df = extract_all_features("proteins.fasta", "features.csv")
print(df.head())
```

### Example 2: Machine Learning Pipeline

```python
from protein_feature_extraction import extract_all_features
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier

# Extract features
X = extract_all_features("proteins.fasta").values

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train classifier
clf = RandomForestClassifier()
clf.fit(X_scaled, y)
```

---

## Performance Notes

| Dataset size | CPU | GPU |
|---|---|---|
| 5 proteins | ~1–2 min | ~1–2 min |
| 270 proteins | ~10–15 min | ~5–10 min |
| 1000+ proteins | 1+ hour | ~30 min |

ESM embeddings dominate runtime — GPU is strongly recommended for large datasets.

Memory usage: ~1–2 MB per sequence for basic features; more for ESM embeddings.

---

## Troubleshooting

**Word2Vec extraction fails**
Provide the correct model path via `word2vec_model_path`, or ensure the default model exists at the expected location.

**ESM extraction fails**
Ensure `fair-esm` and `torch` are installed:
```bash
pip3 install torch fair-esm
```

**Out of memory with ESM**
Use a smaller ESM model, reduce `esm_tokens_per_batch`, or increase available memory.

---

## Citation

If you use this package in your research, please cite the original feature extraction methods:

- Amino Acid Composition: Chou, K. C. (1995)
- Physicochemical Properties: Karchin, R. et al. (2006)
- Conjoint Triad: Mishra, S. et al. (2010)
- Quasi-Sequence Order: Chou, K. C. (1997)
- Word2Vec: Mikolov, T. et al. (2013)
- ESM: Rives, A. et al. (2021)

---

## License

This project is licensed under the MIT License — see the LICENSE file for details.