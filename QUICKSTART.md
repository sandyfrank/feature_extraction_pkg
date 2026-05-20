# Quick Start Guide

## Installation

### 1. Install Dependencies

```bash
# Install all required packages
pip install pandas numpy biopython scipy gensim torch fair-esm
```

### 2. Verify Installation

```bash
python test_installation.py
```

You should see:
```
✓ All tests passed! Package is ready to use.
```

## Quick Test

The package includes `example_proteins.fasta` with 5 test proteins. Try it now:

```bash
python unified_extractor.py example_proteins.fasta example_output.csv
```

**Expected runtime**: ~1-2 minutes

## Basic Usage

### Simple Example (Python)

```python
from unified_extractor import extract_all_features

# Extract all features from proteins.fasta
features_df = extract_all_features(
    fasta_file="proteins.fasta",
    output_file="all_features.csv"
)

print(f"Extracted features for {features_df.shape[0]} proteins")
print(f"Total features: {features_df.shape[1]}")
```

### Command Line Usage

```bash
python unified_extractor.py proteins.fasta all_features.csv
```

## Package Structure

```
feature_extraction_pkg/
├── __init__.py                      # Package initialization
├── unified_extractor.py             # Main unified extraction function
├── example_proteins.fasta           # Example FASTA file (5 proteins, ~1-2 min to process)
├── README.md                        # Full documentation
├── LICENSE                          # MIT License
├── setup.py                         # Setup configuration
├── pyproject.toml                   # Modern Python project config
├── example_usage.py                 # Usage examples
├── test_installation.py             # Installation test script
├── QUICKSTART.md                    # This file
└── Features_Extraction/             # Original extraction modules
    ├── FeaturesExtractionModule.py  # 8 extraction methods
    ├── extract_esm.py               # ESM embeddings
    ├── utils.py                     # Utility functions
    ├── extract_esm_sequences.py     # ESM sequence extraction
    ├── extract_esm_features.sh      # ESM bash script
    ├── joint_features.py            # Join CSV files
    └── word2vec_train.txt           # Word2Vec model (if available)
```

## Features Overview

| Feature Type | Prefix | # Features | Description |
|---|---|---|---|
| Amino Acid Composition | G1_ | 20 | Normalized frequency of each amino acid |
| Dipeptide Composition | G2_ | 400 | Normalized frequency of dipeptide pairs |
| Physicochemical Properties | G3_ | 42 | Properties and transitions |
| Conjoint Triad | G4_ | 343 | 7×7×7 triad matrix |
| Word2Vec | G5_ | 100 | Word2Vec embeddings (3-grams) |
| Entropy | G6_ | 10 | Shannon entropy per segment |
| Protein Parameters | G7_ | 7 | Physical/chemical properties |
| QSO | G8_ | 50 | Quasi-sequence order |
| ESM Embeddings | G9_ | 1280 | Pre-trained language model embeddings |

**Total: ~2,280 combined features** (1280 from ESM + ~1000 from others)

## Common Issues & Solutions

### Issue: "outname variable not defined"
**Solution**: This is automatically handled by the package. No action needed.

### Issue: Word2Vec extraction fails
**Solution**: 
```python
# Provide Word2Vec model path
features_df = extract_all_features(
    fasta_file="proteins.fasta",
    word2vec_model_path="/path/to/your/model.txt"
)
```

### Issue: ESM extraction not working
**Solution**: Ensure `fair-esm` and `torch` are installed correctly. Run `pip install torch fair-esm`. If using CPU, processing will be slower; GPU is recommended.
### Issue: Out of memory with large files
**Solution**:
```python
# Process in batches or reduce ESM batch size
features_df = extract_all_features(
    fasta_file="proteins.fasta",
    esm_tokens_per_batch=2048,  # Reduced from 4096
    temp_dir="/path/with/more/space"
)
```

### Performance
- **Small files** (5 proteins): ~1-2 minutes
- **Medium files** (100+ proteins): 10-15 minutes (CPU) or 5-10 minutes (GPU)
- **Large files** (1000+ proteins): 1+ hours (CPU), 30+ minutes (GPU)

ESM embeddings dominate the execution time. Using a GPU significantly speeds up processing.

## Example FASTA Input

```fasta
>protein_001
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ
VKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSV
>protein_002
MGSVGSGVFGFRRGLLICEHEKAQEAAQQAQQQQQQQQQQQQQQQQQQQQQQQQQ
KGKGKVGGAVGIGMYRDXDXSGYPDXKGKGKGKGK
```

## Example Output CSV Structure

```
,G1_A,G1_G,G1_I,...,G2_AA,G2_AG,...,G3_0,...,G8_QSOSW50
protein_001,0.123,0.098,0.067,...,0.045,0.032,...,0.234,...,0.876
protein_002,0.156,0.087,0.089,...,0.052,0.041,...,0.198,...,0.823
```

## Next Steps

1. **Prepare your data**: Ensure your proteins are in FASTA format
2. **Run extraction**: Use the examples above to extract features
3. **Analyze results**: Load the output CSV in your analysis tool
4. **Use for ML**: Feed the features to your machine learning pipeline

## For More Information

- See [README.md](README.md) for comprehensive documentation
- Check [example_usage.py](example_usage.py) for more code examples
- Review [MIT License](LICENSE) for license information

## Support

If you encounter issues:
1. Check this QUICKSTART guide
2. Review README.md
3. Check test_installation.py output
4. Verify all dependencies are installed with correct versions
