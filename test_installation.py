#!/usr/bin/env python
"""
Quick test script to verify the package installation and functionality.
"""

import sys
import os

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    try:
        import pandas as pd
        print("  ✓ pandas")
        import numpy as np
        print("  ✓ numpy")
        from Bio import SeqIO
        print("  ✓ Biopython")
        import torch
        print("  ✓ torch")
    except ImportError as e:
        print(f"  ✗ Import error: {e}")
        return False
    return True

def test_features_extraction():
    """Test basic feature extraction functions."""
    print("\nTesting feature extraction modules...")
    try:
        from Features_Extraction.FeaturesExtractionModule import (
            extract_aac_features,
            extract_dcp_features,
            extract_pcp_features,
            extract_ctriad_features,
            extract_entropy_features,
            extract_protparams_features,
            extract_qso_features,
            extract_word2vec,
        )
        print("  ✓ All feature extraction functions imported successfully")
        return True
    except ImportError as e:
        print(f"  ✗ Import error: {e}")
        return False

def test_unified_extractor():
    """Test the unified extractor module."""
    print("\nTesting unified extractor...")
    try:
        from unified_extractor import extract_all_features
        print("  ✓ Unified extractor imported successfully")
        return True
    except ImportError as e:
        print(f"  ✗ Import error: {e}")
        return False

def main():
    """Run all tests."""
    print("=" * 70)
    print("PROTEIN FEATURE EXTRACTION - INSTALLATION TEST")
    print("=" * 70)
    
    # Change to package directory
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(pkg_dir)
    
    tests = [
        test_imports,
        test_features_extraction,
        test_unified_extractor,
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"  ✗ Error: {e}")
            results.append(False)
    
    print("\n" + "=" * 70)
    if all(results):
        print("✓ All tests passed! Package is ready to use.")
        print("\nUsage:")
        print("  from protein_feature_extraction import extract_all_features")
        print("  df = extract_all_features('proteins.fasta', 'output.csv')")
        return 0
    else:
        print("✗ Some tests failed. Please check dependencies.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
