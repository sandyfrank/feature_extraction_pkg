"""
Protein Feature Extraction Package

A comprehensive package for extracting 9 different types of features from protein sequences.
"""

__version__ = "1.0.0"
__author__ = "Feature Extraction Team"

from .unified_extractor import extract_all_features, main

__all__ = ['extract_all_features', 'main']
