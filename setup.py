from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="protein-feature-extraction",
    version="1.0.0",
    author="Sandy-Frank NGAHA, Nelle VAROQUAUX, and Sophie ABBY",
    author_email="sandy.ngaha@aims-cameroon.org",
    description=(
        "A comprehensive package for extracting 9 types of protein features "
        "from FASTA files. Includes Word2Vec embeddings trained on NCBI 2023 "
        "protein database (~37,000 genomes, 80+ million proteins)"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sandyfrank/feature_extraction_pkg.git",

    # ── Package discovery ────────────────────────────────────────────────
    # The repo root IS the package (feature_extraction_pkg/__init__.py exists).
    # We expose it under the name 'protein_feature_extraction' so that
    # the entry point and any `import protein_feature_extraction` work.
    packages=["protein_feature_extraction", "protein_feature_extraction.Features_Extraction"],
    package_dir={
        "protein_feature_extraction":                    ".",
        "protein_feature_extraction.Features_Extraction": "Features_Extraction",
    },

    package_data={
        "protein_feature_extraction": [
            "example_proteins.fasta",
            "word2vec_train.txt",
        ],
    },
    include_package_data=True,

    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "pandas>=1.0.0",
        "numpy>=1.19.0",
        "biopython>=1.78",
        "scipy>=1.5.0",
        "gensim>=4.0.0",
        "torch>=1.9.0",
        "fair-esm",
    ],
    extras_require={
        "dev": ["pytest", "pytest-cov"],
    },

    # ── CLI entry point ──────────────────────────────────────────────────
    # After pip install: `extract-proteins proteins.fasta output.csv`
    entry_points={
        "console_scripts": [
            "extract-proteins=protein_feature_extraction.unified_extractor:main",
        ],
    },
)