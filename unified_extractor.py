"""
Unified Feature Extraction Module

This module provides a single function to extract all 9 types of protein features
from a FASTA file and combine them into a single CSV.
"""

import os
import tempfile
import pandas as pd
import pathlib
from typing import Optional

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
from Features_Extraction.extract_esm import extract_esm_embeddings


def extract_all_features(
    fasta_file: str,
    output_file: str = "all_features.csv",
    word2vec_model_path: Optional[str] = None,
    esm_model_name: str = "esm2_t33_650M_UR50D",
    esm_output_dir: Optional[str] = None,
    esm_tokens_per_batch: int = 4096,
    esm_seq_length: int = 1022,
    esm_repr_layers: list = None,
    temp_dir: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Extract all 9 types of protein features from a FASTA file and combine them into a single CSV.

    This function extracts features using the following 9 methods:
    1. Amino Acid Composition (AAC) - G1
    2. Dipeptide Composition (DCP) - G2
    3. Physicochemical Properties (PCP) - G3
    4. Conjoint Triad Composition (CTriad) - G4
    5. Entropy Features - G6
    6. Protein Parameters - G7
    7. Quasi-Sequence Order (QSO) - G8
    8. Word2Vec Features - G5
    9. ESM Embeddings (from pre-trained ESM model) - G9

    Args:
        fasta_file (str): Path to input FASTA file containing protein sequences.
        output_file (str): Path to save the combined features CSV. Defaults to "all_features.csv".
        word2vec_model_path (str, optional): Path to Word2Vec model file. If None, uses default path.
        esm_model_name (str): ESM model name to use. Defaults to "esm2_t33_650M_UR50D".
        esm_output_dir (str, optional): Directory to save ESM embeddings. Defaults to temp directory.
        esm_tokens_per_batch (int): Number of tokens per batch for ESM. Defaults to 4096.
        esm_seq_length (int): Maximum sequence length for ESM. Defaults to 1022.
        esm_repr_layers (list, optional): ESM representation layers to extract. Defaults to [33].
        temp_dir (str, optional): Temporary directory for intermediate CSV files. Defaults to system temp.
        verbose (bool): Print progress messages. Defaults to True.

    Returns:
        pd.DataFrame: DataFrame with all combined features, also saved to output_file.

    Example:
        >>> from protein_feature_extraction import extract_all_features
        >>> features_df = extract_all_features("proteins.fasta", "all_features.csv")
        >>> print(features_df.shape)
        (100, 1024)  # 100 proteins with 1024 combined features
    """
    
    if esm_repr_layers is None:
        esm_repr_layers = [33]
    
    # Create temporary directory for intermediate files
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp(prefix="protein_features_")
    else:
        os.makedirs(temp_dir, exist_ok=True)
    
    # Validate input file
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    fasta_file = os.path.abspath(fasta_file)
    output_file = os.path.abspath(output_file)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    
    try:
        # Set the global outname for extraction functions
        import Features_Extraction.FeaturesExtractionModule as fem
        fem.outname = temp_dir
        
        if verbose:
            print("=" * 70)
            print("PROTEIN FEATURE EXTRACTION")
            print("=" * 70)
            print(f"Input FASTA: {fasta_file}")
            print(f"Output file: {output_file}")
            print(f"Temp directory: {temp_dir}")
            print("=" * 70)
        
        # Extract all features
        all_features = None
        feature_types = []
        
        if verbose:
            print("\n[1/9] Extracting Amino Acid Composition (AAC)...")
        extract_aac_features(fasta_file, "features_aac.csv")
        aac_df = pd.read_csv(os.path.join(temp_dir, "features_aac.csv"), index_col=0)
        all_features = aac_df
        feature_types.append("AAC")
        
        if verbose:
            print("[2/9] Extracting Dipeptide Composition (DCP)...")
        extract_dcp_features(fasta_file, "features_dcp.csv")
        dcp_df = pd.read_csv(os.path.join(temp_dir, "features_dcp.csv"), index_col=0)
        all_features = all_features.join(dcp_df)
        feature_types.append("DCP")
        
        if verbose:
            print("[3/9] Extracting Physicochemical Properties (PCP)...")
        extract_pcp_features(fasta_file, "features_pcp.csv")
        pcp_df = pd.read_csv(os.path.join(temp_dir, "features_pcp.csv"), index_col=0)
        all_features = all_features.join(pcp_df)
        feature_types.append("PCP")
        
        if verbose:
            print("[4/9] Extracting Conjoint Triad Composition (CTriad)...")
        extract_ctriad_features(fasta_file, "features_ctriad.csv")
        ctriad_df = pd.read_csv(os.path.join(temp_dir, "features_ctriad.csv"), index_col=0)
        all_features = all_features.join(ctriad_df)
        feature_types.append("CTriad")
        
        if verbose:
            print("[5/9] Extracting Entropy Features...")
        extract_entropy_features(fasta_file, "features_entropy.csv")
        entropy_df = pd.read_csv(os.path.join(temp_dir, "features_entropy.csv"), index_col=0)
        all_features = all_features.join(entropy_df)
        feature_types.append("Entropy")
        
        if verbose:
            print("[6/9] Extracting Protein Parameters...")
        extract_protparams_features(fasta_file, "features_protparams.csv")
        protparams_df = pd.read_csv(os.path.join(temp_dir, "features_protparams.csv"), index_col=0)
        all_features = all_features.join(protparams_df)
        feature_types.append("ProtParams")
        
        if verbose:
            print("[7/9] Extracting Quasi-Sequence Order (QSO)...")
        extract_qso_features(fasta_file, "features_qso.csv")
        qso_df = pd.read_csv(os.path.join(temp_dir, "features_qso.csv"), index_col=0)
        all_features = all_features.join(qso_df)
        feature_types.append("QSO")
        
        if verbose:
            print("[8/9] Extracting Word2Vec Features...")
        try:
            # Determine Word2Vec model path: use provided path or bundled model
            model_path = word2vec_model_path
            if model_path is None:
                bundled = pathlib.Path(__file__).parent / "Features_Extraction" / "word2vec_train.txt"
                if bundled.exists():
                    model_path = str(bundled)
                else:
                    raise FileNotFoundError("No Word2Vec model provided and bundled model not found")

            extract_word2vec(fasta_file, "features_word2vec.csv", model_filename=model_path)
            w2v_df = pd.read_csv(os.path.join(temp_dir, "features_word2vec.csv"), index_col=0)
            all_features = all_features.join(w2v_df)
            feature_types.append("Word2Vec")
        except Exception as e:
            if verbose:
                print(f"  Warning: Word2Vec extraction skipped - {str(e)}")
        
        if verbose:
            print("[9/9] Extracting ESM Embeddings...")
        
        if esm_output_dir is None:
            esm_output_dir = os.path.join(temp_dir, "esm_embeddings")
        
        extract_esm_embeddings(
            esm_model_name,
            pathlib.Path(fasta_file),
            pathlib.Path(esm_output_dir),
            tokens_per_batch=esm_tokens_per_batch,
            seq_length=esm_seq_length,
            repr_layers=esm_repr_layers
        )
        
        # Load ESM embeddings
        esm_features = _load_esm_embeddings(esm_output_dir, repr_layers=esm_repr_layers)
        if esm_features is not None and len(esm_features) > 0:
            esm_df = pd.DataFrame.from_dict(esm_features, orient='index')
            # Rename columns with G9_ prefix
            esm_df.columns = [f"G9_{i}" for i in range(esm_df.shape[1])]
            # Align with existing sequences
            all_features = all_features.join(esm_df)
            feature_types.append("ESM")
        
        # Save combined features
        all_features.to_csv(output_file, index=True)
        
        if verbose:
            print("=" * 70)
            print(f"✓ Feature extraction complete!")
            print(f"  - Total sequences: {all_features.shape[0]}")
            print(f"  - Total features: {all_features.shape[1]}")
            print(f"  - Feature types extracted: {', '.join(feature_types)}")
            print(f"  - Output saved to: {output_file}")
            print("=" * 70)
        
        return all_features
        
    finally:
        # Cleanup temporary files if not explicitly requested to keep
        if verbose:
            print(f"\nTemporary files in: {temp_dir}")


def _load_esm_embeddings(esm_output_dir: str, repr_layers: list = None) -> dict:
    """
    Load ESM embeddings from PyTorch files.
    
    Args:
        esm_output_dir: Directory containing .pt files with ESM embeddings
        repr_layers: List of representation layers to extract
        
    Returns:
        Dictionary mapping sequence IDs to flattened embeddings
    """
    import torch
    
    if repr_layers is None:
        repr_layers = [33]
    
    esm_features = {}
    esm_output_dir = pathlib.Path(esm_output_dir)
    
    if not esm_output_dir.exists():
        return esm_features
    
    for pt_file in esm_output_dir.glob("*.pt"):
        try:
            data = torch.load(str(pt_file))
            entry_id = data["entry_id"]
            
            # Extract mean representations from specified layers
            mean_reps = data.get("mean_representations", {})
            if mean_reps:
                # Concatenate representations from all layers
                embeddings = []
                for layer in sorted(mean_reps.keys()):
                    if layer in repr_layers or not repr_layers:
                        rep = mean_reps[layer]
                        if isinstance(rep, torch.Tensor):
                            embeddings.append(rep.cpu().numpy().flatten())
                        else:
                            embeddings.append(rep.flatten())
                
                if embeddings:
                    import numpy as np
                    esm_features[entry_id] = np.concatenate(embeddings)
        except Exception as e:
            print(f"Warning: Failed to load {pt_file}: {e}")
    
    return esm_features


def main():
    """Command-line interface for protein feature extraction."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: extract-proteins <fasta_file> [output_file]")
        print("\nExample:")
        print("  extract-proteins proteins.fasta features.csv")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "all_features.csv"
    
    try:
        extract_all_features(fasta_file, output_file, verbose=True)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
