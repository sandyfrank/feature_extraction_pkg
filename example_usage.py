"""
Example script demonstrating how to use the Protein Feature Extraction package.
"""

from unified_extractor import extract_all_features
import os

def main():
    """Main example function."""
    
    # Example 1: Basic usage with minimal parameters
    print("=" * 70)
    print("Example 1: Basic Feature Extraction")
    print("=" * 70)
    
    # Assuming you have a proteins.fasta file
    fasta_file = "proteins.fasta"
    output_file = "protein_features.csv"
    
    if os.path.exists(fasta_file):
        df = extract_all_features(
            fasta_file=fasta_file,
            output_file=output_file,
            verbose=True
        )
        
        print(f"\nResults:")
        print(f"  - Shape: {df.shape}")
        print(f"  - Index (Protein IDs): {df.index.tolist()[:3]}...")
        print(f"  - Columns (Features): {df.columns.tolist()[:5]}...")
    else:
        print(f"Example FASTA file '{fasta_file}' not found.")
        print("Please provide a valid FASTA file path.\n")
    
    # Example 2: Advanced usage with all parameters
    print("\n" + "=" * 70)
    print("Example 2: Advanced Feature Extraction with ESM")
    print("=" * 70)
    
    advanced_example = """
    from unified_extractor import extract_all_features
    
    # Extract features with ESM embeddings
    df = extract_all_features(
        fasta_file="proteins.fasta",
        output_file="features_with_esm.csv",
        word2vec_model_path="/path/to/word2vec_model.txt",
        esm_model_name="esm2_t33_650M_UR50D",
        esm_tokens_per_batch=4096,
        esm_seq_length=1022,
        esm_repr_layers=[33],
        verbose=True
    )
    
    # Access the features
    print(f"Total features extracted: {df.shape[1]}")
    print(f"First protein: {df.index[0]}")
    print(f"First 5 feature values: {df.iloc[0, :5].values}")
    """
    
    print(advanced_example)
    
    # Example 3: Feature processing for machine learning
    print("\n" + "=" * 70)
    print("Example 3: Preparing Features for Machine Learning")
    print("=" * 70)
    
    ml_example = """
    from unified_extractor import extract_all_features
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import RandomForestClassifier
    
    # Extract features
    X = extract_all_features("proteins.fasta").values
    
    # Example labels (you would provide these)
    y = [0, 1, 0, 1, ...]  # Binary classification
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train classifier
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X_train_scaled, y_train)
    
    # Evaluate
    accuracy = clf.score(X_test_scaled, y_test)
    print(f"Accuracy: {accuracy:.2f}")
    """
    
    print(ml_example)
    
    # Example 4: Feature inspection
    print("\n" + "=" * 70)
    print("Example 4: Inspecting Extracted Features")
    print("=" * 70)
    
    inspection_example = """
    from unified_extractor import extract_all_features
    import pandas as pd
    
    # Extract features
    df = extract_all_features("proteins.fasta", "features.csv")
    
    # Inspect feature groups
    aac_features = df.filter(regex='^G1_')  # Amino Acid Composition
    dcp_features = df.filter(regex='^G2_')  # Dipeptide Composition
    pcp_features = df.filter(regex='^G3_')  # Physicochemical Properties
    ctriad_features = df.filter(regex='^G4_')  # Conjoint Triad
    w2v_features = df.filter(regex='^G5_')  # Word2Vec
    entropy_features = df.filter(regex='^G6_')  # Entropy
    protparams_features = df.filter(regex='^G7_')  # Protein Parameters
    qso_features = df.filter(regex='^G8_')  # QSO
    
    print(f"AAC features: {aac_features.shape[1]}")
    print(f"DCP features: {dcp_features.shape[1]}")
    print(f"PCP features: {pcp_features.shape[1]}")
    print(f"CTriad features: {ctriad_features.shape[1]}")
    print(f"Word2Vec features: {w2v_features.shape[1]}")
    print(f"Entropy features: {entropy_features.shape[1]}")
    print(f"ProtParams features: {protparams_features.shape[1]}")
    print(f"QSO features: {qso_features.shape[1]}")
    
    # Statistical summary
    print("\\nFeature Statistics:")
    print(df.describe())
    
    # Check for missing values
    print("\\nMissing Values:")
    print(df.isnull().sum())
    """
    
    print(inspection_example)
    
    print("\n" + "=" * 70)
    print("For more information, see README.md")
    print("=" * 70)


if __name__ == "__main__":
    main()
