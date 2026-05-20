import os
import gzip
import pandas as pd
from Bio import SeqIO
from itertools import product
from contextlib import contextmanager
import collections
from scipy.stats import entropy
from gensim.models import Word2Vec
import csv
import sys
import glob
import os
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Utility functions

# Define all amino acids
all_aas = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', 'D', 'E', 'R',
           'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q']

@contextmanager
def open_file(filename, mode="rt"):
    """Handles opening normal and gzipped files."""
    if filename.endswith(".gz"):
        handle = gzip.open(filename, mode)
    else:
        handle = open(filename, mode)
    
    try:
        yield handle
    finally:
        handle.close()

# Amino Acid Composition (AAC)

def calculate_aac(sequence):
    """Calculate amino acid composition."""
    seq_length = len(sequence)
    return [sequence.count(aa) / seq_length for aa in all_aas]



def extract_aac_features(fasta_file,output_file = "features_aac.csv"):
    
    """Extract Amino Acid Composition (AAC) features from a FASTA file."""
    all_aac = []
    all_sequences_names = []
    
    output_file = os.path.join(outname,output_file)

    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            aac = calculate_aac(sequence)
            all_aac.append(aac)
            all_sequences_names.append(record.id)

    aac_df = pd.DataFrame(all_aac, index=all_sequences_names, columns=[f"G1_{aa}" for aa in all_aas])
    aac_df.to_csv(output_file)
    print(f"AAC Features saved to: {output_file}")

# Dipeptide Composition (DCP)

def calculate_dipeptide_composition(sequence):
    """Calculate dipeptide composition."""
    seq_length = len(sequence) - 1
    dipeptides = ["".join(pair) for pair in product(all_aas, repeat=2)]
    
    dipeptide_counts = {dipep: 0 for dipep in dipeptides}
    for i in range(seq_length):
        dipep = sequence[i:i+2]
        if dipep in dipeptide_counts:
            dipeptide_counts[dipep] += 1
    
    return [dipeptide_counts[dipep] / seq_length for dipep in dipeptides] if seq_length > 0 else [0] * len(dipeptides)



def extract_dcp_features(fasta_file,output_file = "features_dcp.csv"):
    
    """
    Extracts dipeptide composition (DCP) features from a given FASTA file and saves them to a CSV file.
    Args:
        fasta_file (str): Path to the input FASTA file containing sequences.
        output_file (str, optional): Name of the output CSV file where the DCP features will be saved. 
                                        Defaults to "features_dcp.csv".
    Returns:
        None
    Outputs:
        A CSV file containing the DCP features for each sequence in the FASTA file. Each row corresponds 
        to a sequence, and each column represents a dipeptide composition feature.
    Notes:
        - The function assumes that the input FASTA file contains valid sequences.
        - The output file will be saved in the directory specified by the `outname` variable.
        - The dipeptide composition is calculated using the `calculate_dipeptide_composition` function.
        - The column names in the output CSV file are prefixed with "G2_" followed by the dipeptide pair.
    Example:
        extract_dcp_features("input_sequences.fasta", "output_features.csv")
    """

    all_dcp = []
    all_sequences_names = []
    
    output_file = os.path.join(outname,output_file)

    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            dcp = calculate_dipeptide_composition(sequence)
            all_dcp.append(dcp)
            all_sequences_names.append(record.id)

    dcp_df = pd.DataFrame(all_dcp, index=all_sequences_names, columns=[f"G2_{''.join(pair)}" for pair in product(all_aas, repeat=2)])
    dcp_df.to_csv(output_file)
    print(f"DCP Features saved to: {output_file}")


# Physicochemical Properties (PCP)

def extract_pcp_features(fasta_file, output_file="features_pcp.csv"):
    """
    Extracts physicochemical and transition properties (PCP) features from protein sequences 
    in a FASTA file and saves the results to a CSV file.
    The function calculates the following:
    - Physicochemical properties (PCP) based on predefined amino acid attributes.
    - Transition properties, which count transitions between different attributes in the sequence.
    Args:
        fasta_file (str): Path to the input FASTA file containing protein sequences.
        output_file (str, optional): Name of the output CSV file to save the extracted features. 
                                     Defaults to "features_pcp.csv".
    Returns:
        None: The function saves the extracted features to the specified CSV file.
    Notes:
        - The output CSV file contains 42 columns: 21 for normalized physicochemical properties 
          and 21 for normalized transition properties.
        - The index of the CSV file corresponds to the sequence IDs from the FASTA file.
    Example:
        >>> extract_pcp_features("input.fasta", "output_features.csv")
        Processing complete. Results saved in output_features.csv
    """
    p_to_index = {
        'polar': 0, 'neutral': 1, 'apolar': 2,
        'small v': 3, 'medium v': 4, 'large v': 5,  # volume
        'low p': 6, 'medium p': 7, 'high p': 8,  # polarity
        'low pz': 9, 'medium pz': 10, 'high pz': 11,  # polarizability
        'positive': 12, 'no charge': 13, 'negative': 14,
        'helix': 15, 'strand': 16, 'coil': 17,
        'buried': 18, 'exposed': 19, 'intermediate': 20
    }

    amino_acid_map = {
        'A': ['neutral', 'small v', 'medium p', 'low pz', 'no charge', 'helix', 'buried'],
        'C': ['apolar', 'small v', 'low p', 'medium pz', 'no charge', 'strand', 'buried'],
        'D': ['polar', 'small v', 'high p', 'low pz', 'negative', 'coil', 'exposed'],
        'E': ['polar', 'medium v', 'high p', 'medium pz', 'negative', 'helix', 'exposed'],
        'F': ['apolar', 'large v', 'low p', 'high pz', 'no charge', 'strand', 'buried'],
        'G': ['neutral', 'small v', 'medium p', 'low pz', 'no charge', 'coil', 'buried'],
        'H': ['neutral', 'large v', 'high p', 'high pz', 'no charge', 'helix', 'intermediate'],
        'I': ['apolar', 'medium v', 'low p', 'medium pz', 'no charge', 'strand', 'buried'],
        'K': ['polar', 'large v', 'high p', 'high pz', 'positive', 'helix', 'exposed'],
        'L': ['apolar', 'medium v', 'low p', 'medium pz', 'no charge', 'helix', 'buried'],
        'M': ['apolar', 'large v', 'low p', 'high pz', 'no charge', 'helix', 'intermediate'],
        'N': ['polar', 'medium v', 'high p', 'medium pz', 'no charge', 'coil', 'exposed'],
        'P': ['neutral', 'small v', 'medium p', 'medium pz', 'no charge', 'coil', 'intermediate'],
        'Q': ['polar', 'medium v', 'high p', 'medium pz', 'no charge', 'helix', 'exposed'],
        'R': ['polar', 'large v', 'high p', 'high pz', 'positive', 'helix', 'exposed'],
        'S': ['neutral', 'small v', 'medium p', 'low pz', 'no charge', 'coil', 'intermediate'],
        'T': ['neutral', 'small v', 'medium p', 'low pz', 'no charge', 'strand', 'intermediate'],
        'V': ['apolar', 'medium v', 'low p', 'medium pz', 'no charge', 'strand', 'buried'],
        'W': ['apolar', 'large v', 'low p', 'high pz', 'no charge', 'strand', 'buried'],
        'Y': ['neutral', 'large v', 'low p', 'high pz', 'no charge', 'strand', 'intermediate']
    }

    column_names = [f"G3_{i}" for i in range(42)]
    
    data = {}
    
    output_file = os.path.join(outname,output_file)

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)
        
        pc_properties = [0] * 21  # Physicochemical properties
        t_properties = [0] * 21  # Transition properties
        last = []

        for i, aa in enumerate(seq):
            now = amino_acid_map.get(aa, ['neutral', 'small v', 'medium p', 'low pz', 'no charge', 'helix', 'buried'])
            
            for prop in now:
                pc_properties[p_to_index[prop]] += 1
            
            if last:
                for j in range(7):
                    key1, key2 = last[j], now[j]
                    if key1 != key2:
                        t_properties[p_to_index[key1]] += 1

            last = now

        size_prot = len(seq)
        norm_pc_properties = [x / size_prot for x in pc_properties]
        norm_t_properties = [x / (size_prot - 1) for x in t_properties] if size_prot > 1 else [0] * 21

        data[seq_id] = norm_pc_properties + norm_t_properties  # Store in dictionary
     
        

    features_df = pd.DataFrame.from_dict(data, orient='index', columns=column_names)

    # Save to CSV with the index as protein_name
    features_df.to_csv(output_file, index=True)

    print(f"Processing complete. Results saved in {output_file}")


# Conjoint Triad Composition (CTriad)

def calculate_conjoint_triad_composition(sequence):
    """Computes conjoint triad composition from an amino acid sequence."""
    sequence = list(sequence)
    
    # Handle ambiguous amino acids
    for i in range(len(sequence)):
        if sequence[i] in {'U', 'X'}:
            sequence[i] = 'C'
        elif sequence[i] == 'B':
            sequence[i] = 'D'
        elif sequence[i] == 'Z':
            sequence[i] = 'Q'
        elif sequence[i] == 'J':
            sequence[i] = 'L'
    
    sequence = "".join(sequence)
    
    # Initialize the triad matrix (7x7x7)
    ctriad = np.zeros((7,7,7))
    
    # Amino acid group mapping based on Table 4 (S. Mishra et al.)
    aas = {'A':0,'G':0,'I':1,'L':1,'P':1,'V':0,'F':1,'W':3,'Y':2,'D':5,'E':5,
           'R':4,'H':3,'K':4,'S':2,'T':2,'C':6,'M':2,'N':3,'Q':3}
    
    for i in range(len(sequence) - 2):
        ctriad[aas.get(sequence[i]), aas.get(sequence[i+1]), aas.get(sequence[i+2])] += 1
    
    ctriad /= (len(sequence) - 2) if len(sequence) > 2 else 1  # Normalize
    return ctriad.flatten() 


def extract_ctriad_features(fasta_file,output_file="features_ctriad.csv"):

    """
    Extracts Conjoint Triad Composition features from sequences in a FASTA file and saves them to a CSV file.
    The function processes each sequence in the provided FASTA file, computes its Conjoint Triad Composition 
    features, and writes the results to a CSV file. Each sequence is represented by 343 features, labeled 
    as G4_0 to G4_342.
    Args:
        fasta_file (str): Path to the input FASTA file containing sequences.
        output_file (str, optional): Name of the output CSV file to save the extracted features. 
                                        Defaults to "features_ctriad.csv".
    Returns:
        None
    Raises:
        FileNotFoundError: If the specified FASTA file does not exist.
        ValueError: If the sequence data in the FASTA file is invalid or cannot be processed.
    Notes:
        - The output CSV file will contain one row per sequence, with the sequence ID as the index 
            and the computed features as columns.
        - The function assumes that the `calculate_conjoint_triad_composition` function is defined 
            elsewhere and is responsible for computing the features for a given sequence.
        - The `SeqIO` module from Biopython is used to parse the FASTA file.
    Example:
        >>> extract_ctriad_features("input_sequences.fasta", "output_features.csv")
        Conjoint Triad features saved to: output_features.csv
    """

    # Column names: G4_0 to G4_342 (343 features)
    column_names = [f"G4_{i}" for i in range(343)]
    
    data = {}
    
    output_file = os.path.join(outname,output_file)

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)

        ctriad_features = calculate_conjoint_triad_composition(seq)
        data[seq_id] = ctriad_features
        
    features_df = pd.DataFrame.from_dict(data, orient='index', columns=column_names)
    features_df.to_csv(output_file, index=True)
    print(f"Conjoint Triad features saved to: {output_file}")

# Entropy features

def estimate_shannon_entropy(sequence):
    bases = collections.Counter([tmp_base for tmp_base in sequence])
    dist = [x/sum(bases.values()) for x in bases.values()]
    entropy_value = entropy(dist)
    return entropy_value

def calculate_entropies(sequence, segments=10):
    sequence = list(sequence)
    for i in range(0, len(sequence)):
        if(sequence[i] == 'U' or sequence[i] == 'X'):
            sequence[i] = 'C'
        if(sequence[i] == 'B'):
            sequence[i] = 'D'
        if(sequence[i] == 'Z'):
            sequence[i] = 'Q'
        if(sequence[i] == 'J'):
            sequence[i] = 'L'
    sequence = "".join(sequence)
    length = len(sequence)
    window_size = int(length / segments)
    pos = 0
    scores = [0 for _ in range(segments)]
    for i in range(segments-1):
        scores[i] = estimate_shannon_entropy(sequence[pos:pos + window_size])
        pos += window_size
    scores[-1] = estimate_shannon_entropy(sequence[pos:])
    return scores



def extract_entropy_features(fasta_file,output_file="features_entropy.csv"):

    """
    Extract entropy features from sequences in a FASTA file and save them to a CSV file.
    Args:
        fasta_file (str): Path to the input FASTA file containing sequences.
        output_file (str, optional): Name of the output CSV file to save the extracted features. 
                                        Defaults to "features_entropy.csv".
    Returns:
        None
    Outputs:
        A CSV file containing entropy features for each sequence in the FASTA file. 
        The file is saved to the specified output path.
    Notes:
        - The function calculates entropy features for each sequence in the input FASTA file.
        - The output CSV file contains rows indexed by sequence IDs and columns named "G6_0" to "G6_9".
        - The `calculate_entropies` function is assumed to compute entropy values for the input sequence.
    """



    column_names =["G6_%d" % aa for aa in range(10)]
    
    data = {}
    
    output_file = os.path.join(outname,output_file)
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)
        
        entropies = calculate_entropies(seq)
        data[seq_id] = entropies
        
    features_df = pd.DataFrame.from_dict(data, orient='index', columns=column_names)
    features_df.to_csv(output_file, index=True)
    print(f"Entropy features saved to: {output_file}")


# Protparams features

def calculate_protparams(sequence):
    vec = []
    # Need to get rid of weird aminoacids before
    sequence = list(sequence)
    for i in range(0, len(sequence)):
        if(sequence[i] == 'U' or sequence[i] == 'X'):
            sequence[i] = 'C'
        if(sequence[i] == 'B'):
            sequence[i] = 'D'
        if(sequence[i] == 'Z'):
            sequence[i] = 'Q'
        if(sequence[i] == 'J'):
            sequence[i] = 'L'
    sequence = "".join(sequence)
    X = ProteinAnalysis(sequence)
    vec += [X.molecular_weight()]
    vec += [X.aromaticity()]
    vec += [X.instability_index()]

    # huge vector with size depending on the length of the protein
    # vec += [X.flexibility()]

    vec += [X.gravy()]
    vec += [X.isoelectric_point()]
    molar = X.molar_extinction_coefficient()
    vec += [molar[0]]  # with reduced cysteines
    vec += [molar[1]]  # with disulfid bonds
    
    return vec

def extract_protparams_features(fasta_file,output_file="features_protparams.csv"):

    """
    Extracts protein parameters features from a given FASTA file and saves them to a CSV file.
    This function processes sequences from a FASTA file, calculates protein parameters
    for each sequence, and writes the results to a CSV file. The output file contains
    features for each sequence, with sequence IDs as row indices.
    Args:
        fasta_file (str): Path to the input FASTA file containing protein sequences.
        output_file (str, optional): Name of the output CSV file to save the features.
            Defaults to "features_protparams.csv".
    Returns:
        None
    Raises:
        FileNotFoundError: If the specified FASTA file does not exist.
        ValueError: If there is an issue with the sequence data or feature calculation.
    Notes:
        - The function assumes that the `calculate_protparams` function is defined elsewhere
          and is responsible for computing the protein parameters for a given sequence.
        - The `SeqIO` module from Biopython is used to parse the FASTA file.
        - The output CSV file will be saved in the directory specified by `outname`.
    Example:
        >>> extract_protparams_features("input.fasta", "output_features.csv")
        Protparams features saved to: /path/to/output_features.csv
    """
   
    column_names =["G7_%s" % aa for aa in range(7)]
    
    output_file = os.path.join(outname,output_file)
    
    data = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)
        
        protparams = calculate_protparams(seq)
        data[seq_id] = protparams
    
    features_df = pd.DataFrame.from_dict(data, orient='index', columns=column_names)
    features_df.to_csv(output_file, index=True)
    print(f"Protparams features saved to: {output_file}")
    


# QSO features


AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M",
            "F", "P", "S", "T", "W", "Y", "V"]

# Distance is the Schneider-Wrede physicochemical distance matrix used by Chou
# et. al.

_Distance1 = {"GW": 0.923, "GV": 0.464, "GT": 0.272, "GS": 0.158, "GR": 1.0,
              "GQ": 0.467, "GP": 0.323, "GY": 0.728, "GG": 0.0, "GF": 0.727,
              "GE": 0.807, "GD": 0.776, "GC": 0.312, "GA": 0.206, "GN": 0.381,
              "GM": 0.557, "GL": 0.591, "GK": 0.894, "GI": 0.592, "GH": 0.769,
              "ME": 0.879, "MD": 0.932, "MG": 0.569, "MF": 0.182, "MA": 0.383,
              "MC": 0.276, "MM": 0.0, "ML": 0.062, "MN": 0.447, "MI": 0.058,
              "MH": 0.648, "MK": 0.884, "MT": 0.358, "MW": 0.391, "MV": 0.12,
              "MQ": 0.372, "MP": 0.285, "MS": 0.417, "MR": 1.0, "MY": 0.255,
              "FP": 0.42, "FQ": 0.459, "FR": 1.0, "FS": 0.548, "FT": 0.499,
              "FV": 0.252, "FW": 0.207, "FY": 0.179, "FA": 0.508, "FC": 0.405,
              "FD": 0.977, "FE": 0.918, "FF": 0.0, "FG": 0.69, "FH": 0.663,
              "FI": 0.128, "FK": 0.903, "FL": 0.131, "FM": 0.169, "FN": 0.541,
              "SY": 0.615, "SS": 0.0, "SR": 1.0, "SQ": 0.358, "SP": 0.181,
              "SW": 0.827, "SV": 0.342, "ST": 0.174, "SK": 0.883, "SI": 0.478,
              "SH": 0.718, "SN": 0.289, "SM": 0.44, "SL": 0.474, "SC": 0.185,
              "SA": 0.1, "SG": 0.17, "SF": 0.622, "SE": 0.812, "SD": 0.801,
              "YI": 0.23, "YH": 0.678, "YK": 0.904, "YM": 0.268, "YL": 0.219,
              "YN": 0.512, "YA": 0.587, "YC": 0.478, "YE": 0.932, "YD": 1.0,
              "YG": 0.782, "YF": 0.202, "YY": 0.0, "YQ": 0.404, "YP": 0.444,
              "YS": 0.612, "YR": 0.995, "YT": 0.557, "YW": 0.244, "YV": 0.328,
              "LF": 0.139, "LG": 0.596, "LD": 0.944, "LE": 0.892, "LC": 0.296,
              "LA": 0.405, "LN": 0.452, "LL": 0.0, "LM": 0.062, "LK": 0.893,
              "LH": 0.653, "LI": 0.013, "LV": 0.133, "LW": 0.341, "LT": 0.397,
              "LR": 1.0, "LS": 0.443, "LP": 0.309, "LQ": 0.376, "LY": 0.205,
              "RT": 0.808, "RV": 0.914, "RW": 1.0, "RP": 0.796, "RQ": 0.668,
              "RR": 0.0, "RS": 0.86, "RY": 0.859, "RD": 0.305, "RE": 0.225,
              "RF": 0.977, "RG": 0.928, "RA": 0.919, "RC": 0.905, "RL": 0.92,
              "RM": 0.908, "RN": 0.69, "RH": 0.498, "RI": 0.929, "RK": 0.141,
              "VH": 0.649, "VI": 0.135, "EM": 0.83, "EL": 0.854, "EN": 0.599,
              "EI": 0.86, "EH": 0.406, "EK": 0.143, "EE": 0.0, "ED": 0.133,
              "EG": 0.779, "EF": 0.932, "EA": 0.79, "EC": 0.788, "VM": 0.12,
              "EY": 0.837, "VN": 0.38, "ET": 0.682, "EW": 1.0, "EV": 0.824,
              "EQ": 0.598, "EP": 0.688, "ES": 0.726, "ER": 0.234, "VP": 0.212,
              "VQ": 0.339, "VR": 1.0, "VT": 0.305, "VW": 0.472, "KC": 0.871,
              "KA": 0.889, "KG": 0.9, "KF": 0.957, "KE": 0.149, "KD": 0.279,
              "KK": 0.0, "KI": 0.899, "KH": 0.438, "KN": 0.667, "KM": 0.871,
              "KL": 0.892, "KS": 0.825, "KR": 0.154, "KQ": 0.639, "KP": 0.757,
              "KW": 1.0, "KV": 0.882, "KT": 0.759, "KY": 0.848, "DN": 0.56,
              "DL": 0.841, "DM": 0.819, "DK": 0.249, "DH": 0.435, "DI": 0.847,
              "DF": 0.924, "DG": 0.697, "DD": 0.0, "DE": 0.124, "DC": 0.742,
              "DA": 0.729, "DY": 0.836, "DV": 0.797, "DW": 1.0, "DT": 0.649,
              "DR": 0.295, "DS": 0.667, "DP": 0.657, "DQ": 0.584, "QQ": 0.0,
              "QP": 0.272, "QS": 0.461, "QR": 1.0, "QT": 0.389, "QW": 0.831,
              "QV": 0.464, "QY": 0.522, "QA": 0.512, "QC": 0.462, "QE": 0.861,
              "QD": 0.903, "QG": 0.648, "QF": 0.671, "QI": 0.532, "QH": 0.765,
              "QK": 0.881, "QM": 0.505, "QL": 0.518, "QN": 0.181, "WG": 0.829,
              "WF": 0.196, "WE": 0.931, "WD": 1.0, "WC": 0.56, "WA": 0.658,
              "WN": 0.631, "WM": 0.344, "WL": 0.304, "WK": 0.892, "WI": 0.305,
              "WH": 0.678, "WW": 0.0, "WV": 0.418, "WT": 0.638, "WS": 0.689,
              "WR": 0.968, "WQ": 0.538, "WP": 0.555, "WY": 0.204, "PR": 1.0,
              "PS": 0.196, "PP": 0.0, "PQ": 0.228, "PV": 0.244, "PW": 0.72,
              "PT": 0.161, "PY": 0.481, "PC": 0.179, "PA": 0.22, "PF": 0.515,
              "PG": 0.376, "PD": 0.852, "PE": 0.831, "PK": 0.875, "PH": 0.696,
              "PI": 0.363, "PN": 0.231, "PL": 0.357, "PM": 0.326, "CK": 0.887,
              "CI": 0.304, "CH": 0.66, "CN": 0.324, "CM": 0.277, "CL": 0.301,
              "CC": 0.0, "CA": 0.114, "CG": 0.32, "CF": 0.437, "CE": 0.838,
              "CD": 0.847, "CY": 0.457, "CS": 0.176, "CR": 1.0, "CQ": 0.341,
              "CP": 0.157, "CW": 0.639, "CV": 0.167, "CT": 0.233, "IY": 0.213,
              "VA": 0.275, "VC": 0.165, "VD": 0.9, "VE": 0.867, "VF": 0.269,
              "VG": 0.471, "IQ": 0.383, "IP": 0.311, "IS": 0.443, "IR": 1.0,
              "VL": 0.134, "IT": 0.396, "IW": 0.339, "IV": 0.133, "II": 0.0,
              "IH": 0.652, "IK": 0.892, "VS": 0.322, "IM": 0.057, "IL": 0.013,
              "VV": 0.0, "IN": 0.457, "IA": 0.403, "VY": 0.31, "IC": 0.296,
              "IE": 0.891, "ID": 0.942, "IG": 0.592, "IF": 0.134, "HY": 0.821,
              "HR": 0.697, "HS": 0.865, "HP": 0.777, "HQ": 0.716, "HV": 0.831,
              "HW": 0.981, "HT": 0.834, "HK": 0.566, "HH": 0.0, "HI": 0.848,
              "HN": 0.754, "HL": 0.842, "HM": 0.825, "HC": 0.836, "HA": 0.896,
              "HF": 0.907, "HG": 1.0, "HD": 0.629, "HE": 0.547, "NH": 0.78,
              "NI": 0.615, "NK": 0.891, "NL": 0.603, "NM": 0.588, "NN": 0.0,
              "NA": 0.424, "NC": 0.425, "ND": 0.838, "NE": 0.835, "NF": 0.766,
              "NG": 0.512, "NY": 0.641, "NP": 0.266, "NQ": 0.175, "NR": 1.0,
              "NS": 0.361, "NT": 0.368, "NV": 0.503, "NW": 0.945, "TY": 0.596,
              "TV": 0.345, "TW": 0.816, "TT": 0.0, "TR": 1.0, "TS": 0.185,
              "TP": 0.159, "TQ": 0.322, "TN": 0.315, "TL": 0.453, "TM": 0.403,
              "TK": 0.866, "TH": 0.737, "TI": 0.455, "TF": 0.604, "TG": 0.312,
              "TD": 0.83, "TE": 0.812, "TC": 0.261, "TA": 0.251, "AA": 0.0,
              "AC": 0.112, "AE": 0.827, "AD": 0.819, "AG": 0.208, "AF": 0.54,
              "AI": 0.407, "AH": 0.696, "AK": 0.891, "AM": 0.379, "AL": 0.406,
              "AN": 0.318, "AQ": 0.372, "AP": 0.191, "AS": 0.094, "AR": 1.0,
              "AT": 0.22, "AW": 0.739, "AV": 0.273, "AY": 0.552, "VK": 0.889}


def GetAAComposition(ProteinSequence):
    LengthSequence = len(ProteinSequence)
    Result = {}
    for i in AALetter:
        Result[i] = round(float(ProteinSequence.count(i)) / LengthSequence, 3)
    return Result



def GetSequenceOrderCouplingNumber(ProteinSequence, d=1,
                                   distancematrix=_Distance1):
    NumProtein = len(ProteinSequence)
    tau = 0.0
    for i in range(NumProtein - d):
        temp1 = ProteinSequence[i]
        temp2 = ProteinSequence[i + d]
        if(temp1 == 'U' or temp1 == 'X'):
            temp1 = 'C'
        if(temp1 == 'B'):
            temp1 = 'D'
        if(temp1 == 'Z'):
            temp1 = 'Q'
        if(temp1 == 'J'):
            temp1 = 'L'
        if(temp2 == 'U' or temp2 == 'X'):  # selenocysteine
            temp2 = 'C'
        if(temp2 == 'B'):
            temp2 = 'D'
        if(temp2 == 'Z'):
            temp2 = 'Q'
        if(temp2 == 'J'):
            temp2 = 'L'
    return round(tau, 3)

def GetSequenceOrderCouplingNumberSW(ProteinSequence, maxlag=30,
                                     distancematrix=_Distance1):

    Tau = {}
    for i in range(maxlag):
        Tau["tausw" + str(i + 1)] = GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix)
    return Tau



def GetQuasiSequenceOrder1SW(ProteinSequence, maxlag=30,
                             weight=0.1, distancematrix=_Distance1):

    rightpart = 0.0
    for i in range(maxlag):
        rightpart = rightpart + GetSequenceOrderCouplingNumber(
            ProteinSequence, i + 1, distancematrix)
    AAC = GetAAComposition(ProteinSequence)
    result = {}
    temp = 1 + weight * rightpart
    for index, i in enumerate(AALetter):
        result['QSOSW' + str(index + 1)] = round(AAC[i] / temp, 6)

    return result


def GetQuasiSequenceOrder2SW(ProteinSequence, maxlag=30, weight=0.1,
                             distancematrix=_Distance1):

    rightpart = []
    for i in range(maxlag):
        rightpart.append(
            GetSequenceOrderCouplingNumber(
                ProteinSequence, i + 1, distancematrix))
    result = {}
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result['QSOSW' + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp, 6)

    return result




def calculate_qso(ProteinSequence, maxlag=30, weight=0.1):
    result = dict()
    result.update(GetQuasiSequenceOrder1SW(ProteinSequence, maxlag,
                                           weight, _Distance1))
    result.update(
        GetQuasiSequenceOrder2SW(ProteinSequence, maxlag, weight, _Distance1))
    return result

def extract_qso_features(fasta_file,output_file="features_qso.csv"):
    """
    Extract QSO features
    """
    
    output_file = os.path.join(outname,output_file)

    data = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)
        
        qso = calculate_qso(seq)
        data[seq_id] = qso

    column_names = ["G8_%s" % i for i in qso.keys()]    
    
    features_df = pd.DataFrame.from_dict(data, orient='index')
    features_df.columns = column_names
    features_df.to_csv(output_file, index=True)
    print(f"QSO features saved to: {output_file}")
    
    
# Word2Vec features

def calculate_3grams2(sequence):
    sequence = list(sequence)
    for i in range(0, len(sequence)):
        if(sequence[i] == 'U' or sequence[i] == 'X'):
            sequence[i] = 'C'
        if(sequence[i] == 'B'):
            sequence[i] = 'D'
        if(sequence[i] == 'Z'):
            sequence[i] = 'Q'
        if(sequence[i] == 'J'):
            sequence[i] = 'L'
    sequence = "".join(sequence)
    first3gram = []
    for i in range(0, len(sequence)-2):
        first = sequence[i] + sequence[i+1] + sequence[i+2]
        first3gram.append(first)
    return first3gram


def calculate_word2vec(sequence, model):
    grams = calculate_3grams2(sequence)
    protein_vec = None
    for i in grams:
        vector = model.wv[str(i)]
        if protein_vec is None:
            protein_vec = np.array([vector])
        else:
            protein_vec += np.array([vector])
    return protein_vec.flatten()



def extract_word2vec(fasta_file, output_file="features_word2vec.csv",model_filename = "/home/kwamouns/Datas/2023-sf-secretions/scripts/Features_Extraction/word2vec_train.txt"):

    """
    Extract Word2Vec features from sequences in a FASTA file and save them to a CSV file.
    Args:
        fasta_file (str): Path to the input FASTA file containing sequences.
        output_file (str, optional): Name of the output CSV file to save the extracted features. 
            Defaults to "features_word2vec.csv".
        model_filename (str, optional): Path to the pre-trained Word2Vec model file. 
            Defaults to "/home/kwamouns/Datas/2023-sf-secretions/scripts/Features_Extraction/word2vec_train.txt".
    Returns:
        None: The function saves the extracted features to the specified CSV file.
    Notes:
        - Each sequence in the FASTA file is processed to calculate its Word2Vec representation.
        - The Word2Vec features are stored in a pandas DataFrame, with one row per sequence.
        - The output CSV file contains the sequence IDs as the index and the Word2Vec features as columns.
        - The column names are generated as "G5_0", "G5_1", ..., based on the length of the Word2Vec vector.
    Example:
        extract_word2vec("input_sequences.fasta", "output_features.csv", "word2vec_model.txt")
    """


    # results/features/2023-08-22_NCBI-genomes_standardized/word2vec_train.txt ( in morrowind)
    
    output_file = os.path.join(outname,output_file)
    
    data = {}
    model = Word2Vec.load(model_filename)
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)
        
        word2vec = calculate_word2vec(seq, model)
        data[seq_id] = word2vec
    
    column_names = ["G5_%s" % i for i in range(len(word2vec))]
    
    features_df = pd.DataFrame.from_dict(data, orient='index')
    features_df.columns = column_names
    features_df.to_csv(output_file, index=True)
    print(f"Word2Vec features saved to: {output_file}")
    
    

# Only run this code if the module is executed directly, not when imported
if __name__ == '__main__':
    global outname

    #outname = '/home/kwamouns/Datas/2023-sf-secretions/results/script_folder_results/ESKAPE/GCF_000006765.1_ASM676v1/features'

    # Diderm with No SST
    outname = '/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/t7ss_features'
    os.makedirs(outname, exist_ok=True)

    extract_aac_features("/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/uniprotkb_t7ss_AND_reviewed_true_2026_04_14.fasta")
    extract_dcp_features("/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/uniprotkb_t7ss_AND_reviewed_true_2026_04_14.fasta")
    extract_pcp_features("/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/uniprotkb_t7ss_AND_reviewed_true_2026_04_14.fasta")
    extract_ctriad_features("/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/uniprotkb_t7ss_AND_reviewed_true_2026_04_14.fasta")
    extract_entropy_features("/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/uniprotkb_t7ss_AND_reviewed_true_2026_04_14.fasta")
    extract_protparams_features("/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/uniprotkb_t7ss_AND_reviewed_true_2026_04_14.fasta")
    extract_qso_features("/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/uniprotkb_t7ss_AND_reviewed_true_2026_04_14.fasta")
    extract_word2vec("/home/kwamouns/Datas/secret_project_2026/2023-sf-secretions/secret_pipeline1/External_val_data/uniprotkb_t7ss_AND_reviewed_true_2026_04_14.fasta")

