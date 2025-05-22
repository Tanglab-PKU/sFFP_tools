"""
Cross-link Distance Analysis Pipeline

This script calculates spatial distances between cross-linked peptide residues 
using structural data from PDB files. It processes CSV files containing peptide
cross-link information and matches them against protein structures.

Key Features:
- PDB structure parsing with Biopython
- Dynamic chain pair optimization
- Minimum distance calculation between residues
- Multi-position handling (e.g., ambiguous cross-link positions)
- Sequence mapping validation
"""

import pandas as pd
from Bio.PDB import PDBParser, PPBuilder
import os
from itertools import product
from pathlib import Path
import itertools

# Configuration Constants
CONFIG = {
    "csv_root": r'D:\MSdata\250306-BD3\HS90B-BD3(DESTHY)-5C_1FDR',
    "pdb_dir": r"D:\MSdata\pdb\HS90B"
}

def calculate_distance(residue1, residue2, atom1='CA', atom2='CA'):
    """Calculate Euclidean distance between specified atoms in two residues
    
    Args:
        residue1, residue2 (Bio.PDB.Residue): Residue objects
        atom1, atom2 (str): Atom names to compare (default: 'CA')
        
    Returns:
        float: Distance in Ångströms or None if atoms missing
    """
    try:
        return abs(residue1[atom1] - residue2[atom2])
    except KeyError as e:
        print(f"Missing atom: {str(e)}")
        return None

def parse_positions(position_str):
    """Parse position string into list of integers
    
    Handles:
        - Single positions (e.g., '123')
        - Ambiguous positions (e.g., '122/123')
        - Invalid entries
    
    Args:
        position_str (str): Position string from CSV
        
    Returns:
        list[int]: Cleaned positions or [None] for invalid entries
    """
    try:
        if "/" in str(position_str):
            return [int(p) for p in str(position_str).split('/')]
        return [int(position_str)]
    except ValueError:
        return [None]

def get_chain_sequences(model):
    """Extract amino acid sequences from PDB chains
    
    Args:
        model (Bio.PDB.Model): First model in PDB structure
        
    Returns:
        dict: {chain_id: sequence} with I→L substitution
    """
    ppb = PPBuilder()
    chain_sequences = {}
    for chain in model:
        # Build peptides and concatenate sequences
        peptides = ppb.build_peptides(chain)
        seq = ''.join([str(p.get_sequence()) for p in peptides])
        chain_sequences[chain.id] = str(seq).replace('I', 'L')
    return chain_sequences

# File discovery --------------------------------------------------------------
csv_files = []
pdb_files = []

# Find all cross-link CSV files
for root, dirs, _ in os.walk(CONFIG["csv_root"]):
    if "reports" in dirs and "tmp" not in root:
        reports_dir = os.path.join(root, "reports")
        csv_files.extend([
            os.path.join(reports_dir, f) 
            for f in os.listdir(reports_dir) 
            if f.endswith("_XLs_Jnoted.csv")
        ])

# Find all PDB files
pdb_files = [
    os.path.join(CONFIG["pdb_dir"], f) 
    for f in os.listdir(CONFIG["pdb_dir"]) 
    if f.endswith(".pdb")
]

# Main processing loop --------------------------------------------------------
for csv_file in csv_files:
    # Load cross-link data with size limit
    base_df = pd.read_csv(csv_file, low_memory=False)[:50000]
    results = base_df.copy()
    
    # Pre-cache unique peptide sequences for fast lookup
    pep1_seqs = set(base_df['Pep1_seq'].dropna().unique())
    pep2_seqs = set(base_df['Pep2_seq'].dropna().unique())

    # Process each PDB structure
    for pdb_path in pdb_files:
        pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
        structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
        model = structure[0]  # Use first model only
        
        # Map chain IDs to sequences
        chain_sequences = get_chain_sequences(model)
        
        # Optimized chain pairing (critical performance improvement)
        valid_chains1 = [
            cid for cid, seq in chain_sequences.items()
            if any(pep in seq for pep in pep1_seqs)
        ]
        valid_chains2 = [
            cid for cid, seq in chain_sequences.items()
            if any(pep in seq for pep in pep2_seqs)
        ]
        chain_pairs = list(itertools.product(valid_chains1, valid_chains2))

        # Calculate distances for each chain pair
        for chain_a, chain_b in chain_pairs:
            col_name = f"{pdb_id}_{chain_a}-{chain_b}_Dis"
            distance_data = []
            
            for _, row in base_df.iterrows():
                positions1 = parse_positions(row.get("Position1", 0))
                positions2 = parse_positions(row.get("Position2", 0))
                min_dist = None  # Track minimum valid distance
                
                # Sequence presence validation
                if (row["Pep1_seq"] in chain_sequences[chain_a] and 
                    row["Pep2_seq"] in chain_sequences[chain_b]):
                    
                    # Test all position combinations
                    for pos1, pos2 in itertools.product(positions1, positions2):
                        try:
                            res1 = model[chain_a][pos1]
                            res2 = model[chain_b][pos2]
                            dist = calculate_distance(res1, res2)
                            
                            # Update minimum distance
                            if dist is not None:
                                if min_dist is None or dist < min_dist:
                                    min_dist = dist
                        except (KeyError, AttributeError):
                            continue  # Silent skip for missing residues
                
                distance_data.append(min_dist)
            
            results[col_name] = distance_data
    
    # Post-processing ---------------------------------------------------------
    # Merge original data with distance calculations
    merged_df = pd.concat([base_df, results.filter(regex='_Dis$')], axis=1)
    
    # Calculate global minimum distance
    distance_cols = merged_df.filter(regex='_Dis$').columns
    merged_df['Min_Distance'] = merged_df[distance_cols].apply(
        lambda row: min([x for x in row if pd.notnull(x) and x > 0], default=None),
        axis=1
    )
    
    # Save results
    csv_path = Path(csv_file)
    output_path = csv_path.with_name(f"{csv_path.stem}_Dis.csv")
    merged_df.to_csv(output_path, index=False)
    print(f"Distance file generated: {output_path}")

print("Batch processing completed")