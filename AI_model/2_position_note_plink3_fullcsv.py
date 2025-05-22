"""
Cross-linked Peptide Data Processing Pipeline

This script processes mass spectrometry reports containing cross-linked peptides,
extracts key information, and generates standardized output files.

Key Features:
- Recursive directory traversal for batch processing
- Regular expression-based protein ID parsing
- Cross-link key normalization to handle position swaps
- Peptide sequence validation and amino acid residue extraction
- PSM (Peptide Spectrum Match) counting
- Score-based result sorting
"""

import pandas as pd
from pathlib import Path
import os
import re

def process_directory(root_dir):
    """
    Process all qualifying CSV files in directory tree
    
    Args:
        root_dir (str): Root directory containing reports subdirectories
        
    Workflow:
        1. Recursively search for CSV files in 'reports' subdirectories
        2. Process files with numeric suffix in filename (indicating raw data)
        3. Generate standardized output files
    """
    for path in Path(root_dir).rglob('reports/*.csv'):
        # Identify raw data files (with numeric suffix)
        if path.stem[-1].isdigit():
            output_path = process_file(path)
            print("Output generated:", output_path)

def parse_crosslink(sequence):
    """
    Parse cross-linked peptide sequence into standardized format
    
    Args:
        sequence (str): Cross-link sequence in format 'PEPTIDE(pos)-PEPTIDE(pos)/'
        
    Returns:
        tuple: Sorted ((sequence1, position1), (sequence2, position2))
        
    Raises:
        ValueError: For invalid input formats
    """
    # Handle trailing slash and split peptides
    pep1, pep2 = sequence.rstrip('/').split('-')  
    
    # Regular expression matching for sequence components
    pattern = r'(\w+)\((\d+)\)'
    pep1_match = re.match(pattern, pep1)
    pep2_match = re.match(pattern, pep2)
    
    if not pep1_match or not pep2_match:
        raise ValueError("Invalid format. Expected 'SEQUENCE(pos)-SEQUENCE(pos)'")
        
    # Extract sequence components
    pep1_seq, pep1_site = pep1_match.groups()
    pep2_seq, pep2_site = pep2_match.groups()
    
    # Return sorted tuple for consistent ordering
    return tuple(sorted([(pep1_seq, int(pep1_site)), (pep2_seq, int(pep2_site))]))

def process_file(input_path):
    """
    Process individual CSV file containing cross-link data
    
    Args:
        input_path (Path): Path object for input file
        
    Returns:
        Path: Output file path
        
    Processing Steps:
        1. Load CSV with UTF-8 encoding
        2. Filter cross-link entries (Peptide_Type == 3)
        3. Extract protein IDs using regex
        4. Generate normalized cross-link keys
        5. Validate peptide positions
        6. Calculate PSMs per cross-link
        7. Save sorted results
    """
    # Regex pattern for UniProt-style protein IDs
    protein_pattern = r'.*?(sp|lcl)\|([^|]+)\|([^\s]+) \((\d+)\)-(sp|lcl)\|([^|]+)\|([^\s]+) \((\d+).*'
    
    # Load data with RFC4180 CSV standard compliance
    df_original = pd.read_csv(input_path, encoding='utf-8', low_memory=False)
    
    # Filter cross-links (Type 3) and limit to 50k entries
    df = df_original[df_original['Peptide_Type'] == 3].copy()[:50000]
    
    # Extract protein components using regex groups
    extracted = df['Proteins'].str.extract(protein_pattern)
    extracted.columns = [
        'Prefix1', 'ID1', 'Protein1', 'Position1',
        'Prefix2', 'ID2', 'Protein2', 'Position2'
    ]
    
    # Merge extracted data with original dataframe
    result_df = pd.concat([
        df.reset_index(drop=True),
        extracted.reset_index(drop=True)
    ], axis=1)
    
    # Generate normalized cross-link keys (prevents position swapping issues)
    result_df['_key'] = result_df.apply(
        lambda x: tuple(sorted([
            (x['Protein1'], x['Position1']),
            (x['Protein2'], x['Position2'])
        ])),
        axis=1
    )
    
    # Parse peptide sequences and validate ordering
    result_df['_parsed'] = result_df['Peptide'].apply(parse_crosslink)
    result_df[['Pep1', 'Pep2']] = result_df.apply(
        lambda x: x['_parsed'] if x['_key'] == x['_parsed'] else x['_parsed'][::-1],
        axis=1, result_type='expand'
    )
    
    # Extract amino acid residues with position validation
    aa_extraction = lambda pep: (
        pep[0][pep[1]-1] if pep[1] <= len(pep[0]) else 'X'
    )
    result_df['Pep1_AA'] = result_df['Pep1'].apply(aa_extraction)
    result_df['Pep2_AA'] = result_df['Pep2'].apply(aa_extraction)
    
    # Store full peptide sequences
    result_df['Pep1_seq'] = result_df['Pep1'].apply(lambda x: x[0])
    result_df['Pep2_seq'] = result_df['Pep2'].apply(lambda x: x[0])
    
    # Clean intermediate columns
    result_df.drop(columns=['_parsed'], inplace=True)
    
    # Calculate PSMs (Peptide Spectrum Matches) per cross-link
    psm_counts = result_df.groupby('_key').size().reset_index(name='PSMs')
    result_df = result_df.merge(psm_counts, on='_key').drop(columns='_key')
    
    # Save sorted results (descending by Score)
    output_path = input_path.with_name(f"{input_path.stem}_XLs_Jnoted.csv")
    result_df.sort_values(by='Score', ascending=False, inplace=True)
    result_df.to_csv(output_path, index=False)
    
    return output_path


# Execution example (replace with actual data path)
if __name__ == "__main__":
    DEFAULT_PATH = r'D:\MSdata\250306-BD3\HS90B-BD3(DESTHY)-5C_1FDR'
    process_directory(DEFAULT_PATH)
    print("Batch processing completed")