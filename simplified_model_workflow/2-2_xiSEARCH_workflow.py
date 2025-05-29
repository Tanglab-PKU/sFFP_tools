"""
Data Processing Script for Cross-linked Peptide Analysis

This script processes mass spectrometry data to filter and annotate cross-linked peptides
for subsequent fragmentation analysis. It calculates peptide-spectrum matches (PSMs) and
generates ion intensity metrics for structural validation.

Processing Steps:
1. Filter peptides by type and confidence score
2. Calculate peptide-level PSMs
3. Generate cross-linker specific ion intensity columns
"""

import pandas as pd
import re
import bisect
from collections import defaultdict
from numpy import abs
import pandas as pd
import os
from pathlib import Path
import numpy as np

# Read input data
current_dir = Path(__file__).parent
input_file = current_dir / 'xiSEARCH_BSA-SDA_Xi1.8.7.csv'
input_df = pd.read_csv(input_file)
output_xlsx = current_dir / "xiSEARCH_output.xlsx"
MGF_CSV_PATH = current_dir / "test_data" / "mgf_to_csv"
# Create output directory for spectra
SPECTRUM_OUTPUT = current_dir / "xiSEARCH_Jsearch_spectrums"
SPECTRUM_OUTPUT.mkdir(exist_ok=True)


# Step 1: Filter peptides
# Make sure 'Modifications1' and 'Modifications2' are treated as strings to handle possible NaN values
input_df['Modifications1'] = input_df['Modifications1'].astype(str)
input_df['Modifications2'] = input_df['Modifications2'].astype(str)
filter_condition = (input_df['Crosslinker'] == 'SDA') & (~input_df['Modifications1'].str.contains("->")) & (~input_df['Modifications2'].str.contains("->"))
processed_df = input_df.loc[filter_condition].copy().sort_values(by='NormScore', ascending=False)
processed_df = processed_df[:10000]

# Step 2: Calculate peptide PSMs
processed_df['Peptide_PSMs'] = processed_df.groupby(['BasePeptide1', 'BasePeptide2',"Link1","Link2"])['BasePeptide1'].transform('size')

# Step 3: Initialize cross-linker specific columns

# Define atomic masses of common elements (monoisotopic)
ELEMENT_MASSES = {
    'H': 1.00783,    # Hydrogen
    'C': 12.0000,    # Carbon
    'N': 14.0031,    # Nitrogen
    'O': 15.9949,    # Oxygen
    'S': 31.9721     # Sulfur
}

# Define molecular formulas for each standard amino acid
AMINO_ACID_FORMULAS = {
    'A': {'C': 3, 'H': 7, 'N': 1, 'O': 2},      # Alanine
    'R': {'C': 6, 'H': 14, 'N': 4, 'O': 2},     # Arginine
    'N': {'C': 4, 'H': 8, 'N': 2, 'O': 3},      # Asparagine
    'D': {'C': 4, 'H': 7, 'N': 1, 'O': 4},      # Aspartic acid
    'C': {'C': 5, 'H': 10, 'N': 2, 'O': 3, 'S': 1},  # Cysteine # Default Carbamidomethyl[C]: 57.021464 H(3)C(2)N(1)O(1)
    'E': {'C': 5, 'H': 9, 'N': 1, 'O': 4},      # Glutamic acid
    'Q': {'C': 5, 'H': 10, 'N': 2, 'O': 3},     # Glutamine
    'G': {'C': 2, 'H': 5, 'N': 1, 'O': 2},      # Glycine
    'H': {'C': 6, 'H': 9, 'N': 3, 'O': 2},      # Histidine
    'I': {'C': 6, 'H': 13, 'N': 1, 'O': 2},     # Isoleucine
    'L': {'C': 6, 'H': 13, 'N': 1, 'O': 2},     # Leucine
    'K': {'C': 6, 'H': 14, 'N': 2, 'O': 2},     # Lysine
    'M': {'C': 5, 'H': 11, 'N': 1, 'O': 2, 'S': 1},  # Methionine
    'F': {'C': 9, 'H': 11, 'N': 1, 'O': 2},     # Phenylalanine
    'P': {'C': 5, 'H': 9, 'N': 1, 'O': 2},      # Proline
    'S': {'C': 3, 'H': 7, 'N': 1, 'O': 3},      # Serine
    'T': {'C': 4, 'H': 9, 'N': 1, 'O': 3},      # Threonine
    'W': {'C': 11, 'H': 12, 'N': 2, 'O': 2},    # Tryptophan
    'Y': {'C': 9, 'H': 11, 'N': 1, 'O': 3},     # Tyrosine
    'V': {'C': 5, 'H': 11, 'N': 1, 'O': 2}      # Valine
}

# Mass of water (used in peptide bond calculations)
H2O_MASS = 2 * ELEMENT_MASSES['H'] + ELEMENT_MASSES['O']

# Cross-linker specific modifications (BD3 cleavable cross-linker)
XL_MODIFICATION_MASSES = {        
    'sb': -18.011,      # Side-chain fragment b
    'sc': 0.0,          # Side-chain fragment c
    "sz": 82.042,      # Full cross-linker mass
    "sy": 100.052,
}

# Common peptide modifications
MODIFICATION_MASSES = { 
    # "Carbamidomethyl[C]": 57.021464,    # Cysteine alkylation
    "Mox": 15.994915,          # Methionine oxidation
    "ac": 42.010565,     # N-terminal acetylation
}

def calculate_peptide_mass(sequence, modifications=None):
    """
    Calculate the monoisotopic mass of a peptide with optional modifications
    
    Args:
        sequence (str): Amino acid sequence
        modifications (dict): Dictionary of modifications by position
        
    Returns:
        float: Calculated mass rounded to 4 decimal places
    """
    total_mass = 0.0
    
    if len(sequence) > 0:
        # Calculate base mass from amino acid composition
        for aa in sequence:
            for element, count in AMINO_ACID_FORMULAS[aa].items():
                total_mass += count * ELEMENT_MASSES[element]
        
        # Add modification masses if present
        if modifications:
            for pos in modifications:
                if pos <= len(sequence):
                    for mod in modifications[pos]:
                        if "ac" in mod:
                            mod="ac"
                        total_mass += MODIFICATION_MASSES.get(mod, 0)

        # Subtract water molecules lost during peptide bond formation
        total_mass -= (len(sequence) - 1) * H2O_MASS
        return round(total_mass, 4)
    else:
        return round(H2O_MASS, 4)  # Return just water mass for empty sequence
    


def parse_modifications(mod_str_A, mod_pos_A, mod_str_B, mod_pos_B):
    """
    Parse modification string and map modifications to the correct peptide.
    """
    mods = {'A': defaultdict(list), 'B': defaultdict(list)}
    # print(mod_str_A, mod_pos_A, mod_str_B, mod_pos_B)
    # Handle missing/empty modification strings
    if (pd.isna(mod_str_A) or mod_str_A == 'nan') and (pd.isna(mod_str_B) or mod_str_B == 'nan'):
        return mods

    else:
        def parse_mod_string(mod_str, mod_pos, peptide_key):
            if not pd.isna(mod_str) and mod_str != 'nan':
                mod_str_list = str(mod_str).split(";") if ";" in mod_str else [str(mod_str)]
                mod_pos_list = str(mod_pos).split(";") if ";" in mod_pos else [str(mod_pos)]

                for pos_str in mod_pos_list:
                    if pos_str and pos_str != 'nan':
                        pos = int(pos_str)
                        # Check if position is valid
                        if pos < len(mod_str_list):
                            mods[peptide_key][pos].append(mod_str_list[pos])

        parse_mod_string(mod_str_A, mod_pos_A, 'A')
        parse_mod_string(mod_str_B, mod_pos_B, 'B')
        
        return mods

def calculate_ion_masses(peptide, crosslink_site, modifications=None, XL_mods=XL_MODIFICATION_MASSES.keys()):
    """
    Calculate theoretical b/y ion masses for a peptide with cross-linker modifications
    
    Args:
        peptide (str): Amino acid sequence
        crosslink_site (int): Position of cross-linking
        modifications (dict): Peptide modifications
        XL_mods (list): List of cross-linker modifications to consider
        
    Returns:
        tuple: Lists of b ions and y ions with their properties
    """
    b_ions, y_ions = [], []
    peptide_len = len(peptide)

    for i in range(0, peptide_len):
        # Calculate mass for fragment
        sub_mass = calculate_peptide_mass(peptide[:i], modifications)

        # Determine if fragment contains crosslink site
        b_has_crosslink = i >= crosslink_site  # b ions are modified if past crosslink site

        # Base masses for b and y ions
        b_base = sub_mass - H2O_MASS
        y_base = calculate_peptide_mass(peptide, modifications) - b_base

        if b_has_crosslink:
            # Calculate b ions with crosslink modifications (charges 1-3)
            for mod in XL_mods:
                
                if mod != "N":
                    for charge in [1, 2, 3]:
                        b_ions.append({
                            "mz": round((b_base + charge * ELEMENT_MASSES['H'] + XL_MODIFICATION_MASSES[mod]) / charge, 4),
                            "position": i,  # Forward position
                            "charge": charge,
                            "mod": mod,
                            "isotope_mz": round((b_base + charge * ELEMENT_MASSES['H'] + XL_MODIFICATION_MASSES[mod] + 1) / charge, 4),
                        }) 
                elif mod == "N":
                    # Calculate y ions without crosslink modification
                    for charge in [1, 2, 3]:
                        y_ions.append({
                            "mz": round((y_base + charge * ELEMENT_MASSES['H']) / charge, 4),
                            "position": peptide_len - i,  # Reverse position
                            "charge": charge,
                            "mod": "N",
                            "isotope_mz": round((y_base + charge * ELEMENT_MASSES['H'] + 1) / charge, 4),
                        })
        else:
            # Handle fragments before crosslink site
            for mod in XL_mods:
                if mod == "N":
                    # Calculate b ions without crosslink modification
                    for charge in [1, 2, 3]:
                        b_ions.append({
                            "mz": round((b_base + charge * ELEMENT_MASSES['H']) / charge, 4),
                            "position": i,
                            "charge": charge,
                            "mod": "N",
                            "isotope_mz": round((b_base + charge * ELEMENT_MASSES['H'] + 1) / charge, 4),
                        }) 
                elif mod != "N":
                    # Calculate y ions with crosslink modifications
                    for charge in [1, 2, 3]:
                        y_ions.append({
                            "mz": round((y_base + charge * ELEMENT_MASSES['H'] + XL_MODIFICATION_MASSES[mod]) / charge, 4),
                            "position": peptide_len - i,
                            "charge": charge,
                            "mod": mod,
                            "isotope_mz": round((y_base + charge * ELEMENT_MASSES['H'] + XL_MODIFICATION_MASSES[mod] + 1) / charge, 4),
                        })
        
    # Special case for y0 ions (full peptide with crosslink modifications)
    for mod in XL_mods:
        for charge in [1, 2, 3]:
            if mod != "N":
                y_ions.append({
                    "mz": round((calculate_peptide_mass(peptide, modifications) + charge * ELEMENT_MASSES['H'] + XL_MODIFICATION_MASSES[mod]) / charge, 4),
                    "position": 0,
                    "charge": charge,
                    "mod": mod,
                    "isotope_mz": round((calculate_peptide_mass(peptide, modifications) + charge * ELEMENT_MASSES['H'] + XL_MODIFICATION_MASSES[mod] + 1) / charge, 4),
                })
    return b_ions, y_ions

def match_spectrum(theoretical, observed, tolerance=0.02, intensity_range=(0.05, 1.95)):
    """
    Match theoretical ions to observed spectrum peaks with isotope validation
    
    Args:
        theoretical (dict): Theoretical ions by type
        observed (DataFrame): Observed peaks (m/z and intensity)
        tolerance (float): Mass tolerance for matching
        intensity_range (tuple): Valid isotope intensity ratio range
        
    Returns:
        dict: Dictionary mapping observed peak indices to matched ion labels
    """
    matches = defaultdict(list)
    
    # Structure theoretical ions for efficient searching
    ion_list = []
    for ion_type, masses in theoretical.items():
        for mz_data in masses:
            ion_list.append({
                "mz": mz_data['mz'],
                "type": ion_type,
                "pos": mz_data['position'],
                "charge": mz_data['charge'],
                "mod": mz_data['mod'],
                "isotope_mz": mz_data["isotope_mz"]
            })

    # Sort ions by m/z for efficient searching
    sorted_ions = sorted(ion_list, key=lambda x: x["mz"])

    # Sort observed data for binary search
    observed_sorted = observed.sort_values(by='m/z').reset_index(drop=True)

    # Binary search for matches within tolerance
    for idx, row in observed_sorted.iterrows():
        mz = row['m/z']
        intensity = row['intensity']
        
        # Find potential matches within tolerance window
        left = bisect.bisect_left([x["mz"] for x in sorted_ions], mz - tolerance)
        right = bisect.bisect_right([x["mz"] for x in sorted_ions], mz + tolerance)

        for ion in sorted_ions[left:right]:
            if abs(ion["mz"] - mz) <= tolerance:
                # Check for isotope peak if available
                if ion['isotope_mz']:
                    iso_mz = ion['isotope_mz']
                    iso_left = bisect.bisect_left(observed_sorted['m/z'], iso_mz - tolerance)
                    iso_right = bisect.bisect_right(observed_sorted['m/z'], iso_mz + tolerance)
                    
                    # Validate isotope peak intensity ratio
                    for iso_idx in range(iso_left, iso_right):
                        iso_row = observed_sorted.iloc[iso_idx]
                        if abs(iso_row['m/z'] - iso_mz) <= tolerance:
                            observed_iso_intensity = iso_row['intensity']
                            if intensity_range[0] <= observed_iso_intensity / intensity <= intensity_range[1]:
                                label = f"{ion['type']}{ion['pos']}{'+' * ion['charge']}"
                                matches[idx].append(f"{label}[{ion['mod']}]")
                                break
                else:
                    # No isotope check needed
                    label = f"{ion['type']}{ion['pos']}{'+' * ion['charge']}"
                    matches[idx].append(f"{label}[{ion['mod']}]")
    return matches
    
def process_row(row, mgf_csv_path, output_folder):
    """
    Process a single row from the input data file
    
    Args:
        row (Series): DataFrame row containing peptide information
        mgf_csv_path (Path): Path to MGF CSV files
        output_folder (Path): Folder to save output files
        
    Returns:
        dict: Dictionary containing analysis results for the row
    """
    
    try:
        # Load corresponding spectrum data
        csv_file = Path(mgf_csv_path) / f"{row['ScanTitle']}.csv"
        df = pd.read_csv(csv_file, low_memory=False)

        coverage_results = {}

        # Parse cross-linked peptide information
        pepA=row["BasePeptide1"]
        siteA=row["Link1"]
        pepB=row["BasePeptide2"]
        siteB=row["Link2"]
        mods = parse_modifications(row["Modifications1"],row["ModificationPositions1"],row["Modifications2"],row["ModificationPositions2"])
        # print(mods)    

        # Calculate masses for cross-linker attached to full peptides
        XL_MODIFICATION_MASSES["pepA"] = calculate_peptide_mass(pepA, modifications=mods['A']) + 82.042
        XL_MODIFICATION_MASSES["pepB"] = calculate_peptide_mass(pepB, modifications=mods['B']) + 82.042


        # Calculate global maximum intensity for normalization
        max_intensity = df['intensity'].max()
        max_intensity = max_intensity if max_intensity != 0 else 1  # Avoid division by zero

        # Define cross-linker modification configurations to test
        XL_MOD_CONFIGS = [
            {'name': 'base_sc', 'A_mods': ['N','pepB','sc', 'sz'], 
                'B_mods': ['N','pepA','sc', 'sz']},  # All fragments
        ]

        for config in XL_MOD_CONFIGS:
            # Calculate theoretical ions for both peptides
            Ab, Ay = calculate_ion_masses(pepA, siteA, mods['A'], XL_mods=config['A_mods'])
            Bb, By = calculate_ion_masses(pepB, siteB, mods['B'], XL_mods=config['B_mods'])

            # Match theoretical ions to observed spectrum
            matches_A = match_spectrum({'Ab': Ab, 'Ay': Ay}, df[['m/z', 'intensity']])
            matches_B = match_spectrum({'Bb': Bb, 'By': By}, df[['m/z', 'intensity']])
            
            
            # Calculate intensity statistics for each peptide chain
            for chain, matches in {'A': matches_A, 'B': matches_B}.items():
                position_intensity = defaultdict(float)
                peptide_length = len(pepA) if chain == 'A' else len(pepB)
                
                # Record maximum intensity for each ion position
                for idx, ion_labels in matches.items():
                    current_intensity = df.at[idx, 'intensity']
                    for label in ion_labels:
                        match = re.match(rf'^{chain}([by])(\d+)\+*\[.*?]', label)
                        if match:
                            ion_type, pos_str = match.groups()
                            pos = int(pos_str) 
                            # Only consider valid fragmentation sites
                            if 1 <= pos < peptide_length:
                                position_intensity[ion_type + pos_str] = max(
                                    position_intensity[ion_type + pos_str], 
                                    current_intensity
                                )

                # Calculate statistics from valid ions
                valid_intensities = list(position_intensity.values())
                if valid_intensities:
                    total_relative = sum(valid_intensities) / max_intensity
                    avg_relative = total_relative / len(valid_intensities)
                    ion_count = len(valid_intensities)
                    std = np.std(valid_intensities / max_intensity) / avg_relative 
                else:
                    total_relative = avg_relative = ion_count = std = 0

                # Store results
                coverage_results.update({
                    f'Total_Relative_Intensity_{chain}_{config["name"]}': total_relative,
                    f'Average_Relative_Intensity_{chain}_{config["name"]}': round(avg_relative, 4),
                    f'Ion_Count_{chain}_{config["name"]}': ion_count,
                    f'SD_{chain}_{config["name"]}': round(std, 4)
                })
            
            # Combine matches from both peptides for coverage calculation
            df[f'ions_{config["name"]}'] = df.index.map(lambda x: ';'.join(matches_A.get(x, []) + matches_B.get(x, [])))
            
            # Special handling for y0 ions using last configuration
            a_y0_mods = {mod: 0 for mod in config['A_mods'] if mod in XL_MODIFICATION_MASSES}
            b_y0_mods = {mod: 0 for mod in config['B_mods'] if mod in XL_MODIFICATION_MASSES}
            
            for matches_X in [matches_A, matches_B]:
                for idx, ions in matches_X.items():
                    intensity = df.at[idx, 'intensity']
                    for label in ions:
                        match = re.match(r'^([AB])([by])(\d+)\+*\[(.*?)]', label)
                        if match:
                            pep_type, ion, pos_str, mod = match.groups()
                            pos = int(pos_str)
                            
                            # Record maximum intensity for y0 ions
                            if pep_type == 'A' and ion == "y" and pos == 0 and mod in a_y0_mods:
                                a_y0_mods[mod] = max(a_y0_mods[mod], intensity)
                            
                            if pep_type == 'B' and ion == "y" and pos == 0 and mod in b_y0_mods:
                                b_y0_mods[mod] = max(b_y0_mods[mod], intensity)
                                
            # Store y0 ion relative intensities
            for mod, intensity in a_y0_mods.items():
                rel_intensity = round(intensity / max_intensity, 4)
                coverage_results[f'A_y0_{mod}_RelInt_{config["name"]}'] = rel_intensity
                
            for mod, intensity in b_y0_mods.items():
                rel_intensity = round(intensity / max_intensity, 4)
                coverage_results[f'B_y0_{mod}_RelInt_{config["name"]}'] = rel_intensity

        # Save labeled spectrum to CSV
        output_csv = output_folder / f"{csv_file.stem}_labeled.csv"
        df.to_csv(output_csv, index=False)
        
        return coverage_results
    except:
        print("Error:",
              f"{row['ScanTitle']}")
        
# Main processing loop - find and process all report files

results = []    

# Process each row in the input file
for _, row in processed_df.iterrows():
    coverage_data = process_row(row, MGF_CSV_PATH, SPECTRUM_OUTPUT)
    if coverage_data:
        results.append(coverage_data)

# Combine results with original data and save
coverage_df = pd.DataFrame(results).reset_index(drop=True)
processed_df = processed_df.reset_index(drop=True)
final_df = pd.concat([processed_df, coverage_df], axis=1)

# Define a function to calculate the 'human_pred' column based on conditions
def calculate_human_pred(row):
    # Check if the Peptide_PSMs column is greater than or equal to 3
    if row['Peptide_PSMs'] >= 3:
        # Check the conditions: 
        if ((row['Ion_Count_A_base_sc'] >= 3 or row['A_y0_sz_RelInt_base_sc'] > 0) and
            (row['Ion_Count_B_base_sc'] >= 3 or row['B_y0_sz_RelInt_base_sc'] > 0)):
            return -1
        else:
            return 0
    else:
        # When Peptide_PSMs < 3, at least three of the following conditions must be satisfied
        count = sum([
            row['Ion_Count_A_base_sc'] >= 3,
            row['A_y0_sz_RelInt_base_sc'] > 0,
            row['Ion_Count_B_base_sc'] >= 3,
            row['B_y0_sz_RelInt_base_sc'] > 0
        ])
        if count >= 3:
            return -1
        else:
            return 0
# Apply the function to each row in the dataframe
final_df['human_pred'] = final_df.apply(calculate_human_pred, axis=1)

# Save the updated dataframe to an Excel file
final_df.to_excel(output_xlsx, index=False)

