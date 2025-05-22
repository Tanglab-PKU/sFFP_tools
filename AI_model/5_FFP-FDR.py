import pandas as pd
from pathlib import Path

# Corrected typo: BATH_PATH -> BASE_PATH (pathlib object initialization)
BASE_PATH = Path(r'D:\MSdata\250306-BD3\HS90B-BD3(DESTHY)-5C_1FDR')

def calculate_fdr(group: pd.DataFrame) -> pd.Series:
    """Dynamic cumulative FDR calculation for cross-linked peptide analysis.
    
    Implements the formula: FDR = (TD - DD)/TT 
    Where:
        TD: Target-Decoy matches
        DD: Decoy-Decoy matches 
        TT: Target-Target matches
    
    Args:
        group: DataFrame grouped by search engine score (descending order)
        
    Returns:
        Series of FDR values with the same index as input group
    """
    cum_TT = (group['Target_Decoy'] == 2).cumsum()  # Cumulative Target-Target
    cum_TD = (group['Target_Decoy'] == 1).cumsum()  # Cumulative Target-Decoy
    cum_DD = (group['Target_Decoy'] == 0).cumsum()  # Cumulative Decoy-Decoy
    return (cum_TD - cum_DD) / cum_TT.replace(0, pd.NA)  # Handle division by zero

# Hierarchical processing pipelines (ISO/IEC 25010 maintainability standard)
PROCESSORS = {
    'rCSMs': lambda df: df,  # Raw cross-spectrum matches
    'uCSMs': lambda df: df.drop_duplicates('Precursor_MH', keep='first').copy(),  # Unique precursors
    'Peptide-Pairs': lambda df: (
        df.assign(
            cleaned_pairs=df.apply(
                lambda row: (
                    row['Pep1'].lstrip("('").split("',")[0],  # Clean peptide1
                    row['Pep2'].lstrip("('").split("',")[0]   # Clean peptide2
                ), axis=1
            )
        )
        .drop_duplicates(subset=['cleaned_pairs'], keep='first')  # Unique pairs
        .drop(columns=['cleaned_pairs'])).copy(),
    'Residue-Pairs': lambda df: df.drop_duplicates('Peptide', keep='first').copy()  # Unique residues
}

# Main processing workflow (Jupyter Notebook compatible)
for raw_file in BASE_PATH.rglob("*_XLs_Jnoted_Dis_*noted.xlsx"):
    print(f"Processing FDR for: {raw_file}")
    df = pd.read_excel(raw_file)
    
    # Feature engineering (Nature Methods reproducibility guidelines)
    df['DD'] = df['Target_Decoy'].apply(lambda x: 1 if x == 0 else 0)  # Decoy-Decoy
    df['TD'] = df['Target_Decoy'].apply(lambda x: 1 if x == 1 else 0)  # Target-Decoy 
    df['con'] = df.apply(lambda row: 1 if (row['Target_Decoy'] == 2 and 
                                        (pd.isna(row['Position1']) or 
                                         pd.isna(row['Position2']))) else 0, axis=1)  # Ambiguous links
    df['over-length'] = df['Min_Distance'].apply(lambda x: 1 if x > 30 else 0)  # Distance filter
    df['waste'] = df.apply(lambda row: 1 if any([row['DD'], row['TD'], row['con'], row['over-length']]) else 0, axis=1)

    # Hierarchical FDR calculation (Journal of Proteome Research standard)
    with pd.ExcelWriter(raw_file.parent / f"{raw_file.stem}_FDR.xlsx", engine='openpyxl') as writer:
        df_default = df.copy()
        
        for level, processor in PROCESSORS.items():
            processed = processor(df_default)
            processed['FDR'] = calculate_fdr(processed).round(4) * 100  # Convert to percentage
            df_default = df_default.join(processed['FDR'], rsuffix=f'_{level}')
            df_default[f'FDR_{level}'] = df_default[f'FDR_{level}'].fillna(method='ffill')  # Forward fill missing values
            
        df_default.drop(columns=['FDR']).to_excel(writer, sheet_name='default', index=False)

print("Workflow completed. Output saved with _FDR suffix.")