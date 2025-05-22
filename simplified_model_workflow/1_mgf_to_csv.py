import os
from pathlib import Path

def parse_mgf(file_path):
    """Parses an MGF file and yields spectra as dictionaries."""
    with open(file_path, 'r') as mgf_file:
        spectrum = {}
        in_spectrum = False
        for line in mgf_file:
            line = line.strip()
            if line == "BEGIN IONS":
                spectrum = {}
                in_spectrum = True
            elif line == "END IONS":
                if in_spectrum:
                    yield spectrum
                in_spectrum = False
            elif in_spectrum:
                if '=' in line:
                    key, value = line.split('=', 1)
                    spectrum[key] = value
                else:
                    if "peaks" not in spectrum:
                        spectrum["peaks"] = []
                    mz, intensity = line.split()
                    spectrum["peaks"].append((float(mz), float(intensity)))

def spectrum_to_csv(spectrum, output_dir):
    """Writes a single spectrum to a CSV file."""
    # Use the TITLE or a counter for the filename
    title = spectrum.get('TITLE', 'spectrum')
    # Sanitize title to ensure it can be used as a filename
    title = title.replace("/", "_").replace("\\", "_").replace(" ", "_")
    csv_filename = os.path.join(output_dir, f"{title}.csv")

    with open(csv_filename, 'w') as csv_file:
        csv_file.write('m/z,intensity\n')
        for mz, intensity in spectrum["peaks"]:
            csv_file.write(f'{mz},{intensity}\n')

def convert_mgf_to_csv(mgf_file_path, output_dir):
    """Converts an MGF file to multiple CSV files, one per spectrum."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for spectrum in parse_mgf(mgf_file_path):
        spectrum_to_csv(spectrum, output_dir)
    
    print(f"Conversion complete. CSV files are saved in '{output_dir}'.")

if __name__ == "__main__":
# Define the input and output directories
    current_dir = Path(__file__).parent
    input_dir = current_dir / 'test_data' 
    output_dir = current_dir / 'test_data' / 'mgf_to_csv'

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over all .mgf files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.mgf'):
            mgf_file_path = os.path.join(input_dir, file_name)
            convert_mgf_to_csv(mgf_file_path, output_dir)
