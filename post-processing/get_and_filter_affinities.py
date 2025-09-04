import os
import csv
import re
import pandas as pd

# === CONFIGURATION ===
log_folder = 'filtered_ligands_autobox_FMN'   # Folder with *_docking.log
output_csv = 'affinities_autobox_FMN.csv'     # Save all extracted affinities
OUTPUT_FILE = "lig_higher_sub_complete.txt"   # List of ligands above threshold
COMBINED_FOLDER = "combined_complexes_complete"  # Folder containing *_docked_complex.pdb
AFFINITY_THRESHOLD = -6.94                    # Reference affinity to filter against


# === Step 1: Extract affinities from log files ===
results = []
for filename in os.listdir(log_folder):
    if filename.endswith('_docking.log'):
        file_path = os.path.join(log_folder, filename)
        with open(file_path, 'r') as f:
            lines = f.readlines()

            # Extract ID from filename
            ID = filename.split('_docking')[0]

            # Find first docking mode (mode 1)
            for line in lines:
                if re.match(r'^\s*1\s+', line): # get line starting with 1, model 1
                    parts = line.split()
                    if len(parts) >= 2:
                        affinity = float(parts[1]) # get affinity out of column 2
                        results.append((ID, affinity))
                    break


# Save all affinities to CSV
with open(output_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['ID', 'affinity'])
    writer.writerows(results)


# === STEP 2: Filter ligands stronger than threshold ===
df = pd.DataFrame(results, columns=["ID", "affinity"])
filtered = df[df["affinity"] < AFFINITY_THRESHOLD]


# === STEP 3: Write list of complex filenames ===
def write_file_list(output_file, filtered_df):
    with open(output_file, "w") as out:
        for ligand_id in filtered_df["ID"]:
            filename = f"{ligand_id}_docked_complex.pdb"
            out.write(f"{filename}\n")

write_file_list(OUTPUT_FILE, filtered)
