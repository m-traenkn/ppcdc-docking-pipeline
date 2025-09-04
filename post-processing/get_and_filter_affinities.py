import os
import csv
import re
import pandas as pd

# === CONFIGURATION ===
log_folder = 'GninaOut_autobox_FMN'   # Folder with *_docking.log
output_csv_all = 'affinities_autobox_FMN.csv'  # Save all extracted affinities
output_csv_filtered = 'ligands_high_affinity.csv'  # Save filtered ligands with SMILES
smiles_file = 'Components-smiles-stereo-oe.csv'  # CSV with ligand IDs and SMILES
affinity_threshold = -6.94                  # Reference affinity to filter against

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
                if re.match(r'^\s*1\s+', line):
                    parts = line.split()
                    if len(parts) >= 2:
                        affinity = float(parts[1])
                        results.append((ID, affinity))
                    break

# Save all affinities to CSV
df_aff = pd.DataFrame(results, columns=["ID", "affinity"])
df_aff.to_csv(output_csv_all, index=False)

# === Step 2: Filter ligands stronger than threshold ===
df_filtered = df_aff[df_aff["affinity"] < affinity_threshold]

# === Step 3: Merge with SMILES file ===
df_smiles = pd.read_csv(smiles_file)  # Must have columns 'ID' and 'SMILES'
df_high_affinity = df_filtered.merge(df_smiles, on="ID", how="left")

# Reorder columns: ID, SMILES, affinity
df_high_affinity = df_high_affinity[["ID", "SMILES", "affinity"]]

# Save filtered ligands + SMILES to CSV
df_high_affinity.to_csv(output_csv_filtered, index=False)
