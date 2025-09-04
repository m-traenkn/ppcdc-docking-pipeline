"""
Molecule Filtering and Clustering

Steps:
1. Load molecules from CSV (basic SMILES validity + length check)
2. Generate Morgan fingerprints for all molecules
3. Remove near-duplicates using Tanimoto similarity
4. Cluster molecules with Butina clustering
5. Pick a cluster representative (centroid)
6. Write cluster representatives to an output CSV file

Dependencies:
- RDKit
"""

import csv
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina

# === Configuration ===
input_file = "Components-smiles-stereo-oe.csv"
output_file = "filtered_input.csv"

sim_theshold = 0.95     # threshold for near-duplicate removal
cluster_distance_cutoff = 0.65  # cutoff for Butina clustering

# === Data Loading ===
def load_csv(path, smiles_col=0, id_col=1, has_header=True):
    """
    Load molecules from a CSV file.

    Args:
        path (str): Path to CSV file
        smiles_col (int): Column index containing SMILES strings
        id_col (int): Column index containing molecule IDs
        has_header (bool): Skip first row if True

    Returns:
        mols (list[Chem.Mol]): RDKit molecule objects
        ids (list[str]): Molecule identifiers
        smis (list[str]): Original SMILES strings
    """
    mols, ids, smis = [], [], []
    with open(path) as f:
        reader = csv.reader(f)
        if has_header:
            next(reader)  # skip header row
        for row in reader:
            try:
                smi = row[smiles_col].strip()
                # discard overly short/long SMILES
                if 5 <= len(smi) <= 70:
                    mol = Chem.MolFromSmiles(smi)
                    if mol:
                        mols.append(mol)
                        smis.append(smi)
                        ids.append(row[id_col])
            except IndexError:
                # skip rows that don’t have enough columns
                continue
    return mols, ids, smis

# === Fingerprint Calculation ===
def calc_fps(mols, radius=2, nBits=1024):
    return [AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits) for m in mols]


# === Duplicate Removal ===
def filter_similar(fps, smis, ids, sim_thresh=0.95):
    """
    Remove near-duplicate molecules based on Tanimoto similarity

    Args:
        fps (list[ExplicitBitVect]): Molecule fingerprints
        smis (list[str]): SMILES strings
        ids (list[str]): Molecule IDs
        sim_thresh (float): Tanimoto similarity threshold

    Returns:
        keep_fps, keep_smis, keep_ids (lists): Filtered molecules
    """
    keep_fps, keep_smis, keep_ids = [], [], []
    for fp, smi, mid in zip(fps, smis, ids):
        if not any(DataStructs.TanimotoSimilarity(fp, kfp) >= sim_thresh for kfp in keep_fps): # similarity larger than thehsold?
            keep_fps.append(fp)
            keep_smis.append(smi)
            keep_ids.append(mid)
    return keep_fps, keep_smis, keep_ids


# === Representative Selection ===
def pick_representative(fps, cluster, ids, smis):
    """
    Pick representative (centroid) molecule from a cluster

    Args:
        fps (list[ExplicitBitVect]): Fingerprints
        cluster (list[int]): Indices of molecules in cluster
        ids (list[str]): Molecule IDs
        smis (list[str]): SMILES strings

    Returns:
        int: Index of representative molecule
    """
    if len(cluster) == 1: # only one mol in cluster
        return cluster[0]

    best_idx, best_mean = None, -1.0 # start with worse mean sim
    for i in cluster:
        sims = [DataStructs.TanimotoSimilarity(fps[i], fps[j]) for j in cluster if i != j] # Tanimoto sim for each mol i again each mol j
        mean_sim = sum(sims) / len(sims)
        if mean_sim > best_mean:
            best_mean, best_idx = mean_sim, i # keep if mean is bettern than so far
    return best_idx

# === Output Writing ===
def write_csv(outfile, reps):
    """
    Write representative molecules to CSV

    Args:
        outfile (str): Output file path
        reps (list[tuple]): List of (ID, SMILES) pairs
    """
    with open(outfile, 'w', newline='') as w:
        writer = csv.writer(w)
        writer.writerow(['ID', 'SMILES'])
        writer.writerows(reps)
    print(f"Wrote {len(reps)} representatives to {outfile}")
    

# === Workflow ===
def run_filtering():
    """Run the full pipeline: load → filter → cluster → select representatives → save."""
    mols, ids, smis = load_csv(input_file)
    if not mols:
        print("No valid molecules found")
        return

    # Step 1: Generate fingerprints
    fps = calc_fps(mols)

    # Step 2: Remove near-duplicates
    fps_f, smis_f, ids_f = filter_similar(fps, smis, ids, sim_thresh=sim_theshold)

    # Step 3: Cluster molecules
    clusters = butina_cluster(fps_f, cutoff=cluster_distance_cutoff)

    # Step 4: Pick representatives from each cluster
    reps = []
    for cluster in clusters:
        rep_idx = pick_representative(fps_f, cluster, ids_f, smis_f)
        reps.append((ids_f[rep_idx], smis_f[rep_idx]))

    # Step 5: Write representatives to file
    write_csv(output_file, reps)


if __name__ == '__main__':
    run_filtering()