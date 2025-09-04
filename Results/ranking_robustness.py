import pandas as pd
from itertools import combinations
from scipy.stats import spearmanr, kendalltau

# === CONFIG ===
input_file = "final_ligand_summary_complete.csv" # complete ranking with all criteria
output_file = "ranking_robustness_with_order.csv"
criteria = ["rank_affinity", "rank_com", "rank_plip_1MVN", "rank_plip_1QZU_pred"] # criteria to use (without MCS, RMSD)
TOP_K = 100 # only test top 100 molecules

# === LOAD DATA ===
df = pd.read_csv(input_file)

# Compute the full ranking (using all criteria)
df["rank_full"] = df[criteria].sum(axis=1)
df = df.sort_values("rank_full").reset_index(drop=True)

top_full = df.head(TOP_K).copy() # get the top 100
top_full_ids = set(top_full["ID"])  # and IDs

results = []

# Try all combinations of criteria (except the full set, which is baseline)
for r in range(1, len(criteria)+1): # range of criteria from 1 (single) to 3 (+1 for full ranking)
    for subset in combinations(criteria, r):
        subset = list(subset)

        # Rank ligands by the chosen subset
        df["rank_subset"] = df[subset].sum(axis=1) # sum up the ranks of all criteria
        df_subset = df.sort_values("rank_subset").reset_index(drop=True) # new index = rank

        top_subset = df_subset.head(TOP_K).copy()
        top_subset_ids = set(top_subset["ID"])

        # Overlap
        overlap_ids = list(top_full_ids & top_subset_ids) # IDs in both sets
        overlap_count = len(overlap_ids)
        percent_overlap = round(100 * overlap_count / TOP_K, 2) # percentage

        # Spearman (only on overlapping ligands)
        if overlap_count > 1:
            full_ranks = top_full.set_index("ID").loc[overlap_ids]["rank_full"] # which IDs are also in the subset
            subset_ranks = top_subset.set_index("ID").loc[overlap_ids]["rank_subset"] # which IDs are also in the full set

            spearman_corr, _ = spearmanr(full_ranks, subset_ranks)
        else:
            spearman_corr = None

        results.append({
            "criteria_used": ", ".join(subset), # which criteria used
            "percent_overlap": percent_overlap, # overlap %
            "spearman_rho": spearman_corr,      # spearman correlation = correlation between order of both rankings
        })

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, index=False)
