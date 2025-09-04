import pandas as pd

# === Load data ===
df_aff = pd.read_csv("ligands_high_affinity.csv")  # columns: ID, affinity

# COM
df_com = pd.read_csv("com_distances.csv")  # columns: ID, COM
#df_com["ID"] = df_com["ID"].str.replace("_docked.pdb", "", regex=False)

# RMSD
df_rmsd = pd.read_csv("rmsd_results.csv") # columns: ID, RMSD, MCS
#df_rmsd["ID"] = df_rmsd["ID"].str.replace("_docked.pdb", "", regex=False)

# PLIP similarity 1MVN
df_plip_1mvn = pd.read_csv("1MVN_PCO_plip_scores.csv")
df_plip_1mvn.rename(columns={"PLIP_score": "PLIP_1MVN"}, inplace=True)
df_plip_1mvn["ID"] = df_plip_1mvn["ID"].str.replace("_docked_complex", "", regex=False)

# PLIP similarity 1QZU
df_plip_1qzu = pd.read_csv("1QZU_pred_plip_scores.csv")
df_plip_1qzu.rename(columns={"PLIP_score": "PLIP_1QZU_pred"}, inplace=True)
df_plip_1qzu["ID"] = df_plip_1qzu["ID"].str.replace("_docked_complex", "", regex=False)

# === Merge ===
df_final = df_aff.merge(df_com, on="ID", how="left") \
                    .merge(df_rmsd, on="ID", how="left") \
                    .merge(df_plip_1mvn[["ID", "PLIP_1MVN"]], on="ID", how="left") \
                    .merge(df_plip_1qzu[["ID", "PLIP_1QZU_pred"]], on="ID", how="left") # left meaning: keep all rows from df_smiles + adds scores

# === Keep only ligands from ligands_higher_affinity_4A.csv ===
df_final = df_final[df_final["ID"].isin(df_aff["ID"])]

# === Remove ligands with missing or zero COM ===
df_final = df_final[df_final["COM"].notna() & (df_final["COM"] != 0)]
df_final = df_final[df_final["RMSD"].notna() & (df_final["RMSD"] != 0)]

# === Replace NaNs in PLIP similarities with 0 ===
df_final["PLIP_1MVN"] = df_final["PLIP_1MVN"].fillna(0)
df_final["PLIP_1QZU_pred"] = df_final["PLIP_1QZU_pred"].fillna(0)

# === Ranking ===
# sort each criteria
df_final["rank_affinity"] = df_final["affinity"].rank(method="min", ascending=True)                   # lower = better
df_final["rank_com"] = df_final["COM"].rank(method="min", ascending=True)                             # lower = better
df_final["rank_rmsd"] = df_final["RMSD"].rank(method="min", ascending=True)                           # lower = better
df_final["rank_mcs"] = df_final["MCS"].rank(method="min", ascending=False)                            # higher = better
df_final["rank_plip_1MVN"] = df_final["PLIP_1MVN"].rank(method="min", ascending=False)            # higher = better
df_final["rank_plip_1QZU_pred"] = df_final["PLIP_1QZU_pred"].rank(method="min", ascending=False)  # higher = better

# what to include in the ranking -> sum up positions in each ranking
df_final["combined_rank"] = df_final[[
    "rank_affinity", "rank_com", "rank_plip_1MVN", "rank_plip_1QZU_pred"
]].sum(axis=1)

# Sort by combined rank
df_final.sort_values("combined_rank", inplace=True)

# Save result
df_final.to_csv("final_ligand_summary.csv", index=False)

print(df_final.head())

