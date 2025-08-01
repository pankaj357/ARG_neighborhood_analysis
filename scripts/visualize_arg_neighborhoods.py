# scripts/visualize_arg_neighborhoods.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the results
df = pd.read_csv("../results/arg_neighborhoods.csv")

# Convert booleans if needed
df["MGE_Nearby"] = df["MGE_Nearby"].astype(bool)
df["Other_ARGs_Nearby"] = df["Other_ARGs_Nearby"].astype(bool)

# Plot 1: ARGs near MGEs
plt.figure(figsize=(10, 6))
sns.countplot(data=df, x="Genome", hue="MGE_Nearby", palette="Set2")
plt.xticks(rotation=90)
plt.title("ARGs with Nearby MGEs per Genome")
plt.ylabel("Number of ARGs")
plt.tight_layout()
plt.savefig("../results/plot_mge_nearby.png")
plt.close()

# Plot 2: ARGs with Nearby Other ARGs
plt.figure(figsize=(10, 6))
sns.countplot(data=df, x="Genome", hue="Other_ARGs_Nearby", palette="Set1")
plt.xticks(rotation=90)
plt.title("ARGs with Other ARGs Nearby (Clusters)")
plt.ylabel("Number of ARGs")
plt.tight_layout()
plt.savefig("../results/plot_arg_clusters.png")
plt.close()

# Summary table
summary = df.groupby("Genome").agg(
    Total_ARGs=("ARG_Locus_Tag", "count"),
    ARGs_with_MGE=("MGE_Nearby", "sum"),
    ARGs_in_Cluster=("Other_ARGs_Nearby", "sum")
).reset_index()

summary["% MGE-Linked"] = (summary["ARGs_with_MGE"] / summary["Total_ARGs"]) * 100
summary["% Clustered"] = (summary["ARGs_in_Cluster"] / summary["Total_ARGs"]) * 100

summary.to_csv("../results/arg_summary_table.csv", index=False)
print("✅ Plots and summary saved in results/")
