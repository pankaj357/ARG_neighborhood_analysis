# scripts/parse_amr_gff.py

import os
import pandas as pd
from collections import defaultdict

# Constants
AMR_FOLDER = "../amrfinder_output"
GFF_FOLDER = "../prokka_annotations"
OUTPUT_FILE = "../results/arg_neighborhoods.csv"
UPSTREAM_DOWNSTREAM = 10
MGE_KEYWORDS = ['transposase', 'integrase', 'recombinase', 'mobile', 'IS', 'phage']

# Collect all results
all_rows = []

# Iterate over each genome
for amr_file in sorted(os.listdir(AMR_FOLDER)):
    if not amr_file.endswith(".txt"):
        continue

    genome_id = amr_file.replace("_amrfinder.txt", "")
    gff_file = os.path.join(GFF_FOLDER, genome_id + ".gff")
    amr_path = os.path.join(AMR_FOLDER, amr_file)

    print(f"Processing: {genome_id}")

    # Read AMR gene list
    amr_df = pd.read_csv(amr_path, sep="\t", comment='#')
    amr_locus_tags = set(amr_df['Protein id'].fillna("").tolist())
    locus_to_symbol = dict(zip(amr_df['Protein id'], amr_df['Element symbol']))

    # Parse GFF and store gene features in order
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attr_raw = parts[8]
            attrs = {kv.split("=")[0]: kv.split("=")[1] for kv in attr_raw.split(";") if "=" in kv}
            locus_tag = attrs.get("locus_tag", "")
            product = attrs.get("product", "")
            gene = attrs.get("gene", "")
            genes.append({
                "locus_tag": locus_tag,
                "gene": gene,
                "product": product,
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6]
            })

    # Build lookup table
    locus_to_index = {g['locus_tag']: i for i, g in enumerate(genes)}

    # Analyze neighborhoods
    for locus in amr_locus_tags:
        if locus not in locus_to_index:
            continue
        idx = locus_to_index[locus]

        start_idx = max(0, idx - UPSTREAM_DOWNSTREAM)
        end_idx = min(len(genes), idx + UPSTREAM_DOWNSTREAM + 1)
        neighborhood = genes[start_idx:end_idx]

        mge_nearby = any(any(k.lower() in (g['product'] or '').lower() for k in MGE_KEYWORDS) for g in neighborhood)
        other_args = any(g['locus_tag'] in amr_locus_tags and g['locus_tag'] != locus for g in neighborhood)

        all_rows.append({
            "Genome": genome_id,
            "ARG_Locus_Tag": locus,
            "ARG_Symbol": locus_to_symbol.get(locus, ""),
            "MGE_Nearby": mge_nearby,
            "Other_ARGs_Nearby": other_args,
            "Neighbor_Locus_Tags": ",".join([g['locus_tag'] for g in neighborhood]),
            "Neighbor_Products": ",".join([g['product'] for g in neighborhood]),
        })

# Save output
df = pd.DataFrame(all_rows)
df.to_csv(OUTPUT_FILE, index=False)
print(f"\n✅ Output saved: {OUTPUT_FILE}")
