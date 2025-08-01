import os
import csv
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# Define paths
GFF_DIR = "../prokka_annotations"
FASTA_DIR = "../genomes"
AMR_DIR = "../amrfinder_output"
OUT_PATH = "../results/arg_gc_content.csv"

# Helper: Parse GFF to extract gene coords
def parse_gff(gff_path):
    genes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            cols = line.strip().split("\t")
            start, end, strand = int(cols[3]), int(cols[4]), cols[6]
            gene_id = None
            for entry in cols[8].split(";"):
                if entry.startswith("ID="):
                    gene_id = entry[3:]
            if gene_id:
                genes[gene_id] = {"start": start, "end": end, "strand": strand}
    return genes

# Helper: Load genome sequence
def load_genome(fna_path):
    return next(SeqIO.parse(fna_path, "fasta")).seq

# Helper: Parse AMRFinder txt to get ARGs
# ✅ FIXED version — USE this
def parse_amr(amr_path):
    arg_loci = []
    with open(amr_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Element"):
                continue  # Skip comments and column header
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            try:
                gene = cols[5]
                start = int(cols[2])
                end = int(cols[3])
                arg_loci.append((gene, start, end))
            except ValueError:
                continue
    return arg_loci

# Main
os.makedirs("../results", exist_ok=True)
header = ["Strain", "ARG", "Start", "End", "GC_ARG", "GC_Genome", "GC_FlankAvg"]
with open(OUT_PATH, "w", newline="") as out_f:
    writer = csv.writer(out_f)
    writer.writerow(header)

    for filename in os.listdir(AMR_DIR):
        if not filename.endswith("_amrfinder.txt"):
            continue

        strain = filename.replace("_amrfinder.txt", "")
        print(f"🔬 Processing: {strain}")

        amr_path = os.path.join(AMR_DIR, filename)
        gff_path = os.path.join(GFF_DIR, f"{strain}.gff")
        fna_path = os.path.join(FASTA_DIR, f"{strain}.fna")

        if not os.path.exists(gff_path) or not os.path.exists(fna_path):
            print(f"⚠️ Missing files for {strain}")
            continue

        genome = load_genome(fna_path)
        genome_gc = round(gc_fraction(genome) * 100, 2)
        gene_dict = parse_gff(gff_path)
        arg_loci = parse_amr(amr_path)

        for arg, start, end in arg_loci:
            # Get ARG sequence
            arg_seq = genome[start-1:end]
            gc_arg = round(gc_fraction(arg_seq) * 100, 2)

            # Find flanking genes (optional: ±2)
            flanking_gcs = []
            for gene_id, data in gene_dict.items():
                if abs(data["start"] - start) < 3000 or abs(data["end"] - end) < 3000:
                    seq = genome[data["start"]-1:data["end"]]
                    flanking_gcs.append(gc_fraction(seq) * 100)

            gc_flank_avg = round(sum(flanking_gcs)/len(flanking_gcs), 2) if flanking_gcs else "NA"

            writer.writerow([strain, arg, start, end, gc_arg, genome_gc, gc_flank_avg])

