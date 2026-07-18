# project3/scripts/build_field_isolates.py
# Author: Boris Djagou
# Date: July 19, 2026
# PROJECT 3 - Build a realistic field isolate dataset with country metadata
# Simulates MalariaGEN-style kelch13 field isolate sequences

from Bio import Entrez
import random
import os

Entrez.email = "djagouboris@gmail.com"

REFERENCE_ACCESSION = "XP_001350158.1"

# WHO-validated resistance mutations with real geographic prevalence data
RESISTANCE_MUTATIONS = [
    {"name": "F446I", "pos": 446, "wt": "F", "mut": "I", "primary_region": "Myanmar"},
    {"name": "N458Y", "pos": 458, "wt": "N", "mut": "Y", "primary_region": "Southeast Asia"},
    {"name": "M476I", "pos": 476, "wt": "M", "mut": "I", "primary_region": "Southeast Asia"},
    {"name": "Y493H", "pos": 493, "wt": "Y", "mut": "H", "primary_region": "Southeast Asia"},
    {"name": "R539T", "pos": 539, "wt": "R", "mut": "T", "primary_region": "Southeast Asia / emerging Africa"},
    {"name": "I543T", "pos": 543, "wt": "I", "mut": "T", "primary_region": "Southeast Asia"},
    {"name": "P553L", "pos": 553, "wt": "P", "mut": "L", "primary_region": "Southeast Asia"},
    {"name": "R561H", "pos": 561, "wt": "R", "mut": "H", "primary_region": "Rwanda / DR Congo"},
    {"name": "C580Y", "pos": 580, "wt": "C", "mut": "Y", "primary_region": "Global / dominant"},
]

# West African countries for isolate metadata (per Project 21 population panel)
WEST_AFRICAN_COUNTRIES = [
    "Benin", "Burkina Faso", "Ghana", "Nigeria",
    "Mali", "Senegal", "Cote d'Ivoire", "Togo",
]


def fetch_reference():
    """Fetch the verified kelch13 reference sequence."""
    handle = Entrez.efetch(db="protein", id=REFERENCE_ACCESSION,
                            rettype="fasta", retmode="text")
    fasta = handle.read()
    handle.close()
    lines = fasta.strip().split("\n")
    return "".join(lines[1:])


def apply_mutations(reference_seq, mutation_list):
    """Apply a list of mutations to build one isolate sequence."""
    seq = list(reference_seq)
    for mut in mutation_list:
        pos = mut["pos"] - 1
        seq[pos] = mut["mut"]
    return "".join(seq)


def build_isolate_panel(reference_seq, n_isolates=30, seed=42):
    """
    Build a panel of simulated field isolates.
    Most isolates are wild-type (no resistance mutation), consistent
    with real-world low resistance prevalence in West Africa.
    A minority carry 1 resistance mutation, reflecting emerging
    R561H and C580Y detection in the region.
    """
    random.seed(seed)
    isolates = []

    for i in range(n_isolates):
        country = random.choice(WEST_AFRICAN_COUNTRIES)
        isolate_id = f"ISO{i+1:03d}_{country.replace(' ', '_').replace(chr(39), '')}"

        # 80 percent wild-type, 20 percent carrying a resistance mutation
        # This mirrors real surveillance data where resistance is still
        # emerging and not yet widespread in West Africa
        if random.random() < 0.20:
            mutation = random.choice(RESISTANCE_MUTATIONS)
            seq = apply_mutations(reference_seq, [mutation])
            genotype = mutation["name"]
        else:
            seq = reference_seq
            genotype = "wild-type"

        isolates.append({
            "id": isolate_id,
            "country": country,
            "sequence": seq,
            "true_genotype": genotype,
        })

    return isolates


def main():
    os.makedirs("data", exist_ok=True)

    print("\nPROJECT 3 - Building Field Isolate Panel")
    print("=" * 60)

    print(f"\nFetching kelch13 reference: {REFERENCE_ACCESSION}")
    reference_seq = fetch_reference()
    print(f"  Reference length: {len(reference_seq)} aa")

    print(f"\nBuilding panel of 30 field isolates across West Africa...")
    isolates = build_isolate_panel(reference_seq, n_isolates=30)

    # Save FASTA
    fasta_path = "data/field_isolates.fasta"
    with open(fasta_path, "w") as f:
        for iso in isolates:
            f.write(f">{iso['id']}|country={iso['country']}\n")
            f.write(f"{iso['sequence']}\n")

    # Save ground truth metadata (for validation later, not used by the scanner)
    truth_path = "data/ground_truth.tsv"
    with open(truth_path, "w") as f:
        f.write("isolate_id\tcountry\ttrue_genotype\n")
        for iso in isolates:
            f.write(f"{iso['id']}\t{iso['country']}\t{iso['true_genotype']}\n")

    print(f"\nPanel summary:")
    genotype_counts = {}
    for iso in isolates:
        genotype_counts[iso["true_genotype"]] = genotype_counts.get(iso["true_genotype"], 0) + 1

    for genotype, count in sorted(genotype_counts.items()):
        print(f"  {genotype:<12} {count:>3} isolates")

    print(f"\nSaved: {fasta_path}")
    print(f"Saved: {truth_path} (ground truth, for scanner validation)")
    print("=" * 60)


if __name__ == "__main__":
    main()
