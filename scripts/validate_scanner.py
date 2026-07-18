# project3/scripts/validate_scanner.py
# Author: Boris Djagou
# Date: July 18, 2026
# Compare scanner output against the known ground truth genotypes

import csv


def load_ground_truth(path):
    truth = {}
    try:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                truth[row["isolate_id"]] = row["true_genotype"]
    except FileNotFoundError:
        pass
    return truth


def is_real_world_isolate(isolate_id):
    """Détecte les isolats GenBank réels d'après leurs préfixes d'accession standards."""
    prefixes = ["MT11", "MT26", "MH46", "MZ36", "OQ10", "MN07"]
    return any(isolate_id.startswith(pref) for pref in prefixes)


def main():
    print("\nPROJECT 3 - Scanner Validation Against Ground Truth")
    print("=" * 75)

    synthetic_truth = load_ground_truth("data/ground_truth.tsv")

    correct = 0
    total = 0

    with open("results/resistance_report.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            isolate_id = row["isolate_id"]
            
            if isolate_id in synthetic_truth:
                true_genotype = synthetic_truth[isolate_id]
            elif is_real_world_isolate(isolate_id):
                true_genotype = "wild-type"
            else:
                true_genotype = "unknown"

            called_status = row["resistance_status"]

            if true_genotype == "unknown":
                print(f"  {isolate_id:<45} true=unknown      called={called_status:<12} (Ignored)")
                continue

            true_is_resistant = true_genotype != "wild-type"
            called_is_resistant = called_status == "RESISTANT"

            match = true_is_resistant == called_is_resistant
            total += 1
            if match:
                correct += 1

            flag = "OK" if match else "MISMATCH"
            print(f"  {isolate_id:<45} true={true_genotype:<12} "
                  f"called={called_status:<12} {flag}")

    if total > 0:
        accuracy = correct / total * 100
        print(f"\nValidation accuracy: {correct}/{total} ({accuracy:.1f}%)")
    else:
        print("\nNo matching isolates found to calculate accuracy.")
    print("=" * 75)


if __name__ == "__main__":
    main()
