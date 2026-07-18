# project3/scripts/validate_scanner.py
# Author: Boris Djagou
# Date: July 19, 2026
# Compare scanner output against the known ground truth genotypes

import csv


def load_ground_truth(path):
    truth = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            truth[row["isolate_id"]] = row["true_genotype"]
    return truth


def load_scanner_results(path):
    results = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mutations = row["mutations_found"]
            called = row["mutation_positions"]
            results[row["isolate_id"]] = row["resistance_status"], row.get("mutation_details", "")
    return results


def main():
    print("\nPROJECT 3 - Scanner Validation Against Ground Truth")
    print("=" * 60)

    truth = load_ground_truth("data/ground_truth.tsv")

    correct = 0
    total = 0

    with open("results/resistance_report.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            isolate_id = row["isolate_id"]
            true_genotype = truth.get(isolate_id, "unknown")
            called_status = row["resistance_status"]
            mutation_details = row["mutation_details"]

            true_is_resistant = true_genotype != "wild-type"
            called_is_resistant = called_status == "RESISTANT"

            match = true_is_resistant == called_is_resistant
            total += 1
            if match:
                correct += 1

            flag = "OK" if match else "MISMATCH"
            print(f"  {isolate_id:<25} true={true_genotype:<12} "
                  f"called={called_status:<12} {flag}")

    accuracy = correct / total * 100
    print(f"\nValidation accuracy: {correct}/{total} ({accuracy:.1f}%)")
    print("=" * 60)


if __name__ == "__main__":
    main()
