# project3/scripts/scan_resistance.py
# Author: Boris Djagou
# Date: July 18, 2026
# PROJECT 3 - Artemisinin Resistance Mutation Scanner
# BLAST-based pipeline to screen kelch13 field isolates for WHO-validated mutations

import subprocess
import os
import csv
from datetime import datetime

REFERENCE_FASTA = "data/reference.fasta"
ISOLATES_FASTA = "data/real_field_isolates.fasta"
DB_DIR = "results/blast_db"
BLAST_OUT = "results/blast_hits.tsv"

RESISTANCE_POSITIONS = {
    446: {"wt": "F", "name": "F446I", "region": "Myanmar"},
    458: {"wt": "N", "name": "N458Y", "region": "Southeast Asia"},
    476: {"wt": "M", "name": "M476I", "region": "Southeast Asia"},
    493: {"wt": "Y", "name": "Y493H", "region": "Southeast Asia"},
    539: {"wt": "R", "name": "R539T", "region": "Southeast Asia / emerging Africa"},
    543: {"wt": "I", "name": "I543T", "region": "Southeast Asia"},
    553: {"wt": "P", "name": "P553L", "region": "Southeast Asia"},
    561: {"wt": "R", "name": "R561H", "region": "Rwanda / DR Congo"},
    580: {"wt": "C", "name": "C580Y", "region": "Global / dominant"},
}

PROPELLER_DOMAIN = (440, 680)


def prepare_reference():
    """Extract the reference sequence from the isolate panel building script's output."""
    from Bio import Entrez
    Entrez.email = "djagouboris@gmail.com"
    handle = Entrez.efetch(db="protein", id="XP_001350158.1",
                            rettype="fasta", retmode="text")
    fasta = handle.read()
    handle.close()

    os.makedirs("data", exist_ok=True)
    with open(REFERENCE_FASTA, "w") as f:
        f.write(">kelch13_reference_3D7\n")
        lines = fasta.strip().split("\n")
        f.write("".join(lines[1:]) + "\n")

    with open(REFERENCE_FASTA) as f:
        lines = f.read().strip().split("\n")
    return "".join(lines[1:])


def build_blast_db():
    """Build a BLAST database from the reference sequence."""
    os.makedirs(DB_DIR, exist_ok=True)
    cmd = ["makeblastdb", "-in", REFERENCE_FASTA, "-dbtype", "prot",
           "-out", f"{DB_DIR}/kelch13_ref"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def run_blast():
    """
    BLAST all isolates against the reference.
    We request custom output formats including qseq (query alignment) 
    and sseq (subject alignment) to perform position-specific coordinate mapping.
    """
    cmd = [
        "blastp", "-query", ISOLATES_FASTA, "-db", f"{DB_DIR}/kelch13_ref",
        "-out", BLAST_OUT, "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qseq sseq",
        "-evalue", "1e-5"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stderr


def parse_isolate_fasta(filepath):
    """Parse the field isolate FASTA, extracting country from header."""
    isolates = {}
    current_id = None
    current_country = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    isolates[current_id] = {
                        "country": current_country,
                        "sequence": "".join(current_seq),
                    }
                header = line[1:]
                current_id = header.split("|")[0]
                current_country = header.split("country=")[1] if "country=" in header else "unknown"
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            isolates[current_id] = {
                "country": current_country,
                "sequence": "".join(current_seq),
            }

    return isolates


def map_and_scan_mutations(blast_hsp):
    """
    Map reference positions to query coordinates using BLAST HSP alignment string traversal.
    This resolves issues caused by sequence truncation, insertions, or deletions (indels).
    """
    mutations_found = []
    
    # Extraction des coordonnées et chaînes d'alignement de l'HSP (High-scoring Segment Pair)
    sstart = int(blast_hsp["sstart"])
    send = int(blast_hsp["send"])
    qseq = blast_hsp["qseq"]
    sseq = blast_hsp["sseq"]
    
    for ref_pos, info in RESISTANCE_POSITIONS.items():
        # Vérification si la position d'intérêt est couverte par l'alignement
        if not (sstart <= ref_pos <= send):
            continue  # Position non couverte par le fragment séquencé

        # Parcours des chaînes alignées pour mapper la position de la référence (subject)
        current_ref_pos = sstart
        observed_residue = None

        for i in range(len(sseq)):
            char_s = sseq[i]
            char_q = qseq[i]

            if char_s != '-':
                if current_ref_pos == ref_pos:
                    observed_residue = char_q
                    break
                current_ref_pos += 1
            # Les gaps dans la référence (insertions dans la requête) ne font pas avancer la position de la référence

        # Vérification si l'acide aminé observé dévie du type sauvage (et n'est pas un gap)
        expected_wt = info["wt"]
        if observed_residue and observed_residue != expected_wt and observed_residue != "-":
            mutations_found.append({
                "position": ref_pos,
                "wild_type": expected_wt,
                "observed": observed_residue,
                "mutation_name": info["name"],
                "known_region": info["region"],
            })

    return mutations_found


def generate_csv_report(results, output_path):
    """Write the per-isolate resistance report as CSV."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "isolate_id", "country", "blast_identity_pct",
            "blast_evalue", "resistance_status",
            "mutations_found", "mutation_positions", "mutation_details"
        ])
        for r in results:
            mutation_names = ";".join(m["mutation_name"] for m in r["mutations"])
            mutation_positions = ";".join(str(m["position"]) for m in r["mutations"])
            mutation_details = ";".join(
                f"{m['wild_type']}{m['position']}{m['observed']}({m['known_region']})"
                for m in r["mutations"]
            )
            status = "RESISTANT" if r["mutations"] else "wild-type"

            writer.writerow([
                r["isolate_id"], r["country"], r["identity"],
                r["evalue"], status,
                len(r["mutations"]), mutation_positions, mutation_details
            ])


def generate_html_report(results, output_path):
    """Write a per-isolate resistance report as a styled HTML file."""
    n_total = len(results)
    n_resistant = sum(1 for r in results if r["mutations"])
    n_wildtype = n_total - n_resistant

    country_counts = {}
    for r in results:
        country_counts[r["country"]] = country_counts.get(r["country"], 0) + 1

    rows_html = ""
    for r in sorted(results, key=lambda x: (len(x["mutations"]) == 0, x["isolate_id"])):
        status = "RESISTANT" if r["mutations"] else "wild-type"
        status_class = "resistant" if r["mutations"] else "wildtype"
        mutations_str = ", ".join(
            f"{m['wild_type']}{m['position']}{m['observed']}"
            for m in r["mutations"]
        ) if r["mutations"] else "-"

        identity_val = f"{r['identity']}%" if r['identity'] != "N/A" else "N/A"

        rows_html += f"""
        <tr class="{status_class}">
            <td>{r['isolate_id']}</td>
            <td>{r['country']}</td>
            <td>{identity_val}</td>
            <td>{r['evalue']}</td>
            <td class="status">{status}</td>
            <td>{mutations_str}</td>
        </tr>"""

    country_rows = "".join(
        f"<tr><td>{c}</td><td>{n}</td></tr>"
        for c, n in sorted(country_counts.items())
    )

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>Artemisinin Resistance Scan Report</title>
<style>
    body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #0D1B2A; color: #E2E8F0; margin: 0; padding: 40px; }}
    .container {{ max-width: 1100px; margin: 0 auto; }}
    h1 {{ color: #FFFFFF; border-bottom: 3px solid #0D9488; padding-bottom: 12px; }}
    .subtitle {{ color: #94A3B8; margin-bottom: 30px; }}
    .stats {{ display: flex; gap: 20px; margin-bottom: 30px; }}
    .stat-box {{ background: #1E293B; border-radius: 10px; padding: 20px; flex: 1; text-align: center; border: 1px solid #334155; }}
    .stat-num {{ font-size: 36px; font-weight: bold; color: #0D9488; }}
    .stat-label {{ font-size: 13px; color: #94A3B8; text-transform: uppercase; letter-spacing: 1px; margin-top: 5px; }}
    .stat-box.resistant .stat-num {{ color: #DC2626; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 20px; background: #1E293B; border-radius: 8px; overflow: hidden; }}
    th {{ background: #0D9488; color: white; padding: 12px 15px; text-align: left; font-size: 13px; text-transform: uppercase; }}
    td {{ padding: 10px 15px; border-bottom: 1px solid #334155; font-size: 14px; }}
    tr.resistant {{ background: rgba(220, 38, 38, 0.08); }}
    tr.resistant .status {{ color: #EF4444; font-weight: bold; }}
    tr.wildtype .status {{ color: #10B981; }}
    .section-title {{ color: #0D9488; margin-top: 40px; font-size: 18px; }}
    .footer {{ margin-top: 40px; color: #64748B; font-size: 12px; border-top: 1px solid #334155; padding-top: 20px; }}
</style>
</head>
<body>
<div class="container">
    <h1>Artemisinin Resistance Mutation Scanner</h1>
    <div class="subtitle">Project 3 - kelch13 Field Isolate Screening Report</div>

    <div class="stats">
        <div class="stat-box">
            <div class="stat-num">{n_total}</div>
            <div class="stat-label">Isolates Screened</div>
        </div>
        <div class="stat-box resistant">
            <div class="stat-num">{n_resistant}</div>
            <div class="stat-label">Resistant Isolates</div>
        </div>
        <div class="stat-box">
            <div class="stat-num">{n_wildtype}</div>
            <div class="stat-label">Wild-type Isolates</div>
        </div>
        <div class="stat-box">
            <div class="stat-num">{len(country_counts)}</div>
            <div class="stat-label">Countries Represented</div>
        </div>
    </div>

    <div class="section-title">Per-Isolate Results</div>
    <table>
        <tr>
            <th>Isolate ID</th>
            <th>Country</th>
            <th>BLAST Identity</th>
            <th>E-value</th>
            <th>Status</th>
            <th>Mutations Found</th>
        </tr>
        {rows_html}
    </table>

    <div class="section-title">Isolates per Country</div>
    <table>
        <tr><th>Country</th><th>Isolate Count</th></tr>
        {country_rows}
    </table>

    <div class="footer">
        Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}<br>
        Igbega X Bioinformatics Group - MSc Bioinformatics Self-Training Program<br>
        Reference: kelch13 P. falciparum 3D7 (XP_001350158.1)<br>
        Screened positions: 9 WHO-validated resistance mutations (propeller domain 440-680)
    </div>
</div>
</body>
</html>"""

    with open(output_path, "w") as f:
        f.write(html)


def main():
    print("\nPROJECT 3 - Artemisinin Resistance Mutation Scanner (Alignment-Aware)")
    print("=" * 65)

    print("\n[1/5] Preparing reference sequence...")
    reference_seq = prepare_reference()
    print(f"  Reference: {len(reference_seq)} aa")

    print("\n[2/5] Building BLAST database...")
    if not build_blast_db():
        print("  Failed to build BLAST database")
        return
    print(f"  Database built: {DB_DIR}/kelch13_ref")

    print("\n[3/5] Running BLAST search on all field isolates...")
    ok, err = run_blast()
    if not ok:
        print(f"  BLAST failed: {err[:300]}")
        return
    print(f"  BLAST hits saved: {BLAST_OUT}")

    print("\n[4/5] Parsing isolates and scanning for resistance mutations...")
    isolates = parse_isolate_fasta(ISOLATES_FASTA)

    # Lecture des données BLAST avec l'inclusion de qseq et sseq (champs 12 et 13)
    blast_data = {}
    with open(BLAST_OUT) as f:
        for line in f:
            fields = line.strip().split("\t")
            qseqid = fields[0].split("|")[0]
            pident = fields[2]
            evalue = fields[10]
            
            # Stockage de toutes les informations HSP nécessaires à la cartographie des coordonnées
            if qseqid not in blast_data:
                blast_data[qseqid] = {
                    "identity": pident,
                    "evalue": evalue,
                    "qstart": fields[6],
                    "qend": fields[7],
                    "sstart": fields[8],
                    "send": fields[9],
                    "qseq": fields[12],
                    "sseq": fields[13]
                }

    results = []
    for isolate_id, data in isolates.items():
        # Cartographie des mutations basée sur l'HSP BLAST ou repli sur l'approche par défaut si aucun hit
        blast_info = blast_data.get(isolate_id)
        
        if blast_info:
            mutations = map_and_scan_mutations(blast_info)
            identity = blast_info["identity"]
            evalue = blast_info["evalue"]
        else:
            # Cas d'exclusion : l'isolat ne s'aligne pas du tout sur kelch13 (hors cible / mauvaise qualité)
            mutations = []
            identity = "N/A"
            evalue = "N/A"

        results.append({
            "isolate_id": isolate_id,
            "country": data["country"],
            "sequence": data["sequence"],
            "mutations": mutations,
            "identity": identity,
            "evalue": evalue,
        })

        status = "RESISTANT" if mutations else "wild-type"
        mut_str = ", ".join(m["mutation_name"] for m in mutations) if mutations else "-"
        print(f"  {isolate_id:<25} {data['country']:<15} {status:<10} {mut_str}")

    print("\n[5/5] Generating reports...")
    csv_path = "results/resistance_report.csv"
    html_path = "results/resistance_report.html"

    generate_csv_report(results, csv_path)
    generate_html_report(results, html_path)

    n_resistant = sum(1 for r in results if r["mutations"])
    print(f"\n  Total isolates screened: {len(results)}")
    print(f"  Resistant isolates found: {n_resistant}")
    print(f"  CSV report: {csv_path}")
    print(f"  HTML report: {html_path}")

    print("\nPROJECT 3 pipeline complete")
    print("=" * 65)


if __name__ == "__main__":
    main()
