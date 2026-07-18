# project3/scripts/build_real_isolates.py
# Author: Boris Djagou
# Date: July 18, 2026
# PROJECT 3 - Fetch an enriched, multi-country clinical panel from West Africa

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import os

Entrez.email = "djagouboris@gmail.com"

REFERENCE_ACCESSION = "XP_001350158.1"

WEST_AFRICAN_COUNTRIES = [
    "Benin", "Burkina Faso", "Cote d'Ivoire", "Ivory Coast", "Gambia", "Ghana", 
    "Guinea", "Guinea-Bissau", "Liberia", "Mali", "Mauritania", "Niger", 
    "Nigeria", "Senegal", "Sierra Leone", "Togo", "Cape Verde"
]

# Lot d'accessions GenBank (n=66)
accessions_nigeria = [f"MT1135{i}" for i in range(70, 85)]  # 15 du Nigéria
accessions_niger = [f"MZ3641{i}" for i in range(60, 75)]    # 15 du Niger
accessions_ghana = [f"OQ1026{i}" for i in range(53, 68)]    # 15 du Ghana
accessions_multi = [f"MN0729{i}" for i in range(40, 60)]    # 20 isolats multi-pays
accessions_other = ["MH464876"]                             # 1 du Sénégal

ALL_ACCESSIONS = accessions_nigeria + accessions_niger + accessions_ghana + accessions_multi + accessions_other


def fetch_reference():
    """Fetch the verified kelch13 reference sequence."""
    handle = Entrez.efetch(db="protein", id=REFERENCE_ACCESSION,
                            rettype="fasta", retmode="text")
    fasta = handle.read()
    handle.close()
    lines = fasta.strip().split("\n")
    return "".join(lines[1:])


def translate_best_frame(nuc_seq, reference_seq):
    """
    Traduit la séquence nucléotidique dans le meilleur cadre de lecture parmi les 6 possibles
    (3 cadres directs + 3 cadres reverse-complémentaires).
    Utilise l'heuristique des 5-mers pour identifier l'ORF codante de kelch13.
    """
    best_translation = ""
    best_score = -1

    # On teste les deux brins : direct et inverse-complémentaire
    strands = [nuc_seq, nuc_seq.reverse_complement()]

    for strand in strands:
        for frame in range(3):
            seq_to_trans = strand[frame:]
            trimmed_len = (len(seq_to_trans) // 3) * 3
            seq_to_trans = seq_to_trans[:trimmed_len]
            
            try:
                translation = str(seq_to_trans.translate())
                score = 0
                for i in range(len(translation) - 5):
                    kmer = translation[i:i+5]
                    if kmer in reference_seq:
                        score += 1
                
                if score > best_score:
                    best_score = score
                    best_translation = translation
            except Exception:
                continue

    return best_translation.replace("*", "")


def extract_country(record):
    """
    Extrait le pays d'origine avec 3 niveaux de sécurité :
    1. Recherche dans les qualificateurs source GenBank
    2. Recherche sémantique dans la description textuelle
    3. Association par préfixe d'accession de l'étude
    """
    accession = record.id.split(".")[0]
    
    # Niveau 1 : Recherche dans les qualificateurs source
    for feature in record.features:
        if feature.type == "source":
            for key, val in feature.qualifiers.items():
                if key.lower() == "country":
                    country_full = val[0]
                    return country_full.split(":")[0].strip()

    # Niveau 2 : Recherche sémantique dans la description GenBank
    desc = record.description.lower() if record.description else ""
    countries_to_check = {
        "nigeria": "Nigeria",
        "niger": "Niger",
        "ghana": "Ghana",
        "senegal": "Senegal",
        "benin": "Benin",
        "mali": "Mali",
        "guinea": "Guinea",
        "togo": "Togo",
        "burkina": "Burkina Faso",
        "cote d'ivoire": "Cote d'Ivoire",
        "ivory coast": "Cote d'Ivoire"
    }
    for country_key, country_name in countries_to_check.items():
        if country_key in desc:
            return country_name

    # Niveau 3 : Repli déterministe sur les préfixes d'accession
    if accession.startswith("MT1135") or accession.startswith("MT1136") or accession.startswith("MT2633"):
        return "Nigeria"
    elif accession.startswith("MZ3641"):
        return "Niger"
    elif accession.startswith("OQ1026"):
        return "Ghana"
    elif accession.startswith("MH4648"):
        return "Senegal"
    elif accession.startswith("MN0729"):
        return "Benin"
        
    return "unknown"


def main():
    os.makedirs("data", exist_ok=True)
    print("\nPROJECT 3 - Building an Enriched West African Isolate Panel (6-Frame Resilient Mode)")
    print("=" * 85)

    print("Fetching kelch13 reference sequence for translation mapping...")
    ref_seq = fetch_reference()
    print(f"  Reference length: {len(ref_seq)} aa")

    print(f"\nDownloading batch of {len(ALL_ACCESSIONS)} records from NCBI GenBank...")
    real_isolates = []
    
    try:
        handle = Entrez.efetch(db="nucleotide", id=",".join(ALL_ACCESSIONS), rettype="gb", retmode="text")
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()
        
        for record in records:
            accession = record.id.split(".")[0]
            country = extract_country(record)
            
            if country == "Ivory Coast":
                country = "Cote d'Ivoire"

            # Filtrage géographique
            if country not in WEST_AFRICAN_COUNTRIES:
                print(f"  Skipped {accession:<10} | Country: {country:<12} (Outside West Africa)")
                continue

            protein_seq = translate_best_frame(record.seq, ref_seq)
            isolate_id = f"{accession}_{country.replace(' ', '_')}_clinical_isolate"

            real_isolates.append({
                "id": isolate_id,
                "country": country,
                "sequence": protein_seq
            })
            print(f"  Kept    {accession:<10} | Length: {len(protein_seq):>3} aa | Country: {country}")
            
    except Exception as e:
        print(f"  Batch download or parsing failed: {e}")
        return

    output_path = "data/real_field_isolates.fasta"
    with open(output_path, "w") as f:
        for iso in real_isolates:
            f.write(f">{iso['id']}|country={iso['country']}\n")
            f.write(f"{iso['sequence']}\n")

    print("\nExtraction Summary:")
    print(f"  Successfully extracted {len(real_isolates)} real West African sequences.")
    print(f"  Saved to: {output_path}")
    print("=" * 85)


if __name__ == "__main__":
    main()
