from Bio import Entrez, SeqIO
import time
from pathlib import Path
import pandas as pd


Entrez.email = "runchen.yang@mail.utoronto.ca"


def read_species_list(file_path: Path) -> list[dict]:
    """Read the species spreadsheet and return species records."""
    df = pd.read_excel(file_path)

    print("Columns found:")
    print(df.columns.tolist())

    species_records = []

    for _, row in df.iterrows():
        species_records.append({
            "taxonomic_group": row["Taxonomic group"],
            "common_name": row["Species (Common Name)"],
            "scientific_name": row["Species (Scientific Name)"],
        })

    return species_records

def search_ncbi_candidates(gene: str, species: str, max_results: int = 10) -> list[dict]:
    """Search NCBI Nucleotide for mRNA records matching one gene and one species."""
    query = f'{gene}[Gene Name] AND "{species}"[Organism] AND mRNA'

    print(f"Searching: {query}")

    with Entrez.esearch(db="nucleotide", term=query, retmax=max_results) as handle:
        search_record = Entrez.read(handle)

    ids = search_record["IdList"]

    if not ids:
        return []

    with Entrez.esummary(db="nucleotide", id=",".join(ids)) as handle:
        summaries = Entrez.read(handle)

    candidates = []

    for item in summaries:
        accession = str(item.get("AccessionVersion", ""))
        title = str(item.get("Title", ""))
        length = int(item.get("Length", 0))

        # Keep only RefSeq mRNA records for now
        if not accession.startswith(("NM_", "XM_")):
            continue

        if accession.startswith("NM_"):
            record_type = "RefSeq curated mRNA"
        else:
            record_type = "RefSeq predicted mRNA"

        candidates.append({
            "gene": gene,
            "species": species,
            "accession": accession,
            "title": title,
            "mRNA_length": length,
            "record_type": record_type,
            "notes": "",
        })

    time.sleep(0.4)
    return candidates

def fetch_genbank(accession: str):
    """Fetch one NCBI nucleotide record in GenBank format."""
    print(f"Fetching GenBank record: {accession}")

    with Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text",
    ) as handle:
        record = SeqIO.read(handle, "genbank")

    time.sleep(0.4)
    return record

def extract_cds(record, expected_gene: str) -> list[dict]:
    """Extract CDS features from a GenBank record that match the expected gene."""
    cds_list = []

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene_name = feature.qualifiers.get("gene", [""])[0]
        product = feature.qualifiers.get("product", [""])[0]
        protein_id = feature.qualifiers.get("protein_id", [""])[0]
        db_xrefs = feature.qualifiers.get("db_xref", [])

        gene_matches = gene_name.upper() == expected_gene.upper()

        cds_seq = feature.extract(record.seq)
        location = format_cds_location(feature)

        if gene_matches:
            cds_list.append({
                "accession": record.id,
                "gene": gene_name,
                "product": product,
                "protein_id": protein_id,
                "db_xrefs": db_xrefs,
                "location": location,
                "cds_length": len(cds_seq),
                "sequence": str(cds_seq),
            })

    return cds_list

def build_ncbi_style_header(common_name: str, cds: dict, variant_note: str = "") -> str:
    """Build a FASTA header similar to NCBI's downloaded CDS FASTA header."""
    clean_common_name = clean_name(common_name)

    # Add transcript variant to the common name if present.
    # Example: Mouse_Transcript_Variant_2
    if variant_note:
        clean_variant = clean_name(variant_note)
        display_name = f"{clean_common_name}_{clean_variant}"
    else:
        display_name = clean_common_name

    accession = cds["accession"]
    protein_id = cds["protein_id"]
    gene = cds["gene"]
    product = cds["product"]
    location = cds["location"]

    # db_xref part
    db_xref_text = ""
    if cds["db_xrefs"]:
        db_xref_text = f" [db_xref={','.join(cds['db_xrefs'])}]"

    header = (
        f">{display_name} "
        f"lcl|{accession}_cds_{protein_id}_1 "
        f"[gene={gene}]"
        f"{db_xref_text} "
        f"[protein={product}] "
        f"[protein_id={protein_id}] "
        f"[location={location}] "
        f"[gbkey=CDS]"
    )

    return header

def clean_name(name: str) -> str:
    """Make a name safe for FASTA headers."""
    return (
        str(name)
        .strip()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
        .replace(",", "")
    )


def get_variant_note(title: str) -> str:
    """Extract transcript variant information from the NCBI title if present."""
    title_lower = title.lower()

    if "transcript variant" not in title_lower:
        return ""

    # Example:
    # "Mus musculus calmodulin 1 (Calm1), transcript variant 2, mRNA"
    parts = title.split(",")

    for part in parts:
        if "transcript variant" in part.lower():
            return part.strip()

    return "transcript variant mentioned"

def format_cds_location(feature) -> str:
    """Convert a CDS feature location to a 1-based NCBI-like location string."""
    location = feature.location

    # Simple location, e.g. 195..644
    if not hasattr(location, "parts"):
        start = int(location.start) + 1
        end = int(location.end)
        return f"{start}..{end}"

    # Compound location, e.g. join(100..200,300..400)
    parts = []
    for part in location.parts:
        start = int(part.start) + 1
        end = int(part.end)
        parts.append(f"{start}..{end}")

    return "join(" + ",".join(parts) + ")"

def write_combined_sequence_txt(fasta_records: list[dict], output_path: Path) -> None:
    """Write all selected CDS sequences into one FASTA file."""
    with open(output_path, "w") as fasta_file:
        for record in fasta_records:
            fasta_file.write(record["header"] + "\n")

            sequence = record["sequence"]
            for i in range(0, len(sequence), 70):
                fasta_file.write(sequence[i:i + 70] + "\n")

            fasta_file.write("\n")

def main():
    target_gene = input("Enter target gene: ").strip().upper()

    input_dir = Path("input")
    base_output_dir = Path("output")
    gene_output_dir = base_output_dir / target_gene
    gene_output_dir.mkdir(parents=True, exist_ok=True)

    species_file = input_dir / "Additional_sampling_species.xlsx"
    species_records = read_species_list(species_file)

    tracking_rows = []
    fasta_records = []

    for species_record in species_records:
        taxonomic_group = species_record["taxonomic_group"]
        common_name = species_record["common_name"]
        scientific_name = species_record["scientific_name"]

        print(f"\nSearching {target_gene} in {common_name} ({scientific_name})")

        candidates = search_ncbi_candidates(
            gene=target_gene,
            species=scientific_name,
            max_results=10,
        )

        if not candidates:
            tracking_rows.append({
                "Taxonomic group": taxonomic_group,
                "Species (Common Name)": common_name,
                "Species (Scientific Name)": scientific_name,
                "Accession Numbers": "-",
                "Sequence Type": "-",
                "CDS Length (bp)": "-",
                "Notes": "no NM_/XM_ RefSeq mRNA candidates found",
                "Title": "-",
            })
            continue

        for candidate in candidates:
            accession = candidate["accession"]
            title = candidate["title"]
            sequence_type = candidate["record_type"] + " CDS"

            try:
                record = fetch_genbank(accession)
                cds_entries = extract_cds(record, target_gene)

            except Exception as error:
                tracking_rows.append({
                    "Taxonomic group": taxonomic_group,
                    "Species (Common Name)": common_name,
                    "Species (Scientific Name)": scientific_name,
                    "Accession Numbers": accession,
                    "Sequence Type": sequence_type,
                    "CDS Length (bp)": "-",
                    "Notes": f"error fetching/extracting CDS: {error}",
                    "Title": title,
                })
                continue

            if not cds_entries:
                tracking_rows.append({
                    "Taxonomic group": taxonomic_group,
                    "Species (Common Name)": common_name,
                    "Species (Scientific Name)": scientific_name,
                    "Accession Numbers": accession,
                    "Sequence Type": sequence_type,
                    "CDS Length (bp)": "-",
                    "Notes": "no matching CDS feature found",
                    "Title": title,
                })
                continue

            for cds in cds_entries:
                notes = get_variant_note(title)

                fasta_header = build_ncbi_style_header(
                    common_name=common_name,
                    cds=cds,
                    variant_note=notes,
                )

                fasta_records.append({
                    "header": fasta_header,
                    "sequence": cds["sequence"],
                })

                tracking_rows.append({
                    "Taxonomic group": taxonomic_group,
                    "Species (Common Name)": common_name,
                    "Species (Scientific Name)": scientific_name,
                    "Accession Numbers": accession,
                    "Sequence Type": sequence_type,
                    "CDS Length (bp)": cds["cds_length"],
                    "Notes": notes,
                    "Title": title,
                })

    # Keep column order clean.
    columns = [
        "Taxonomic group",
        "Species (Common Name)",
        "Species (Scientific Name)",
        "Accession Numbers",
        "Sequence Type",
        "CDS Length (bp)",
        "Notes",
        "Title",
    ]

    tracking_df = pd.DataFrame(tracking_rows, columns=columns)

    tracking_output_path = gene_output_dir / f"{target_gene}_tracking_sheet.csv"
    sequence_output_path = gene_output_dir / f"{target_gene}.txt"

    tracking_df.to_csv(tracking_output_path, index=False)
    write_combined_sequence_txt(fasta_records, sequence_output_path)

    print("\nDone.")
    print(f"Tracking sheet saved to: {tracking_output_path}")
    print(f"Combined sequence text file saved to: {sequence_output_path}")

if __name__ == "__main__":
    main()
