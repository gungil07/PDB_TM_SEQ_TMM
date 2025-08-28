from pdb_utils import fetch_released_pdb_ids, fetch_metadata, save_to_csv
from fasta_fetch import fetch_from_csv
import os
from datetime import datetime


def main():
    since_date = input("Enter a release date (YYYY-MM-DD): ").strip()
    try:
        datetime.strptime(since_date, "%Y-%m-%d")
    except ValueError:
        print("Invalid date format. Use YYYY-MM-DD.")
        return

    today_str = datetime.now().strftime('%Y-%m-%d')
    output_dir = f"uniprot_entries_{since_date.replace('-', '')}"
    os.makedirs(output_dir, exist_ok=True)

    print("\n[1] Fetching released PDB IDs...")
    pdb_ids = fetch_released_pdb_ids(since_date)
    if not pdb_ids:
        return

    print("\n[2] Fetching metadata for each PDB ID...")
    all_data = []
    filtered_data = []
    uniprot_cache = {}
    for i, pdb_id in enumerate(pdb_ids, 1):
        print(f"Fetching metadata for {pdb_id} ({i}/{len(pdb_ids)})...")
        metadata = fetch_metadata(pdb_id, output_dir, uniprot_cache)
        if metadata:
            all_data.append(metadata)
            if (
                metadata.get("Virus") == "No"
                and metadata.get("Transmembrane") == "Yes"
                and metadata.get("Journal_Abbrev") not in ("", None, "To Be Published", "Biorxiv")
            ):
                filtered_data.append(metadata)

    # Save metadata CSVs
    original_file = f"pdb_metadata_all_since_{since_date}_saved_{today_str}.csv"
    filtered_file = f"pdb_metadata_filtered_since_{since_date}_saved_{today_str}.csv"

    save_to_csv(all_data, original_file)
    save_to_csv(filtered_data, filtered_file)

    print("\n[3] Fetching FASTA sequences for filtered PDBs...")
    fetch_from_csv(filtered_file)


if __name__ == "__main__":
    main()
