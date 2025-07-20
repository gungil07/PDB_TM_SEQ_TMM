import requests
import os
import csv
import re
from datetime import datetime

def fetch_pdb_fasta_by_chain(pdb_id):
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(url)
    return response.text if response.status_code == 200 and response.text.strip() else None

def parse_fasta_entries(fasta_text):
    entries = []
    current_header = ""
    current_seq = []
    for line in fasta_text.strip().splitlines():
        if line.startswith(">"):
            if current_header:
                entries.append((current_header, ''.join(current_seq)))
            current_header = line.strip()
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_header:
        entries.append((current_header, ''.join(current_seq)))
    return entries

def extract_info_from_header(header):
    match = re.match(r'^>(\w{4})_([A-Za-z0-9])\s+.*?Chain\s+[A-Za-z0-9],\s*(.*)', header)
    return match.groups() if match else ("UNKNOWN", "?", header[1:])

def fetch_from_csv(csv_path, save_dir="fasta_output"):
    os.makedirs(save_dir, exist_ok=True)
    date_str = datetime.now().strftime("%Y%m%d")
    fasta_path = os.path.join(save_dir, f"pdb_fasta_by_chain_{date_str}.fasta")
    txt_path = os.path.join(save_dir, f"pdb_fasta_by_chain_{date_str}.txt")
    csv_path_out = os.path.join(save_dir, f"pdb_fasta_by_chain_{date_str}.csv")
    failed_log_path = os.path.join(save_dir, f"failed_ids_{date_str}.txt")
    pdb_ids = []
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if "PDB_ID" in reader.fieldnames:
            pdb_ids = [row["PDB_ID"].strip() for row in reader if row["PDB_ID"].strip()]
    failures = []
    all_entries = []
    for pdb_id in pdb_ids:
        fasta_text = fetch_pdb_fasta_by_chain(pdb_id)
        if not fasta_text:
            failures.append(pdb_id)
            continue
        parsed = parse_fasta_entries(fasta_text)
        for header, sequence in parsed:
            pdb, chain, desc = extract_info_from_header(header)
            all_entries.append({
                "PDB_ID": pdb,
                "Chain_ID": chain,
                "Description": desc,
                "Sequence": sequence,
                "Header": header
            })
    with open(fasta_path, "w") as fasta_out, open(txt_path, "w") as txt_out:
        for entry in all_entries:
            fasta_block = f"{entry['Header']}\n{entry['Sequence']}\n"
            fasta_out.write(fasta_block)
            txt_out.write(fasta_block)
    with open(csv_path_out, "w", newline='') as csv_out:
        writer = csv.DictWriter(csv_out, fieldnames=["PDB_ID", "Chain_ID", "Description", "Sequence"])
        writer.writeheader()
        for entry in all_entries:
            writer.writerow({
                "PDB_ID": entry["PDB_ID"],
                "Chain_ID": entry["Chain_ID"],
                "Description": entry["Description"],
                "Sequence": entry["Sequence"]
            })
    if failures:
        with open(failed_log_path, "w") as fail_file:
            fail_file.write("\n".join(failures))