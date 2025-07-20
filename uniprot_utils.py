import requests
import os
import re

def get_uniprot_ids_from_pdb(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    response = requests.get(url)
    if response.status_code != 200:
        return []
    data = response.json()
    if pdb_id.lower() not in data or 'UniProt' not in data[pdb_id.lower()]:
        return []
    return list(data[pdb_id.lower()]['UniProt'].keys())

def fetch_uniprot_entry(uniprot_id, save_dir, cache):
    if uniprot_id in cache:
        return cache[uniprot_id]
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
    response = requests.get(url)
    if response.status_code != 200:
        cache[uniprot_id] = None
        return None
    entry_text = response.text
    filepath = os.path.join(save_dir, f"{uniprot_id}.txt")
    with open(filepath, "w") as f:
        f.write(entry_text)
    cache[uniprot_id] = entry_text
    return entry_text

def is_transmembrane_from_text(entry_text):
    return entry_text and any("FT   TRANSMEM" in line or "Transmembrane" in line for line in entry_text.splitlines())

def is_virus_from_text(entry_text):
    if not entry_text:
        return False
    for line in entry_text.splitlines():
        if line.startswith("OX") and "NCBI_TaxID=10239" in line:
            return True
    for line in entry_text.splitlines():
        if line.startswith("OS") and "virus" in line.lower():
            return True
    return False

def extract_go_terms_from_entry(entry_text):
    go_terms = set()
    if not entry_text:
        return go_terms
    for line in entry_text.splitlines():
        if line.startswith("DR   GO;"):
            match = re.search(r"GO:(\d+);\s+([CFP]):([^;]+)", line)
            if match:
                go_id, aspect, term = match.groups()
                go_terms.add(f"GO:{go_id} [{aspect}:{term.strip()}]")
    return go_terms

def extract_classification_from_entry(entry_text):
    classifications = set()
    if not entry_text:
        return classifications
    for line in entry_text.splitlines():
        if line.startswith("KW   "):
            keywords = line[5:].split(";")
            classifications.update(kw.strip() for kw in keywords if kw.strip())
    return classifications