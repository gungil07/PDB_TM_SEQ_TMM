import requests
from uniprot_utils import get_uniprot_ids_from_pdb, fetch_uniprot_entry, is_transmembrane_from_text, is_virus_from_text, extract_go_terms_from_entry, extract_classification_from_entry
import os
import csv

def fetch_released_pdb_ids(since_date):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_accession_info.initial_release_date",
                "operator": "greater_or_equal",
                "value": since_date
            }
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    response = requests.post(url, json=query)
    if response.status_code != 200:
        print(f"Failed to fetch from RCSB ({response.status_code})")
        return []
    return [item["identifier"] for item in response.json().get("result_set", [])]

def fetch_metadata(pdb_id, save_dir, cache):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to fetch metadata for {pdb_id}")
        return None
    data = response.json()
    title = data.get("struct", {}).get("title", "")
    pubmed_id = data.get("rcsb_primary_citation", {}).get("pdbx_database_id_pub_med", "")
    citations = data.get("citation", [])
    journal_abbrev = year = volume = page_first = page_last = ""
    for cit in citations:
        if cit.get("rcsb_citation", {}).get("primary", False):
            journal_abbrev = cit.get("journal_abbrev", "")
            year = cit.get("year", "")
            volume = cit.get("volume", "")
            page_first = cit.get("page_first", "")
            page_last = cit.get("page_last", "")
            break
    if not journal_abbrev and citations:
        cit = citations[0]
        journal_abbrev = cit.get("journal_abbrev", "")
        year = cit.get("year", "")
        volume = cit.get("volume", "")
        page_first = cit.get("page_first", "")
        page_last = cit.get("page_last", "")
    method = data.get("exptl", [{}])[0].get("method", "")
    resolution = data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]

    uniprot_ids = get_uniprot_ids_from_pdb(pdb_id)
    transmembrane = is_viral = "No"
    all_go_terms = set()
    all_classifications = set()

    for uid in uniprot_ids:
        entry_text = fetch_uniprot_entry(uid, save_dir, cache)
        if is_transmembrane_from_text(entry_text):
            transmembrane = "Yes"
        if is_virus_from_text(entry_text):
            is_viral = "Yes"
        all_go_terms.update(extract_go_terms_from_entry(entry_text))
        all_classifications.update(extract_classification_from_entry(entry_text))
        if transmembrane == "Yes" and is_viral == "Yes":
            break

    return {
        "PDB_ID": pdb_id,
        "Title": title,
        "PubMed_ID": pubmed_id,
        "Journal_Abbrev": journal_abbrev,
        "Year": year,
        "Volume": volume,
        "Page_First": page_first,
        "Page_Last": page_last,
        "Resolution (Å)": resolution,
        "Experimental Method": method,
        "UniProt_IDs": ", ".join(uniprot_ids),
        "Transmembrane": transmembrane,
        "Virus": is_viral,
        "Classification": ", ".join(sorted(all_classifications)),
        "GO_Terms": ", ".join(sorted(all_go_terms))
    }

def save_to_csv(pdb_data_list, output_file):
    fieldnames = [
        "PDB_ID", "Title", "PubMed_ID", "Journal_Abbrev", "Year", "Volume",
        "Page_First", "Page_Last", "Resolution (Å)", "Experimental Method",
        "UniProt_IDs", "Transmembrane", "Virus", "Classification", "GO_Terms"
    ]
    with open(output_file, "w", newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in pdb_data_list:
            writer.writerow(row)
