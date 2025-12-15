import os
import uuid
import time
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from Tic_sys import FOLDER_PATH

from Bio import Entrez
import xml.etree.ElementTree as ET

# ✅ Foldseek API endpoint
FOLDSEEK_URL = "https://search.foldseek.com/api/ticket"

# ---------- STEP 1: Scan folder for PDB/CIF files ----------
def get_pdb_files(folder_path):
    """
    Scans the given folder and returns a list of .pdb or .cif files.
    """
    folder = Path(folder_path)
    if not folder.exists():
        raise FileNotFoundError(f"❌ Folder '{folder_path}' does not exist.")

    pdb_files = [
        str(file.resolve())
        for file in folder.iterdir()
        if file.is_file() and file.suffix.lower() in [".pdb", ".cif"]
    ]

    print(f" Found {len(pdb_files)} structure files in '{folder_path}'")
    for f in pdb_files:
        print(" -", os.path.basename(f))

    return pdb_files


# ---------- STEP 2: Submit one structure to Foldseek ----------
def submit_to_foldseek(pdb_path):
    """
    Submits a single PDB/CIF file to Foldseek and returns a ticket dict.
    """
    try:
        with open(pdb_path, "rb") as f:
            files = {"q": f}
            data = [
                ("mode", "3diaa"),
                ("database[]", "afdb-proteome"),  # ✅ Only search AFDB50
            ]
            r = requests.post(FOLDSEEK_URL, files=files, data=data)

        print("🔎 Response status:", r.status_code)
        print("🔎 Response body:", r.text)  # <-- add this line
        job_info = r.json()

        ticket = {
            "ticket_id": str(uuid.uuid4()),  # our internal ID
            "file": pdb_path,  # input file
            "job_id": job_info.get("id"),  # Foldseek's job ID
            "status": job_info.get("status", "unknown"),
            "created_at": time.time(),
        }

        print(f" Submitted {os.path.basename(pdb_path)} → job_id: {ticket['job_id']}")
        return ticket

    except Exception as e:
        print(f"Failed to submit {pdb_path}: {e}")
        return {
            "ticket_id": str(uuid.uuid4()),
            "file": pdb_path,
            "job_id": None,
            "status": "error",
            "error": str(e),
        }


# ---------- STEP 3: Submit all files concurrently ----------
def submit_all(folder_path, max_workers=5):
    """
    Scans a folder and submits all found PDB/CIF files to Foldseek concurrently.
    Returns a list of ticket dictionaries.
    """
    pdb_files = get_pdb_files(folder_path)
    tickets = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(submit_to_foldseek, pdb): pdb for pdb in pdb_files}
        for future in as_completed(futures):
            ticket = future.result()
            tickets.append(ticket)

    print("\n Submission complete. Generated tickets:")
    for t in tickets:
        print(
            f" - {os.path.basename(t['file'])}: "
            f"job_id={t['job_id']} status={t['status']}"
        )

    return tickets


# ---------- STEP 4: Fetch results ----------
def fetch_html_result(job_id, output_folder="results_html"):
    print(f"🔗 View results in browser:")
    print(f"https://search.foldseek.com/result/{ticket['job_id']}/0")


"""
Goal: retrieve the Uniprot IDs (field 'target') for each pbd file (only for hits from the AFDB-proteome database)
and sort them according to e-value (field 'eval')
"""

def save_dict_to_txt(dictionary: dict, filename: str):
    """
    Save a Python dictionary to a text file in a reusable format.
    """
    with open(filename, "w", encoding="utf-8") as f:
        f.write(repr(dictionary))
    print(f"Saved {filename}")

"""
First, a function to get the complete Foldseek output from the ticket ID of completed jobs:
"""

def get_full_Foldseek_results(folder_path, max_workers=5):
    tickets = submit_all(folder_path,
                         max_workers=max_workers)  # list of dictionaries with ticket_id, file path, job_id, status, error for each submission

    all_results = {}  # key: pdb file name, value: Foldseek results

    for t in tickets:
        # Select only completed tickets ("status": "COMPLETE")
        if t['status'] == "COMPLETE":
            job_id = t.get("job_id")
            # Use the function 'os.path.basename' to extract name of pdb file from path
            pdb_file_name = os.path.basename(t.get("file"))
            entry = 0  # 1 pdf file at a time
            url = f"https://search.foldseek.com/api/result/{job_id}/{entry}"
            results = requests.get(url)
            all_results[pdb_file_name] = results.json()

    return all_results


"""
A function to extract the Uniprot IDs of the top 10 hits from the AFDB-proteome db with an e-value < 1e-5 and 
sort them from lowest to highest (most significant hit on top)
"""

def get_top10_Uniprot_ids(all_results, evalue):
    Uniprot_ids = {}  # key: pdb file, value: list of Uniprot IDs of top 10 most significant hits

    for pdb_file, foldseek_result in all_results.items():

        significant_hits = []

        for alignment_list in foldseek_result["results"][0]["alignments"]:
            """Check example output: foldseek_result["results"] is a list of dictionaries (one for each database).
            We select index 0, since we are only searching against AFDB-proteome db

            "alignments" is a list of dictionaries with keys (like query, target, eval etc.), 
            each representing a separate alignment/hit/homolog (1000 in total)
            """
            for alignment in alignment_list:
                if alignment["eval"] < evalue:
                    significant_hits.append(alignment)

        # Smallest e-value on top, select top 10 hits
        significant_hits.sort(key=lambda hit: hit["eval"])
        top_10_hits = significant_hits[:]  # was 10

        # We only need the Uniprot IDs ("target"), stored as a list under pdb file key
        pdb_uniprot_ids = [hit["target"].split("-")[1] for hit in top_10_hits]  # list comprehension
        Uniprot_ids[pdb_file] = pdb_uniprot_ids

    return Uniprot_ids


def foldseek_main(input, evalue):
    evalue = float(evalue)
    all_results = get_full_Foldseek_results(folder_path=input, max_workers=5)
    top10_ids = get_top10_Uniprot_ids(all_results, evalue)

    # Print top 10 hits (Uniprot IDs) for each PDB file
    # for pdb_file, ids in top10_ids.items():
    # print(f"\nTop {len(ids)} hits for {pdb_file}:")
    # for uniprot_id in ids:
    # print(f" - {uniprot_id}")
    return top10_ids


# Connection function for ID-mapping
def map_uniprot_to_ncbi(uniprot_ids, mapping_file):
    uniprot_to_name = {uid: name for name, ids in uniprot_ids.items() for uid in ids}
    # Prepare a dictionary for the output with NCBI IDs
    NCBI_ids = {file_name: [] for file_name in uniprot_ids.keys()}
    with open(mapping_file, 'r') as f:
        for line in f:
            uniprot_id, db, mapped_id = line.strip().split('\t')

            # We only keep NCBI Gene IDs
            if uniprot_id in uniprot_to_name:
                name = uniprot_to_name[uniprot_id]
                NCBI_ids[name].append(mapped_id)
    return NCBI_ids

def ID_mapping_main(mapping, top10_ids):
    parsed_map = mapping
    ncbi_ids = map_uniprot_to_ncbi(top10_ids, parsed_map)
    return ncbi_ids


def get_genomic_coordinates(gene_id: int) -> tuple:
    Entrez.email = 'A.N.Other@example.com'

    # Fetch gene record (XML format)
    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
    record = handle.read()

    # Parse the XML data
    root = ET.fromstring(record)

    found_current = False
    start_values = []
    end_values = []

    # Select the most recent annotation release (status "current") for genes with multiple annotation releases on NCBI.
    commentaries = root.findall(".//Gene-commentary")
    for comment in commentaries:
        text_tag = comment.find(".//Gene-commentary_text")
        if text_tag is not None and text_tag.text == "current":  # must be an exact match; e.g. for gene ID 873 'current' is part of the string 'concurrent' in the text under the first <Gene-commentary_text> tag
            found_current = True
            # print(f"Found first 'current' annotation release for gene ID {gene_id}")
            seqs_tag = comment.find(".//Gene-commentary_seqs")
            if seqs_tag is not None:
                seq_intervals = seqs_tag.findall(".//Seq-interval")  # in case there are multiple <Seq-interval> tags
                for seq_interval in seq_intervals:
                    start = seq_interval.find('./Seq-interval_from')
                    end = seq_interval.find('./Seq-interval_to')
                    if start is not None and end is not None:
                        start_values.append(int(start.text))
                        end_values.append(int(end.text))

            break  # Exit loop once first "current" annotation release is found;
            # for some genes, e.g. gene ID 26063, 'current' appears twice as annotation release version on NCBI

    if not found_current:
        # print(f"Only one annotation release found for gene ID {gene_id}.")
        # Check if <Seq-interval_from/to> tags exist somewhere in entire gene record
        start = root.findall(".//Seq-interval_from")
        end = root.findall(".//Seq-interval_to")
        if start:
            for i in start:
                start_values.append(int(i.text))
        if end:
            for i in end:
                end_values.append(int(i.text))
        # if not start_values or not end_values:
        # print(f"Missing coordinates for gene ID {gene_id}")
    if start_values == []:
        # print(f"Genomic coordinates for gene ID don't exist")
        return (0, 0)
    else:
        lowest = min(start_values)
        highest = max(end_values)
        # print(f"Genomic coordinates for gene ID {gene_id}: ({lowest}, {highest})")
    return (lowest, highest)


def main_coordinates(ncbi_ids):
    coordinates = {}
    count = 0
    for name in ncbi_ids:

        ID = {}
        for NCBI_ID in ncbi_ids[name]:
            """
            #add this count-stopper part when testing
            if count == 5:
                return coordinates
            """
            ID[NCBI_ID] = get_genomic_coordinates(NCBI_ID)
            count += 1
        coordinates[name] = ID
    return coordinates


"""
Step 1a. Retrieve scaffold IDs for each gene ID using efetch. 
"""

def get_scaffold_ID(gene_ID: int) -> str:  # takes about 1 sec per gene ID
    Entrez.email = 'A.N.Other@example.com'
    handle = Entrez.efetch(db="gene", id=gene_ID, retmode="xml")
    record = handle.read()
    # Parse the XML data (convert string to ElementTree object)
    root = ET.fromstring(record)
    # Check if <Gene-commentary_accession> tag exists
    scaffold_IDs = root.findall(".//Gene-commentary_accession")
    NC = []
    NZ = []

    if scaffold_IDs:
        count = 0
        for scaffold_ID in scaffold_IDs:
            """
            #remove this count-stopper part when testing
            if count == 5:
                break
            """
            """
            print(scaffold_ID.text) gives a list of scaffold_IDs
            e.g. for gene ID '22220378': NC_001659, NP_042842, NC_001659, NP_042842, KM102979, AJZ77157...
            """
            if scaffold_ID.text.startswith("NC"):
                NC.append(scaffold_ID.text)  # e.g. NC_001659 appears multiple times
            elif scaffold_ID.text.startswith("NZ"):
                NZ.append(scaffold_ID.text)
            count += 1

    if NC:
        first_NC_scaffold_ID = NC[0]
        # print(f"Scaffold ID for gene {gene_ID}: {first_NC_scaffold_ID}")
        return first_NC_scaffold_ID
    elif NZ:
        first_NZ_scaffold_ID = NZ[0]
        # print(f"Scaffold ID for gene {gene_ID}: {first_NZ_scaffold_ID}")
        return first_NZ_scaffold_ID
    elif scaffold_IDs:
        first_other_scaffold_ID = scaffold_IDs[0].text
        # print(f"Scaffold ID for gene {gene_ID}: {first_other_scaffold_ID}")
        return first_other_scaffold_ID
    # else:
    # print(f"No scaffold IDs found for gene {gene_ID}")


"""
Step 1b. Store results in scaffold_dictionary under gene ID as key.
Use ThreadPoolExecutor to fetch scaffold IDs for multiple gene IDs in parallel.
"""

def create_scaffold_dict(gene_IDS: list, max_workers=2) -> dict:
    # 5 or 3 workers gives HTTP Error 429: Too Many Requests
    # too many requests to the NCBI Entrez API
    scaffold_dict = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:  
        futures = {}
        count = 0
        for gene_ID in gene_IDS:
            """
            #remove this count-stopper part when testing
            if count == 5:
                break
            """
            future = executor.submit(get_scaffold_ID, gene_ID)
            # submit 'get_scaffold_ID' function with argument gene_ID to executor; returns future object with result
            futures[future] = gene_ID  # store future object in dictionary and link to corresponding gene_ID
            count += 1

        count = 0
        for future in as_completed(futures):
            """
            #remove this count-stopper part when not testing
            if count == 5:
                return scaffold_dict
            """
            gene_ID = futures[future]
            scaffold_ID = future.result()  # result of get_scaffold_ID
            scaffold_dict[gene_ID] = scaffold_ID
            count += 1

    return scaffold_dict


"""
Step 2. Link the scaffold ID to the corresponding assembly record ID with elink.
"""

def get_assembly_record_ID(scaffold_ID: str) -> int:
    Entrez.email = 'A.N.Other@example.com'

    handle = Entrez.elink(dbfrom="nucleotide", id=scaffold_ID, db="assembly")
    record = handle.read()
    # print(record) #assembly record ID should be somewhere under <LinkSetDb>\n...<Link> <Id>

    root = ET.fromstring(record)
    LinkSetDb = root.find(".//LinkSetDb")

    assembly_record_IDs = []

    if LinkSetDb is not None:
        for link in LinkSetDb.findall(".//Link"):
            Id = link.find("Id")
            if Id is not None:  # sometimes more than 1 Id: e.g. for NC_003070 Assembly record IDs: '1733481', '237408', '3248'
                assembly_record_ID = Id.text
                assembly_record_IDs.append(assembly_record_ID)
            # else:
            # print(f"LinkSetDb for scaffold ID {scaffold_ID} has no corresponding Id")
    # else:
    # print(f"Scaffold {scaffold_ID} is not linked to any assembly record.")
    if not assembly_record_IDs:  # <---- IMPORTANT FIX
        return None
    # print(f"Assembly record ID: {assembly_record_IDs} for scaffold ID {scaffold_ID}")
    return assembly_record_IDs[0]  # first value if multiple Ids


"""
Step 3a. Get assembly ID from assembly record ID with esummary.
"""

def get_assembly_ID(assembly_record_ID: int) -> str:
    Entrez.email = 'A.N.Other@example.com'

    handle = Entrez.esummary(db="assembly", id=assembly_record_ID, report="full")
    record = Entrez.read(handle)

    assembly_accession = record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
    if assembly_accession:
        return assembly_accession.split('.')[0]  # only base for easier comparison during cluster mapping
    # else:
    # print(f"No assembly ID found for assembly record {assembly_record_ID}.")


"""
Step 3b. Store results in assembly_dictionary under gene ID as key.
Use ThreadPoolExecutor for parallelisation.
"""

def create_assembly_dictionary(scaffold_dictionary: dict, max_workers=2) -> dict:
    # Helper function to link scaffold ID to assembly ID
    def link_scaffold_to_assembly_ID(scaffold_ID: str) -> int:
        if scaffold_ID:
            assembly_record_ID = get_assembly_record_ID(scaffold_ID)
            if assembly_record_ID:
                return get_assembly_ID(assembly_record_ID)
            # else:
            # print(f"No assembly record found.")
        # else:
        # print("No scaffold ID found.")

    unlinked_count = 0
    total_scaffolds = 0
    assembly_dictionary = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        count = 0
        for gene_ID, scaffold_ID in scaffold_dictionary.items():
            """
            #remove this count-stopper part when not testing
            if count == 5:
                break
            """
            if scaffold_ID:
                future = executor.submit(link_scaffold_to_assembly_ID, scaffold_ID)
                futures[future] = gene_ID
                total_scaffolds += 1
            count += 1

        count = 0
        for future in as_completed(futures):
            """
            #remove this count-stopper part when not testing
            if count == 5:
                return assembly_dictionary
            """
            gene_ID = futures[future]
            assembly_ID = future.result()
            if assembly_ID:
                assembly_dictionary[gene_ID] = assembly_ID
                # print(f" Genomic assembly ID for Gene ID {gene_ID} is {assembly_ID}")
            else:
                unlinked_count += 1
            count += 1

    print(f"{unlinked_count} out of {total_scaffolds} scaffolds are not linked to an assembly record")
    return assembly_dictionary


def main_scaffold_assembly_ID(ncbi_ids):
    final_scaffold_dict = {}  # dictionary with key: pdb file, value: corresponding scaffold dictionary (maps gene ID to scaffold ID)
    final_assembly_dict = {}  # dictionary with key: pdb file, value: corresponding assembly dictionary (maps gene ID to assembly ID)
    for pdb_file in ncbi_ids:
        gene_ids = ncbi_ids[pdb_file]
        final_scaffold_dict[pdb_file] = create_scaffold_dict(gene_ids)
        print(f"Scaffold IDs for {pdb_file}:")
        print(final_scaffold_dict[pdb_file])

        scaffold_ids = final_scaffold_dict[pdb_file]
        final_assembly_dict[pdb_file] = create_assembly_dictionary(scaffold_ids)
        print(f"Genomic assembly IDs for {pdb_file}:")
        print(final_assembly_dict[pdb_file])

    return final_scaffold_dict, final_assembly_dict


def run(input, mapping, evalue):
    # foldseek_main() is only giving top10 results, for testing, this can be changed later easily
    top_ids = foldseek_main(input, evalue)
    print("✔ Search for homologs using Foldseek complete!")
    # The ID mapping part
    ncbi_ids = ID_mapping_main(mapping, top_ids)
    print("✔ ID mapping complete!")
    # print(ncbi_ids)
    # for pdb_file, ids in ncbi_ids.items():
    # print(f"\nTop {len(ids)} hits for {pdb_file}:")
    # for Gene_id in ids:
    # print(f" - {Gene_id}")
    coordinates = main_coordinates(ncbi_ids)
    print("✔ coordinates retrieved!")
    save_dict_to_txt(coordinates, "Coordinates_dictionary.txt")

    # scaffold and genomic assembly ID
    final_scaffold_dict, final_assembly_dict = main_scaffold_assembly_ID(
        ncbi_ids)  # sometimes interrupted: urllib.error.HTTPError: HTTP Error 429: Too Many Requests
    print("✔ Scaffold and assembly IDs retrieved!")
    save_dict_to_txt(final_scaffold_dict, "Scaffold_dictionary.txt")
    save_dict_to_txt(final_assembly_dict, "Assembly_dictionary.txt")


if __name__ == "__main__":
    run()
