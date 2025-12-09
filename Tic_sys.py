import os
import uuid
import time
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# ✅ Foldseek API endpoint
FOLDSEEK_URL = "https://search.foldseek.com/api/ticket"

#  Change your folder path here!
FOLDER_PATH = r"C:\Users\aniru\Documents\KUL\Third Sem\Integrated Project\PDB_files"

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
        print("🔎 Response body:", r.text)   # <-- add this line
        job_info = r.json()

        ticket = {
            "ticket_id": str(uuid.uuid4()),         # our internal ID
            "file": pdb_path,                       # input file
            "job_id": job_info.get("id"),        # Foldseek's job ID
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



# ---------- MAIN ----------
if __name__ == "__main__":
    tickets = submit_all(FOLDER_PATH, max_workers=5)

    print("\n📥 Fetching results for completed jobs...\n")
    for ticket in tickets:
        if ticket["job_id"]:
            try:
                fetch_html_result(ticket["job_id"])
            except Exception as e:
                print(f"⚠️ Could not fetch results for {ticket['job_id']}: {e}")
        else:
            print(f"⚠️ No job_id for {ticket['file']} — skipping result fetch.")
