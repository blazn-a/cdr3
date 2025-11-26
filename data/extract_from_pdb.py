# -*- coding: ascii -*-
import os
import re
import pandas as pd
import concurrent.futures
from functools import partial

# ---------------- CONFIG ----------------
IMGTPDB_DIR = "data/imgt_pdbs"
OUT_CSV = "cdr3_antigen_dataset.csv"
CDR3_START = 105
CDR3_END = 117

# ---------------- AA MAP ----------------
AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O'
}

# Pre-compile Regex
PAIRED_HL_PATTERN = re.compile(
    r"PAIRED_HL\s+HCHAIN=([A-Za-z0-9])\s+LCHAIN=([A-Za-z0-9])\s+AGCHAIN=([\w,]+)"
)

def three_to_one(res):
    return AA_MAP.get(res, 'X')

# ---------------- CORE PROCESSOR ----------------
def process_single_pdb(filepath):
    """
    Parses a single PDB file in one pass.
    Returns a list of row dictionaries.
    """
    pdb_id = os.path.splitext(os.path.basename(filepath))[0].upper()
    
    try:
        with open(filepath, "r") as f:
            # Read all lines once
            lines = f.readlines()
    except Exception as e:
        return []

    # Data containers
    pairs = []
    seqres_data = {} # chain -> list of chars
    
    # ATOM data containers
    # chain -> list of chars
    atom_full_seq = {} 
    # chain -> list of chars (filtered by IMGT 105-117 range)
    atom_cdr3_parts = {} 
    
    # State trackers for parsing ATOMs efficiently
    # chain -> last_residue_id_string (e.g. "112A") to handle deduplication
    last_seen_res = {} 

    for line in lines:
        # 1. Parse REMARK 5 PAIRED_HL
        if line.startswith("REMARK   5") and "PAIRED_HL" in line:
            m = PAIRED_HL_PATTERN.search(line)
            if m:
                pairs.append({
                    "H": m.group(1),
                    "L": m.group(2),
                    "AG": m.group(3).split(",")
                })
            continue

        # 2. Parse SEQRES
        if line.startswith("SEQRES"):
            if len(line) < 19: continue
            chain = line[11] # Column 12 (index 11) is chain ID
            residues = line[19:].split()
            if chain not in seqres_data:
                seqres_data[chain] = []
            for r in residues:
                seqres_data[chain].append(AA_MAP.get(r, 'X'))
            continue

        # 3. Parse ATOM
        if line.startswith("ATOM"):
            # Fixed width slicing is faster than strip() or split()
            chain = line[21]
            res_name = line[17:20]
            
            # IMGT numbering: Cols 23-26 (index 22:26) is number, Col 27 (index 26) is insertion
            # Combined ID for deduplication (e.g., " 112 A")
            res_id_full = line[22:27] 

            # Quick check: Only process if we moved to a new residue
            if last_seen_res.get(chain) == res_id_full:
                continue
            
            # Update state
            last_seen_res[chain] = res_id_full

            # Convert AA
            aa = AA_MAP.get(res_name, 'X') # res_name is usually upper in PDB
            
            # Append to full ATOM sequence
            if chain not in atom_full_seq:
                atom_full_seq[chain] = []
                atom_cdr3_parts[chain] = []
            
            atom_full_seq[chain].append(aa)

            # Check CDR3 Range (105-117)
            # Parse the integer part. PDB format guarantees this is numeric or space.
            try:
                res_num = int(line[22:26])
                if CDR3_START <= res_num <= CDR3_END:
                    atom_cdr3_parts[chain].append(aa)
            except ValueError:
                # Handle edge cases where res number is empty/malformed
                pass

    if not pairs:
        return []

    # Build Rows
    rows = []
    
    # Join SEQRES lists into strings once
    seqres_strings = {k: "".join(v) for k, v in seqres_data.items()}
    atom_strings = {k: "".join(v) for k, v in atom_full_seq.items()}
    cdr3_strings = {k: "".join(v) for k, v in atom_cdr3_parts.items()}

    for pair in pairs:
        H = pair["H"]
        L = pair["L"]
        AG_list = pair["AG"]

        # Require heavy and light CDR3s to exist
        if H not in cdr3_strings or L not in cdr3_strings:
            continue

        heavy_cdr3 = cdr3_strings[H]
        light_cdr3 = cdr3_strings[L]

        for ag in AG_list:
            # Priority: SEQRES -> ATOM -> Skip
            if ag in seqres_strings:
                antigen_seq = seqres_strings[ag]
            elif ag in atom_strings:
                antigen_seq = atom_strings[ag]
            else:
                continue

            rows.append({
                "pdb_id": pdb_id,
                "heavy_cdr3": heavy_cdr3,
                "light_cdr3": light_cdr3,
                "antigen_chain": ag,
                "antigen_seq": antigen_seq
            })

    return rows

# ---------------- RUN MANAGER ----------------
def run_all(imgt_dir=IMGTPDB_DIR, out_csv=OUT_CSV):
    if not os.path.isdir(imgt_dir):
        raise RuntimeError(f"IMGTPDB_DIR not found: {imgt_dir}")

    # Gather all files
    files = [
        os.path.join(imgt_dir, f) 
        for f in os.listdir(imgt_dir) 
        if f.lower().endswith(".pdb")
    ]
    
    if not files:
        print("No PDB files found.")
        return

    all_rows = []
    
    # Parallel Execution
    # Adjust max_workers based on your CPU; None defaults to number of processors
    print(f"Processing {len(files)} files...")
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Map returns results in order
        results = executor.map(process_single_pdb, files)
        
        for res in results:
            if res:
                all_rows.extend(res)

    df = pd.DataFrame(all_rows)
    df.to_csv(out_csv, index=False)
    print(f"Wrote {len(df)} rows to {out_csv}")

# ---------------- MAIN ----------------
if __name__ == "__main__":
    run_all()