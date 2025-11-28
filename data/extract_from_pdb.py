# -*- coding: ascii -*-
import os
import pandas as pd
import concurrent.futures
from collections import Counter, defaultdict

# ---------------- CONFIGURATION ----------------
# NOTE: Ensure these directories exist and contain the relevant PDB files
IMGTPDB_DIR = "data/imgt_pdbs"       # Source for CDR3s and Headers
ORIGINAL_PDB_DIR = "data/original_pdbs" # Source for Antigen SEQRES

# Output CSV files will be saved in the script's current directory
OUT_CSV = os.path.join(os.getcwd(), "cdr3_antigen_dataset.csv") 
SIMPLIFIED_OUT_CSV = os.path.join(os.getcwd(), "antigen_cdr3_pairs.csv") 

CDR3_START = 105
CDR3_END = 117

# ---------------- AMINO ACID MAP ----------------
AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'MSE': 'M', 'UNK': 'X'
}

# ---------------- CORE PARSING FUNCTIONS ----------------

def parse_header_flexible(line):
    """
    Parses headers to detect:
    1. PAIRED (H + L)
    2. SINGLE (H only)
    3. SINGLE (L only)
    """
    is_single = "SINGLE" in line
    
    # Clean noise
    clean_line = line.replace("REMARK", "").replace("5", "").replace("PAIRED_HL", "").replace("SINGLE", "")
    tokens = clean_line.split()
    
    data = {}
    for token in tokens:
        if "=" in token:
            key, val = token.split("=", 1)
            data[key.strip()] = val.strip()
    
    # Must have an Antigen Chain to be useful
    if "AGCHAIN" not in data:
        return None

    # Handle mixed delimiters in Antigen Chains
    ag_raw = data["AGCHAIN"].replace(";", ",")
    ag_list = [x.strip().upper() for x in ag_raw.split(",") if x.strip()]

    h_chain = data.get("HCHAIN", "").upper() if "HCHAIN" in data else None
    l_chain = data.get("LCHAIN", "").upper() if "LCHAIN" in data else None

    # --- DETERMINE MODE ---
    if is_single:
        if h_chain:
            return {"mode": "SINGLE_H", "H": h_chain, "L": None, "AG": ag_list}
        elif l_chain:
            return {"mode": "SINGLE_L", "H": None, "L": l_chain, "AG": ag_list}
        else:
            return None # SINGLE tag but no chain ID found
    else:
        # PAIRED requires both H and L
        if h_chain and l_chain:
            return {"mode": "PAIRED", "H": h_chain, "L": l_chain, "AG": ag_list}
        
    return None

def get_pdb_id_from_filepath(filepath):
     name = os.path.splitext(os.path.basename(filepath))[0].upper()
     return name.replace('PDB', '')

def parse_seqres_from_original_pdb(pdb_id, original_pdb_dir):
    """Extracts all SEQRES lines from the original PDB file."""
    seqres_data = {}
    
    original_filepaths = [
        os.path.join(original_pdb_dir, f"{pdb_id}.PDB"),
        os.path.join(original_pdb_dir, f"PDB{pdb_id}.ENT"),
        os.path.join(original_pdb_dir, f"{pdb_id}.ENT")
    ]
    
    found_path = next((p for p in original_filepaths if os.path.exists(p)), None)
    
    if not found_path:
        return seqres_data

    try:
        with open(found_path, "r") as f:
            for line in f:
                if line.startswith("SEQRES"):
                    if len(line) < 19: continue
                    chain = line[11].upper()
                    residues = line[19:].split()
                    if chain not in seqres_data:
                        seqres_data[chain] = []
                    for r in residues:
                        seqres_data[chain].append(AA_MAP.get(r, 'X'))
    except Exception:
        pass
    
    return {k: "".join(v) for k, v in seqres_data.items()}

# ---------------- MAIN PROCESS FUNCTION ----------------

def process_single_pdb(imgt_filepath):
    pdb_id = get_pdb_id_from_filepath(imgt_filepath) 
    
    # --- STEP 1: Parse IMGT file for CDR3 and Headers ---
    try:
        with open(imgt_filepath, "r") as f:
            lines = f.readlines()
    except Exception:
        return [], "ReadError_IMGT", pdb_id

    header_info = []
    atom_cdr3_parts = defaultdict(list)
    last_seen_res = {}
    found_header_line = False
    
    for line in lines:
        if line.startswith("REMARK") and ("PAIRED_HL" in line or "SINGLE" in line):
            info = parse_header_flexible(line)
            if info:
                header_info.append(info)
                found_header_line = True
            continue

        if line.startswith("ATOM"):
            chain = line[21].upper()
            res_id_full = line[22:27]
            
            if last_seen_res.get(chain) == res_id_full:
                continue
            last_seen_res[chain] = res_id_full

            res_name = line[17:20].strip()
            aa = AA_MAP.get(res_name, 'X')

            try:
                res_num = int(line[22:26])
                if CDR3_START <= res_num <= CDR3_END:
                    atom_cdr3_parts[chain].append(aa)
            except ValueError:
                pass
    
    # --- STEP 2: Get full SEQRES from Original PDB ---
    seqres_strings = parse_seqres_from_original_pdb(pdb_id, ORIGINAL_PDB_DIR)

    # --- STEP 3: Build Final Rows ---
    if not found_header_line or not header_info:
        return [], "No_Header_Found", pdb_id

    cdr3_strings = {k: "".join(v) for k, v in atom_cdr3_parts.items()}

    rows = []
    success = False
    
    for info in header_info:
        H, L, AG_list, mode = info["H"], info["L"], info["AG"], info["mode"]

        # Retrieve sequences (safely handle None for H or L)
        heavy_cdr3 = "".join(cdr3_strings.get(H, [])) if H else ""
        light_cdr3 = "".join(cdr3_strings.get(L, [])) if L else ""

        # --- VALIDATION LOGIC ---
        valid_entry = False

        if mode == "PAIRED":
            # Needs both valid
            if len(heavy_cdr3) >= 3 and len(light_cdr3) >= 3:
                valid_entry = True
        elif mode == "SINGLE_H":
            # Needs valid Heavy
            if len(heavy_cdr3) >= 3:
                valid_entry = True
        elif mode == "SINGLE_L":
            # Needs valid Light
            if len(light_cdr3) >= 3:
                valid_entry = True

        if not valid_entry:
            continue

        # --- ANTIGEN EXTRACTION ---
        for ag in AG_list:
            antigen_seq = seqres_strings.get(ag, "")
            
            if not antigen_seq: 
                continue 

            rows.append({
                "pdb_id": pdb_id,
                "type": mode,
                "heavy_cdr3": heavy_cdr3,
                "light_cdr3": light_cdr3,
                "antigen_chain": ag,
                "antigen_seq": antigen_seq 
            })
            success = True

    if success:
        return rows, "Success", pdb_id
    else:
        return [], "Antigen_SEQRES_Missing_or_Invalid_CDR3", pdb_id

# ---------------- RUNNER ----------------

def run_all(imgt_dir=IMGTPDB_DIR, original_dir=ORIGINAL_PDB_DIR, out_csv=OUT_CSV, simplified_out_csv=SIMPLIFIED_OUT_CSV):
    if not os.path.isdir(imgt_dir) or not os.path.isdir(original_dir):
        print(f"Error: Ensure both directories exist: {imgt_dir} and {original_dir}")
        return

    files = [os.path.join(imgt_dir, f) for f in os.listdir(imgt_dir) if f.lower().endswith(".pdb")]
    print(f"Processing {len(files)} IMGT files in parallel...")

    all_rows = []
    stats = Counter()
    error_samples = defaultdict(list)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(process_single_pdb, files)
        
        for rows, status, pdb_id in results:
            stats[status] += 1
            error_samples[status].append(pdb_id)
            if rows:
                all_rows.extend(rows)

    print("\n" + "="*60)
    print("FINAL PROCESSING REPORT")
    print("="*60)
    for status, count in stats.most_common():
        print(f"[{status}]: {count}")
        if status != "Success":
             print(f"   Sample IDs: {', '.join(error_samples[status][:3])}")

    # --- 1. WRITE FULL OUTPUT CSV ---
    df_full = pd.DataFrame(all_rows)
    df_full.to_csv(out_csv, index=False)
    print("\n" + "="*60)
    print(f"SUCCESS: Wrote {len(df_full)} rows to FULL output ({out_csv})")

    # --- 2. CREATE AND WRITE SIMPLIFIED ML OUTPUT CSV ---
    if not df_full.empty:
        # Create a copy for ML input
        df_ml = df_full.copy()
        
        # KEY LOGIC: "cdr3_seq" will be Heavy if present, otherwise Light (for SINGLE_L)
        # 1. Initialize with Heavy
        df_ml['cdr3_seq'] = df_ml['heavy_cdr3']
        
        # 2. Fill missing/empty Heavy spots with Light (specifically for SINGLE_L cases)
        # Note: In SINGLE_H or PAIRED, heavy_cdr3 is populated. In SINGLE_L, heavy_cdr3 is "".
        df_ml.loc[df_ml['type'] == 'SINGLE_L', 'cdr3_seq'] = df_ml.loc[df_ml['type'] == 'SINGLE_L', 'light_cdr3']
        
        df_ml = df_ml.rename(columns={'pdb_id': 'pdb'})
        df_ml['label'] = 1
        
        # Select final columns
        df_simplified = df_ml[['pdb', 'antigen_seq', 'cdr3_seq', 'label']]
        
        df_simplified.to_csv(simplified_out_csv, index=False)
        print(f"SUCCESS: Wrote {len(df_simplified)} rows to SIMPLIFIED ML output ({simplified_out_csv})")
    else:
        print("WARNING: No rows generated, skipping simplified CSV.")
    
    print("="*60)

if __name__ == "__main__":
    run_all()