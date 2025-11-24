# -*- coding: utf-8 -*-
import os
import re
import pandas as pd


# ---------- CONFIG ----------
IMGTPDB_DIR = "data/imgt_pdbs"
OUT_CSV = "cdr3_antigen_dataset.csv"

# ---------- UTILS ----------
# mapping from 3-letter to 1-letter amino acids
AA_MAP = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'SEC':'U','PYL':'O'
}

def three_to_one(res):
    return AA_MAP.get(res.upper(), 'X')

# ---------- PARSING FUNCTIONS ----------
def parse_paired_hl_lines(lines):
    """
    Extracts all antibody-antigen pairs from REMARK 5 PAIRED_HL lines.

    Returns list of dicts like:
    [
       {"H": "B", "L": "A", "AG": ["C"]},
       {"H": "E", "L": "D", "AG": ["F"]}
    ]
    """
    pairs = []

    pattern = re.compile(
        r"PAIRED_HL\s+HCHAIN=([A-Za-z0-9])\s+LCHAIN=([A-Za-z0-9])\s+AGCHAIN=([\w,]+)"
    )

    for line in lines:
        if "PAIRED_HL" in line:
            m = pattern.search(line)
            if m:
                h = m.group(1)
                l = m.group(2)
                ag_raw = m.group(3)
                ag_list = ag_raw.split(",")  # Option 1: one row per antigen chain
                pairs.append({"H": h, "L": l, "AG": ag_list})

    return pairs


def parse_sequences_from_pdb(lines):
    """
    Reads ATOM records but collapses all atoms for the same residue number
    into ONE amino acid.
    seqs[chain] = [(res_num, aa)] sorted by res_num
    """
    seqs = {}
    seen = set()   # (chain, res_num) to avoid duplicates

    for line in lines:
        if not line.startswith("ATOM"):
            continue
        
        res_name = line[17:20].strip()
        chain = line[21].strip()
        res_num = line[22:27].strip()

        if not chain:
            continue
        if res_name.upper() not in AA_MAP:
            continue

        try:
            pos = int(res_num)
        except:
            continue
        
        key = (chain, pos)
        if key in seen:
            continue  # skip repeated atoms for same residue

        seen.add(key)

        aa = three_to_one(res_name)
        if chain not in seqs:
            seqs[chain] = []
        seqs[chain].append((pos, aa))

    # sort by IMGT numbering
    for ch in seqs:
        seqs[ch] = sorted(seqs[ch], key=lambda x: x[0])

    return seqs



def extract_cdr3(seq_list, start=105, end=117):
    """Extract IMGT CDR3: residues with positions in [105,117]"""
    cdr = [aa for (pos, aa) in seq_list if start <= pos <= end]
    return "".join(cdr)


# ---------- MAIN EXTRACTION ----------
def process_pdb(filepath):
    pdb_id = os.path.splitext(os.path.basename(filepath))[0].upper()

    with open(filepath, "r") as f:
        lines = f.readlines()

    # 1) Parse chain pairing info
    pairs = parse_paired_hl_lines(lines)
    if not pairs:
        return []  # skip PDBs with no PAIRED_HL

    # 2) Parse amino acid sequences
    seqs = parse_sequences_from_pdb(lines)

    rows = []

    for pair in pairs:
        H = pair["H"]
        L = pair["L"]
        AG_list = pair["AG"]

        if H not in seqs or L not in seqs:
            continue

        heavy_seq = seqs[H]
        light_seq = seqs[L]

        heavy_cdr3 = extract_cdr3(heavy_seq)
        light_cdr3 = extract_cdr3(light_seq)

        # One row per antigen chain
        for ag in AG_list:
            if ag not in seqs:
                continue

            antigen_seq = "".join([aa for (pos, aa) in seqs[ag]])

            rows.append({
                "pdb_id": pdb_id,
                "heavy_cdr3": heavy_cdr3,
                "light_cdr3": light_cdr3,
                "antigen_chain": ag,
                "antigen_seq": antigen_seq
            })

    return rows


# ---------- RUN ----------
all_rows = []

for fname in os.listdir(IMGTPDB_DIR):
    if fname.lower().endswith(".pdb"):
        fpath = os.path.join(IMGTPDB_DIR, fname)
        rows = process_pdb(fpath)
        all_rows.extend(rows)

df = pd.DataFrame(all_rows)
df.to_csv(OUT_CSV, index=False)

print(f"Wrote {len(df)} rows to {OUT_CSV}")