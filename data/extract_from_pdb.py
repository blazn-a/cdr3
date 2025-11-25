# -*- coding: ascii -*-
import os
import re
import pandas as pd

# ---------------- CONFIG ----------------
IMGTPDB_DIR = "data/imgt_pdbs"
OUT_CSV = "cdr3_antigen_dataset.csv"

# ---------------- AA MAP ----------------
AA_MAP = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'SEC':'U','PYL':'O'
}

def three_to_one(res):
    return AA_MAP.get(res.upper(), 'X')

# ---------------- PAIRED_HL PARSER ----------------
def parse_paired_hl_lines(lines):
    """
    Parse REMARK 5 PAIRED_HL lines. Returns list of dicts:
      {"H": "B", "L": "A", "AG": ["C","D"]}
    """
    pairs = []
    pattern = re.compile(
        r"PAIRED_HL\s+HCHAIN=([A-Za-z0-9])\s+LCHAIN=([A-Za-z0-9])\s+AGCHAIN=([\w,]+)"
    )
    for line in lines:
        if "PAIRED_HL" in line:
            m = pattern.search(line)
            if m:
                pairs.append({
                    "H": m.group(1),
                    "L": m.group(2),
                    "AG": m.group(3).split(",")
                })
    return pairs

# ---------------- SEQRES PARSER ----------------
def parse_seqres_sequences(lines):
    """
    Return dict: chain -> full sequence string from SEQRES.
    """
    seqres = {}
    for line in lines:
        if not line.startswith("SEQRES"):
            continue
        # chain id typically at column 11 (0-indexed)
        if len(line) < 12:
            continue
        chain = line[11].strip()
        if not chain:
            continue
        residues_part = line[19:].strip()
        if not residues_part:
            continue
        three_letter = residues_part.split()
        one_letter = [three_to_one(x) for x in three_letter]
        seqres.setdefault(chain, []).extend(one_letter)
    for ch in list(seqres.keys()):
        seqres[ch] = "".join(seqres[ch])
    return seqres

# ---------------- ATOM PARSER (PDB/IMGT order) ----------------
def parse_sequences_from_pdb(lines):
    """
    Parse ATOM records **in file order** and build seqs[chain] = [(pos_raw, aa), ...].
    pos_raw is taken from columns 22:27 (IMGT renumbered position string, may include letters).
    Deduplicate only exact pos_raw (keep first occurrence).
    """
    seqs_raw = {}

    for line in lines:
        if not line.startswith("ATOM"):
            continue
        # safe slicing
        res_name = line[17:20].strip()
        chain = line[21].strip()
        pos_raw = line[22:27].strip()  # e.g., "112", "112A", "112AA", or ""
        if not chain:
            continue
        if not res_name:
            continue
        if res_name.upper() not in AA_MAP:
            continue
        aa = three_to_one(res_name)
        seqs_raw.setdefault(chain, []).append((pos_raw, aa))

    # Deduplicate exact pos_raw strings while preserving first occurrence
    seqs = {}
    for ch, items in seqs_raw.items():
        seen = set()
        out = []
        for pos_raw, aa in items:
            # use full pos_raw as key; empty pos_raw allowed but will be deduped too
            key = pos_raw
            if key in seen:
                continue
            seen.add(key)
            out.append((pos_raw, aa))
        seqs[ch] = out

    return seqs

# ---------------- HELPER: numeric prefix ----------------
def numeric_prefix(pos_raw):
    """
    Return the numeric prefix of an IMGT position string.
    e.g., "112" -> 112, "112A" -> 112, "112AA" -> 112
    Returns None if no leading digits found.
    """
    if not pos_raw:
        return None
    i = 0
    while i < len(pos_raw) and pos_raw[i].isdigit():
        i += 1
    if i == 0:
        return None
    try:
        return int(pos_raw[:i])
    except:
        return None

# ---------------- CDR3 EXTRACTION ----------------
def extract_cdr3_from_seq_list(seq_list, start=105, end=117):
    """
    Given seq_list = [(pos_raw, aa), ...] in IMGT-provided order,
    return CDR3 string consisting of all aa whose numeric_prefix in [start,end].
    Keeps order exactly as in seq_list. Includes lettered positions.
    """
    out = []
    for pos_raw, aa in seq_list:
        num = numeric_prefix(pos_raw)
        if num is None:
            continue
        if start <= num <= end:
            out.append(aa)
    return "".join(out)

# ---------------- MAIN PROCESS PDB ----------------
def process_pdb(filepath):
    """
    Process a single IMGT-renumbered PDB file and return rows:
    [{"pdb_id":..., "heavy_cdr3":..., "light_cdr3":..., "antigen_chain":..., "antigen_seq":...}, ...]
    """
    pdb_id = os.path.splitext(os.path.basename(filepath))[0].upper()
    with open(filepath, "r") as f:
        lines = f.readlines()

    pairs = parse_paired_hl_lines(lines)
    if not pairs:
        return []

    seqs_atom = parse_sequences_from_pdb(lines)
    seqres = parse_seqres_sequences(lines)

    rows = []
    for pair in pairs:
        H = pair["H"]
        L = pair["L"]
        AG_list = pair["AG"]

        # require heavy and light present in ATOM to extract CDR3
        if H not in seqs_atom or L not in seqs_atom:
            continue

        heavy_cdr3 = extract_cdr3_from_seq_list(seqs_atom[H])
        light_cdr3 = extract_cdr3_from_seq_list(seqs_atom[L])

        for ag in AG_list:
            # prefer SEQRES for antigen full sequence, else fallback to ATOM order
            if ag in seqres:
                antigen_seq = seqres[ag]
            else:
                if ag in seqs_atom:
                    antigen_seq = "".join([aa for (pos, aa) in seqs_atom[ag]])
                else:
                    # antigen chain missing; skip
                    continue

            rows.append({
                "pdb_id": pdb_id,
                "heavy_cdr3": heavy_cdr3,
                "light_cdr3": light_cdr3,
                "antigen_chain": ag,
                "antigen_seq": antigen_seq
            })

    return rows

# ---------------- RUN OVER DIRECTORY ----------------
def run_all(imgt_dir=IMGTPDB_DIR, out_csv=OUT_CSV):
    all_rows = []
    if not os.path.isdir(imgt_dir):
        raise RuntimeError("IMGTPDB_DIR not found: " + str(imgt_dir))

    for fname in sorted(os.listdir(imgt_dir)):
        if not fname.lower().endswith(".pdb"):
            continue
        fpath = os.path.join(imgt_dir, fname)
        try:
            rows = process_pdb(fpath)
            all_rows.extend(rows)
        except Exception as e:
            # print minimal error and continue
            print("Error processing", fname, ":", repr(e))

    df = pd.DataFrame(all_rows)
    df.to_csv(out_csv, index=False)
    print("Wrote {} rows to {}".format(len(df), out_csv))

# ---------------- MAIN ----------------
if __name__ == "__main__":
    run_all()
