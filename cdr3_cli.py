import os, sys, json, time, subprocess, csv, random
from pathlib import Path

# -------- colored output --------
try:
    from colorama import init as colorama_init, Fore, Style
    colorama_init()
except Exception:
    class _N: RESET_ALL=""; CYAN=""; GREEN=""; YELLOW=""; RED=""
    Fore = Style = _N()

# -------- state I/O --------
STATE_FILE = "ui_state.json"

DEFAULTS = {
    # core
    "antigen": "",
    "seed_cdr3": "CARDRSTGYVYFDYW",
    "esm": "t12_35M",           # default to the stronger backbone
    "checkpoint": "score_model.pt",

    # generation / filter / cluster
    "steps": 200,
    "beam": 50,
    "k_mut": 1,
    "keep_top": 50,
    "embed_batch": 64,
    "clusters": 12,
    "out_folder": "runs",

    # training
    "pairs_csv": "data/antigen_cdr3_pairs.csv",
    "epochs": 30,
    "train_batch": 128,
    "threads": 2,
    "use_gpu": False,

    # NEW: training upgrades
    "hard_neg": True,           # enable hard-negative mining (antigen x mismatched CDR3)
    "hn_factor": 1,             # how many mismatched negatives per original pair
    "pos_mult": 1,              # positive oversampling factor for class weighting (1 = none)

    # last outputs
    "last_raw": "",
    "last_filtered": "",
    "last_panel": ""
}

def load_state():
    if Path(STATE_FILE).exists():
        try:
            return {**DEFAULTS, **json.loads(Path(STATE_FILE).read_text())}
        except Exception:
            pass
    return DEFAULTS.copy()

def save_state(S):
    Path(STATE_FILE).write_text(json.dumps(S, indent=2))

def ensure_out_folder(S):
    out = S.get("out_folder") or "runs"
    S["out_folder"] = out
    os.makedirs(out, exist_ok=True)
    save_state(S)
    return out

# -------- UI helpers --------
def clear_screen():
    os.system("cls" if os.name == "nt" else "clear")

def header(title):
    print(Fore.CYAN + "=" * 30 + Style.RESET_ALL)
    print(Fore.CYAN + f" {title}" + Style.RESET_ALL)
    print(Fore.CYAN + "=" * 30 + Style.RESET_ALL)

def press_enter():
    try:
        input("\nPress Enter to continue…")
    except EOFError:
        pass

def show_config(S):
    print(Fore.YELLOW + "Current config:" + Style.RESET_ALL)
    print(f"  Antigen length : {len(S.get('antigen',''))}")
    print(f"  Seed CDR3      : {S['seed_cdr3']}")
    print(f"  ESM backbone   : {S['esm']}")
    print(f"  Checkpoint     : {S['checkpoint']}")
    print(f"  Steps/Beam/k   : {S['steps']} / {S['beam']} / {S['k_mut']}")
    print(f"  Keep top       : {S['keep_top']}")
    print(f"  Embed batch    : {S['embed_batch']}")
    print(f"  Clusters (k)   : {S['clusters']}")
    print(f"  Output folder  : {S['out_folder']}")
    print(f"  Pairs CSV      : {S['pairs_csv']}")
    print(f"  Train epochs   : {S['epochs']}")
    print(f"  Train batch    : {S['train_batch']}")
    print(f"  Threads        : {S['threads']}")
    print(f"  Use GPU        : {S['use_gpu']}")
    print(Fore.YELLOW + "Training upgrades:" + Style.RESET_ALL)
    print(f"  Hard negatives : {S['hard_neg']}  (factor={S['hn_factor']})")
    print(f"  Pos oversample : x{S['pos_mult']}")

def run_py(script, *args, check=True, capture=False):
    cmd = [sys.executable, script, *map(str, args)]
    print(Fore.CYAN + "$ " + " ".join(cmd) + Style.RESET_ALL)
    if capture:
        return subprocess.run(cmd, check=check, text=True, capture_output=True)
    return subprocess.run(cmd, check=check)

# -------- FASTA  --------
def write_fasta_from_csv(panel_csv: str, fasta_path: str):
    """Fallback FASTA writer if export_fasta.py isn't present."""
    n = 0
    with open(panel_csv, newline="", encoding="utf-8") as f_in, \
         open(fasta_path, "w", encoding="utf-8") as f_out:
        reader = csv.reader(f_in)
        header_row = next(reader, None)
        cdr3_idx = None
        if header_row:
            for i, h in enumerate(header_row):
                if str(h).strip().lower() in ("cdr3", "sequence", "seq"):
                    cdr3_idx = i; break
            if cdr3_idx is None and len(header_row) >= 2:
                cdr3_idx = 1
        if cdr3_idx is None:
            f_in.seek(0); reader = csv.reader(f_in); cdr3_idx = 1
        for i, row in enumerate(reader, start=1):
            if not row or len(row) <= cdr3_idx: continue
            seq = row[cdr3_idx].strip()
            if not seq: continue
            f_out.write(f">cand_{i}\n{seq}\n"); n += 1
    return n

# -------- actions --------
def action_edit(S):
    clear_screen(); header("Edit config"); show_config(S)
    print("\nEdit which field?")
    print("  1) Antigen sequence (paste)")
    print("  2) Seed CDR3")
    print("  3) ESM model (t6_8M or t12_35M)")
    print("  4) Checkpoint path")
    print("  5) Steps")
    print("  6) Beam width")
    print("  7) Mutations per step (k_mut)")
    print("  8) Keep top")
    print("  9) Embed batch size")
    print(" 10) Clusters (k)")
    print(" 11) Output folder")
    print(" 12) Pairs CSV (for training)")
    print(" 13) Train epochs")
    print(" 14) Train batch size")
    print(" 15) Threads")
    print(" 16) Toggle GPU (True/False)")
    print(" 17) Toggle hard negatives")
    print(" 18) Set hard-neg factor")
    print(" 19) Set positive oversample (class weighting)")
    print(" 20) Done")
    while True:
        ch = input("> ").strip()
        if ch == "1":
            S["antigen"] = input("Paste antigen (single line, no spaces):\n> ").strip()
        elif ch == "2":
            S["seed_cdr3"] = input("Seed CDR3:\n> ").strip()
        elif ch == "3":
            v = input("ESM backbone [t6_8M/t12_35M]: ").strip()
            if v in ("t6_8M","t12_35M"): S["esm"] = v
        elif ch == "4":
            S["checkpoint"] = input("Checkpoint path:\n> ").strip() or S["checkpoint"]
        elif ch == "5":
            S["steps"] = int(input("Steps: ").strip() or S["steps"])
        elif ch == "6":
            S["beam"] = int(input("Beam width: ").strip() or S["beam"])
        elif ch == "7":
            S["k_mut"] = int(input("Mutations per step (k_mut): ").strip() or S["k_mut"])
        elif ch == "8":
            S["keep_top"] = int(input("Keep top: ").strip() or S["keep_top"])
        elif ch == "9":
            S["embed_batch"] = int(input("Embed batch: ").strip() or S["embed_batch"])
        elif ch == "10":
            S["clusters"] = int(input("Clusters (k): ").strip() or S["clusters"])
        elif ch == "11":
            S["out_folder"] = input("Output folder: ").strip() or S["out_folder"]
        elif ch == "12":
            S["pairs_csv"] = input("Pairs CSV path:\n> ").strip() or S["pairs_csv"]
        elif ch == "13":
            S["epochs"] = int(input("Train epochs: ").strip() or S["epochs"])
        elif ch == "14":
            S["train_batch"] = int(input("Train batch: ").strip() or S["train_batch"])
        elif ch == "15":
            S["threads"] = int(input("Threads: ").strip() or S["threads"])
        elif ch == "16":
            S["use_gpu"] = not S["use_gpu"]
            print(f"use_gpu -> {S['use_gpu']}")
            time.sleep(0.7)
        elif ch == "17":
            S["hard_neg"] = not S["hard_neg"]
            print(f"hard_neg -> {S['hard_neg']}")
            time.sleep(0.7)
        elif ch == "18":
            S["hn_factor"] = max(0, int(input("Hard-neg factor (0..): ").strip() or S["hn_factor"]))
        elif ch == "19":
            S["pos_mult"]  = max(1, int(input("Positive oversample (>=1): ").strip() or S["pos_mult"]))
        elif ch == "20":
            break
        else:
            print("Enter 1–20")
        save_state(S); clear_screen(); header("Edit config"); show_config(S)
    save_state(S)

def action_optimize(S):
    clear_screen(); header("Optimize (generate candidates)")
    if not S.get("antigen"):
        print(Fore.RED + "No antigen set. Edit config first." + Style.RESET_ALL)
        press_enter(); return
    out_dir = ensure_out_folder(S)
    out_csv = str(Path(out_dir, f"opt_{time.strftime('%Y%m%d_%H%M%S')}.csv"))
    try:
        run_py("optimize_cdr3.py",
               "--antigen", S["antigen"],
               "--start_cdr3", S["seed_cdr3"],
               "--steps", S["steps"],
               "--beam", S["beam"],
               "--topk", 20,
               "--k_mut", S["k_mut"],
               "--esm", S["esm"],
               "--ckpt", S["checkpoint"],
               "--out_csv", out_csv)
        print(Fore.GREEN + f"Wrote {out_csv}" + Style.RESET_ALL)
        S["last_raw"] = out_csv; save_state(S)
    except subprocess.CalledProcessError:
        print(Fore.RED + "Optimization failed." + Style.RESET_ALL)
    press_enter()

def action_filter(S):
    clear_screen(); header("Filter candidates")
    if not S.get("last_raw"):
        print(Fore.RED + "No raw candidate CSV. Run optimize first." + Style.RESET_ALL)
        press_enter(); return
    inp = S["last_raw"]
    out_dir = ensure_out_folder(S)
    out_csv = str(Path(out_dir, f"{Path(inp).stem}_filtered.csv"))
    try:
        run_py("filter_cdr3s.py", "--inp", inp, "--out", out_csv, "--keep_top", S["keep_top"])
        print(Fore.GREEN + f"Wrote {out_csv}" + Style.RESET_ALL)
        S["last_filtered"] = out_csv; save_state(S)
    except subprocess.CalledProcessError:
        print(Fore.RED + "Filtering failed." + Style.RESET_ALL)
    press_enter()

def action_cluster(S):
    clear_screen(); header("Cluster into panel + export FASTA")
    if not S.get("last_filtered"):
        print(Fore.RED + "No filtered CSV. Run filtering first." + Style.RESET_ALL)
        press_enter(); return
    inp = S["last_filtered"]
    out_dir = ensure_out_folder(S)
    out_csv  = str(Path(out_dir, f"{Path(inp).stem}_panel.csv"))
    out_fa   = str(Path(out_dir, f"{Path(inp).stem}_panel.fasta"))
    try:
        run_py("cluster_candidates.py",
               "--inp", inp, "--out", out_csv,
               "--k", S["clusters"], "--esm", S["esm"])
        print(Fore.GREEN + f"Wrote {out_csv}" + Style.RESET_ALL)
        S["last_panel"] = out_csv; save_state(S)

        if Path("export_fasta.py").exists():
            run_py("export_fasta.py", "--inp", out_csv, "--out", out_fa)
            print(Fore.GREEN + f"Wrote {out_fa}" + Style.RESET_ALL)
        else:
            n = write_fasta_from_csv(out_csv, out_fa)
            print(Fore.GREEN + f"Wrote {out_fa} ({n} sequences)" + Style.RESET_ALL)
    except subprocess.CalledProcessError:
        print(Fore.RED + "Clustering or FASTA export failed." + Style.RESET_ALL)
    press_enter()

def action_quick(S):
    # 2 -> 3 -> 4 pipeline
    action_optimize(S)
    if S.get("last_raw"): action_filter(S)
    if S.get("last_filtered"): action_cluster(S)

def action_outputs(S):
    clear_screen(); header("Outputs")
    print("Last candidates CSV : ", S.get("last_raw") or "(none)")
    print("Last filtered CSV   : ", S.get("last_filtered") or "(none)")
    print("Last panel CSV      : ", S.get("last_panel") or "(none)")
    press_enter()

# ---------- training data aug. helpers (a little broken) ----------
def _augment_with_hard_negatives(emb_path: str, lab_path: str, factor: int):
    import numpy as np
    X = np.load(emb_path)     # shape (N, 2, D)
    y = np.load(lab_path)     # shape (N,), 1/0

    if factor <= 0 or X.shape[0] < 2:
        return 0

    N, two, D = X.shape
    assert two == 2, "Expected embeddings in shape (N,2,D) [antigen, cdr3]."

    A = X[:,0,:]     # (N,D) antigen
    C = X[:,1,:]     # (N,D) cdr3

    added = 0
    parts = [X]
    labels = [y]

    for _ in range(factor):
        perm = np.arange(N)
        # random permutation w/no fixed points (or at least minimize)
        while True:
            np.random.shuffle(perm)
            if not np.all(perm == np.arange(N)):
                break
        # mismatched pairs do antigen i with cdr3 perm[i]
        C_perm = C[perm]
        X_hn = np.stack([A, C_perm], axis=1)     # (N,2,D)
        y_hn = np.zeros(N, dtype=y.dtype)
        parts.append(X_hn); labels.append(y_hn)
        added += N

    X_aug = np.concatenate(parts, axis=0)
    y_aug = np.concatenate(labels, axis=0)

    # shuffle
    idx = np.random.permutation(X_aug.shape[0])
    np.save(emb_path, X_aug[idx])
    np.save(lab_path, y_aug[idx])
    return added

def _oversample_positives(emb_path: str, lab_path: str, mult: int):
    import numpy as np
    if mult <= 1: return 0
    X = np.load(emb_path)     # (N,2,D)
    y = np.load(lab_path)     # (N,)
    pos_idx = np.where(y == 1)[0]
    if pos_idx.size == 0: return 0

    # repeat positives (mult-1 additional times)
    extra = []
    for _ in range(mult - 1):
        extra.append(X[pos_idx])
    if extra:
        X_aug = np.concatenate([X] + extra, axis=0)
        y_aug = np.concatenate([y] + [np.ones_like(pos_idx, dtype=y.dtype) for _ in range(mult-1)], axis=0)
        # shuffle
        idx = np.random.permutation(X_aug.shape[0])
        np.save(emb_path, X_aug[idx])
        np.save(lab_path, y_aug[idx])
        return (mult - 1) * pos_idx.size
    return 0

def action_train(S):
    clear_screen(); header("Train (embeddings + scorer)")
    pairs = Path(S["pairs_csv"])
    if not pairs.exists():
        print(Fore.RED + f"Missing pairs CSV: {pairs}" + Style.RESET_ALL)
        print("Tip: place antigen/CDR3 pairs at data/antigen_cdr3_pairs.csv")
        press_enter(); return

    emb_out = "data/pair_embeddings.npy"
    lab_out = "data/labels.npy"

    # 1) embed pairs
    print(Fore.YELLOW + "\n[1/3] Embedding pairs with ESM…" + Style.RESET_ALL)
    try:
        args = [
            "--csv", str(pairs),
            "--emb_out", emb_out,
            "--lab_out", lab_out,
            "--batch", S["embed_batch"],
            "--threads", S["threads"],
            "--model", S["esm"]
        ]
        if not S["use_gpu"]:
            args.append("--force_cpu")
        run_py("esm_embeddings.py", *args)
    except subprocess.CalledProcessError:
        print(Fore.RED + "Embedding failed." + Style.RESET_ALL)
        press_enter(); return

    # 2) optional: data-level upgrades
    print(Fore.YELLOW + "\n[2/3] Applying training upgrades…" + Style.RESET_ALL)
    augmented = 0
    try:
        if S.get("hard_neg", True) and S.get("hn_factor", 1) > 0:
            added_hn = _augment_with_hard_negatives(emb_out, lab_out, int(S["hn_factor"]))
            augmented += added_hn
            if added_hn:
                print(Fore.GREEN + f"  + Hard negatives: {added_hn} rows added" + Style.RESET_ALL)
        if S.get("pos_mult", 1) > 1:
            added_pos = _oversample_positives(emb_out, lab_out, int(S["pos_mult"]))
            augmented += added_pos
            if added_pos:
                print(Fore.GREEN + f"  + Pos oversample x{S['pos_mult']}: +{added_pos} rows" + Style.RESET_ALL)
        if augmented == 0:
            print("  (no augmentation applied)")
    except Exception as e:
        print(Fore.RED + f"Augmentation step failed (continuing without): {e}" + Style.RESET_ALL)

    # 3) train scorer
    print(Fore.YELLOW + "\n[3/3] Training scorer…" + Style.RESET_ALL)
    try:
        run_py("score_model.py", "--epochs", S["epochs"], "--batch", S["train_batch"])
        # assume score_model.py writes score_model.pt in cwd
        S["checkpoint"] = "score_model.pt"
        save_state(S)
        print(Fore.GREEN + "Training complete. Using checkpoint: score_model.pt" + Style.RESET_ALL)
    except subprocess.CalledProcessError:
        print(Fore.RED + "Training failed." + Style.RESET_ALL)

    # quick AUC readout
    try:
        out = run_py("eval_auc.py", check=False, capture=True)
        if out and out.stdout:
            print("\n" + out.stdout.strip())
    except Exception:
        pass
    press_enter()

def action_evaluate(S):
    clear_screen(); header("Evaluate")
    # Global AUC
    try:
        print(Fore.YELLOW + "[AUC] eval_auc.py" + Style.RESET_ALL)
        out = subprocess.check_output([sys.executable, "eval_auc.py"], text=True)
        print(out.strip())
    except subprocess.CalledProcessError as e:
        print(Fore.RED + "eval_auc.py failed:\n" + (e.output or "") + Style.RESET_ALL)

    # Seed vs Panel
    if not S.get("last_panel"):
        print(Fore.RED + "Missing panel CSV (run clustering first)." + Style.RESET_ALL)
    elif not S.get("antigen") or not S.get("seed_cdr3"):
        print(Fore.RED + "Need antigen and seed_cdr3 set in config." + Style.RESET_ALL)
    else:
        try:
            print(Fore.YELLOW + "[Seed vs Panel] seed_vs_panel.py" + Style.RESET_ALL)
            out = subprocess.check_output(
                [sys.executable, "seed_vs_panel.py",
                 S["antigen"], S["seed_cdr3"], S["last_panel"]],
                text=True
            )
            print(out.strip())
        except subprocess.CalledProcessError as e:
            print(Fore.RED + "seed_vs_panel.py failed:\n" + (e.output or "") + Style.RESET_ALL)
    press_enter()

# -------- main loop --------
def main():
    S = load_state()
    ensure_out_folder(S)

    while True:
        clear_screen()
        header("CDR3 Optimizer - CLI UI")
        show_config(S)
        print("\nMenu:")
        print(" 1) Edit config")
        print(" 2) Optimize (generate candidates)")
        print(" 3) Filter candidates")
        print(" 4) Cluster into panel + export FASTA")
        print(" 5) Quick pipeline (2 -> 3 -> 4)")
        print(" 6) Show last outputs")
        print(" 7) Train (embed + fit scorer)")
        print(" 8) Evaluate (global AUC + seed vs panel)")
        print(" 9) Quit")
        choice = input("> ").strip()

        if   choice == "1": action_edit(S)
        elif choice == "2": action_optimize(S)
        elif choice == "3": action_filter(S)
        elif choice == "4": action_cluster(S)
        elif choice == "5": action_quick(S)
        elif choice == "6": action_outputs(S)
        elif choice == "7": action_train(S)
        elif choice == "8": action_evaluate(S)
        elif choice == "9": print("Bye."); return
        else:
            print(Fore.RED + "Invalid choice." + Style.RESET_ALL); time.sleep(0.8)

if __name__ == "__main__":
    main()
