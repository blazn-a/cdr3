import os, argparse, numpy as np, pandas as pd, torch, esm
from tqdm import tqdm

CSV = "data/antigen_cdr3_pairs.csv"
EMB_OUT = "data/pair_embeddings.npy"
LAB_OUT = "data/labels.npy"

def load_model(name):
    if name == "t6_8M":
        return esm.pretrained.esm2_t6_8M_UR50D()
    elif name == "t12_35M":
        return esm.pretrained.esm2_t12_35M_UR50D()
    else:
        raise ValueError("model must be one of: t6_8M, t12_35M")

def embed_many(model, alphabet, seqs, batch_size=4, device="cpu"):
    model.eval()
    converter = alphabet.get_batch_converter()
    out = np.zeros((len(seqs), model.embed_dim), dtype=np.float32)
    with torch.no_grad():
        for i in tqdm(range(0, len(seqs), batch_size), desc="Embedding"):
            batch = seqs[i:i+batch_size]
            batch = [(str(i+k), s) for k, s in enumerate(batch)]
            _, _, toks = converter(batch)
            toks = toks.to(device, non_blocking=True)
            with torch.cuda.amp.autocast(enabled=(device!="cpu")):
                rep = model(toks, repr_layers=[model.num_layers])["representations"][model.num_layers]
            for j, s in enumerate(seqs[i:i+batch_size]):
                out[i+j] = rep[j, 1:1+len(s)].mean(0).detach().cpu().numpy()
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=CSV)
    ap.add_argument("--emb_out", default=EMB_OUT)
    ap.add_argument("--lab_out", default=LAB_OUT)
    ap.add_argument("--model", default="t6_8M", choices=["t6_8M","t12_35M"])
    ap.add_argument("--batch", type=int, default=4)          # keep small on CPU
    ap.add_argument("--threads", type=int, default=2)        # limit BLAS threads
    ap.add_argument("--force_cpu", action="store_true")
    args = ap.parse_args()

    # limit CPU thread explosion (so my pc doesnt die)
    torch.set_num_threads(max(1, args.threads))
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)

    print(f"Loading {args.model}…")
    model, alphabet = load_model(args.model)
    device = "cpu" if args.force_cpu or (not torch.cuda.is_available()) else "cuda"
    if device == "cuda": model = model.to("cuda")
    print(f"Device: {device}, batch={args.batch}, threads={args.threads}")

    df = pd.read_csv(args.csv)
    ant = df["antigen_seq"].astype(str).tolist()
    cdr = df["cdr3_seq"].astype(str).tolist()

    # embed in two passes to reduce peak RAM
    print("Embedding antigens…")
    A = embed_many(model, alphabet, ant, batch_size=args.batch, device=device)
    print("Embedding CDR3s…")
    C = embed_many(model, alphabet, cdr, batch_size=args.batch, device=device)

    pairs = np.stack([A, C], axis=1).astype(np.float32)
    os.makedirs(os.path.dirname(args.emb_out), exist_ok=True)
    np.save(args.emb_out, pairs)
    print(f"Saved {args.emb_out} with shape {pairs.shape}")

    # save labels if present. otherwise make dummy 1s
    if "label" in df.columns:
        y = df["label"].astype(np.float32).values
    else:
        y = np.ones(len(df), dtype=np.float32)
    np.save(args.lab_out, y)
    print(f"Saved {args.lab_out} with shape {y.shape}")

if __name__ == "__main__":
    main()
