# optimize_cdr3.py  (3rd patch xd)
import argparse, random, numpy as np, torch, torch.nn as nn
from tqdm import tqdm
import esm, csv, os

AA = "ACDEFGHIKLMNPQRSTVWY"

class MLP(nn.Module):
    def __init__(self, in_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 512), nn.ReLU(),
            nn.Linear(512, 128), nn.ReLU(),
            nn.Linear(128, 1), nn.Sigmoid()
        )
    def forward(self, x): return self.net(x)

def _to_numpy(x):
    if isinstance(x, np.ndarray): return x
    if torch.is_tensor(x): return x.detach().cpu().numpy()
    raise TypeError(f"Unsupported type for mu/sd: {type(x)}")

def load_scorer(ckpt_path: str, in_dim: int):
    print(f">>> loading checkpoint {ckpt_path} with weights_only=False")
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    mu = _to_numpy(ckpt["mu"])
    sd = _to_numpy(ckpt["sd"])
    model = MLP(in_dim)
    model.load_state_dict(ckpt["state_dict"])
    model.eval()
    return model, mu, sd

def load_esm(model_name="t6_8M"):
    if model_name == "t6_8M":
        model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    elif model_name == "t12_35M":
        model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
    else:
        raise ValueError("model_name must be t6_8M or t12_35M")
    model.eval()
    return model, alphabet

def embed_mean(model, alphabet, seqs, batch=64):
    conv = alphabet.get_batch_converter()
    outs = []
    with torch.no_grad():
        for i in range(0, len(seqs), batch):
            batch_items = [(str(j), s) for j, s in enumerate(seqs[i:i+batch])]
            _, _, toks = conv(batch_items)
            rep = model(toks, repr_layers=[model.num_layers])["representations"][model.num_layers]
            for k, s in enumerate(seqs[i:i+batch]):
                outs.append(rep[k, 1:1+len(s)].mean(0).detach().cpu().numpy())
    return np.vstack(outs).astype(np.float32)

def score_pairs(scorer, mu, sd, antigen_emb, cdr_embs):
    X = np.concatenate(
        [np.repeat(antigen_emb[None, :], len(cdr_embs), axis=0), cdr_embs],
        axis=1
    ).astype(np.float32)
    X = (X - mu) / (sd + 1e-8)
    with torch.no_grad():
        scores = scorer(torch.from_numpy(X)).squeeze(1).numpy()
    return scores

# ---- 8/23/25 patched: keep leading 'C' and never mutate index 0
def mutate(seq, k=1):
    s = list(seq)
    L = len(s)
    for _ in range(k):
        i = random.randrange(1, L)  # never mutate first position
        s[i] = random.choice(AA)
    s[0] = 'C'  # enforce leading C
    return "".join(s)

def optimize(antigen_seq, start_cdr3, steps=200, beam=20, k_mut=1,
             esm_model_name="t6_8M", ckpt_path="score_model.pt",
             embed_batch=64, topk=5):
    esm_model, alphabet = load_esm(esm_model_name)
    antigen_emb = embed_mean(esm_model, alphabet, [antigen_seq], batch=embed_batch)[0]

    scorer, mu, sd = load_scorer(ckpt_path, in_dim=antigen_emb.size * 2)

    # enforce leading C on the seed
    if not start_cdr3.startswith('C'):
        start_cdr3 = 'C' + start_cdr3[1:]

    pool = [start_cdr3]
    pool_scores = np.array([0.0])

    for _ in tqdm(range(steps), desc="Optimizing"):
        candidates = set(pool)
        for s in pool:
            for _ in range(3):
                candidates.add(mutate(s, k=k_mut))
        candidates = list(candidates)
        cdr_embs = embed_mean(esm_model, alphabet, candidates, batch=embed_batch)
        scores = score_pairs(scorer, mu, sd, antigen_emb, cdr_embs)
        idx = np.argsort(-scores)[:beam]
        pool = [candidates[i] for i in idx]
        pool_scores = scores[idx]

    best_idx = np.argsort(-pool_scores)[:topk]
    best = [(pool[i], float(pool_scores[i])) for i in best_idx]
    return best, pool, pool_scores

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--antigen", required=True)
    ap.add_argument("--start_cdr3", default="CASSIRSSYEQYF")
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--beam", type=int, default=20)
    ap.add_argument("--k_mut", type=int, default=1)
    ap.add_argument("--esm", default="t6_8M", choices=["t6_8M","t12_35M"])
    ap.add_argument("--ckpt", default="score_model.pt")
    ap.add_argument("--embed_batch", type=int, default=64)
    ap.add_argument("--topk", type=int, default=5)
    ap.add_argument("--out_csv", default="")
    args = ap.parse_args()

    best, pool, pool_scores = optimize(
        antigen_seq=args.antigen,
        start_cdr3=args.start_cdr3,
        steps=args.steps,
        beam=args.beam,
        k_mut=args.k_mut,
        esm_model_name=args.esm,
        ckpt_path=args.ckpt,
        embed_batch=args.embed_batch,
        topk=args.topk,
    )

    print("\nTop suggestions:")
    for seq, sc in best:
        print(f"{sc:.3f}\t{seq}")

    if args.out_csv:
        os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
        import csv as _csv
        with open(args.out_csv, "w", newline="") as f:
            w = _csv.writer(f); w.writerow(["antigen_seq","cdr3_seq","score"])
            order = np.argsort(-pool_scores)
            for i in order:
                w.writerow([args.antigen, pool[i], f"{float(pool_scores[i]):.6f}"])
        print(f"Wrote {len(pool)} rows to {args.out_csv}")
