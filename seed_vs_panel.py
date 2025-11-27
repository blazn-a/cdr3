import sys, csv, numpy as np, torch, torch.nn as nn, esm

class MLP(nn.Module):
    def __init__(self,d): super().__init__(); self.net=nn.Sequential(
        nn.Linear(d,512),nn.ReLU(),nn.Linear(512,128),nn.ReLU(),nn.Linear(128,1),nn.Sigmoid())
    def forward(self,x): return self.net(x)

def load_ckpt(p,in_dim):
    ck=torch.load(p,map_location="cpu",weights_only=False)
    m=MLP(in_dim); m.load_state_dict(ck["state_dict"]); m.eval()
    return m, np.array(ck["mu"]), np.array(ck["sd"])

def embed_mean(model,alphabet,seq):
    conv=alphabet.get_batch_converter()
    with torch.no_grad():
        _,_,tok=conv([("x",seq)])
        rep=model(tok,repr_layers=[model.num_layers])["representations"][model.num_layers]
        return rep[0,1:1+len(seq)].mean(0).cpu().numpy().astype("float32")

antigen = sys.argv[1].strip()
seed    = sys.argv[2].strip()
panel_csv = sys.argv[3].strip()

esm_model, alphabet = esm.pretrained.esm2_t6_8M_UR50D(); esm_model.eval()
a = embed_mean(esm_model, alphabet, antigen)
model, mu, sd = load_ckpt("score_model.pt", a.shape[0]*2)

def score(c):
    x = np.concatenate([a,c]).astype("float32")
    xn= (x - mu)/(sd+1e-8)
    with torch.no_grad():
        return float(model(torch.from_numpy(xn).unsqueeze(0)).item())

seed_c = embed_mean(esm_model, alphabet, seed)
seed_p = score(seed_c)

rows=[]; 
with open(panel_csv, newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        s = r.get("cdr3") or r.get("cdr3_seq") or r.get("sequence")
        if not s: continue
        pc = embed_mean(esm_model, alphabet, s)
        rows.append((s, score(pc)))

rows.sort(key=lambda t: -t[1])
top = rows[:5]
print(f"Seed score: {seed_p:.3f}  ({seed})")
print("Top panel candidates:")
for s,p in top: print(f"  {p:.3f}  {s}")
better = sum(p>seed_p for _,p in rows)
print(f"{better} / {len(rows)} candidates scored higher than the seed.")
