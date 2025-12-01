import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score

EMB = "data/pair_embeddings.npy"
LAB = "data/labels.npy"
CKPT = "score_model.pt"   # or "score_model_safe.pt" if you convert later (DONT FORGET ANDREW)

# ---- model must match training ----
class MLP(nn.Module):
    def __init__(self, in_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 512), nn.ReLU(),
            nn.Linear(512, 128), nn.ReLU(),
            nn.Linear(128, 1), nn.Sigmoid()
        )
    def forward(self, x): return self.net(x)

def to_numpy(x):
    """Support mu/sd saved as numpy arrays or torch tensors."""
    if isinstance(x, np.ndarray):
        return x
    if torch.is_tensor(x):
        return x.detach().cpu().numpy()
    raise TypeError(f"Unsupported type for mu/sd: {type(x)}")

def main():
    pairs = np.load(EMB)                 # (N,2,D)
    y = np.load(LAB).astype(np.int32)    # (N,)
    A, C = pairs[:,0,:], pairs[:,1,:]
    X = np.concatenate([A, C], axis=1).astype(np.float32)

    # note to future self: allow non-tensor objects in ckpt (pytorch 2.6 default changed xd)
    ckpt = torch.load(CKPT, map_location="cpu", weights_only=False)

    mu = to_numpy(ckpt["mu"])
    sd = to_numpy(ckpt["sd"])
    X = (X - mu) / (sd + 1e-8)

    model = MLP(X.shape[1])
    model.load_state_dict(ckpt["state_dict"])
    model.eval()

    with torch.no_grad():
        s = model(torch.from_numpy(X)).squeeze(1).numpy()

    print("ROC-AUC:", roc_auc_score(y, s))

if __name__ == "__main__":
    main()
