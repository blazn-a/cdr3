import argparse, numpy as np, torch, torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from tqdm import tqdm

EMB="data/pair_embeddings.npy"
LAB="data/labels.npy"

class MLP(nn.Module):
    def __init__(self, in_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 512), nn.ReLU(),
            nn.Linear(512, 128), nn.ReLU(),
            nn.Linear(128, 1), nn.Sigmoid()
        )
    def forward(self, x): return self.net(x)

def eval_metrics(model, loader, device):
    model.eval(); ys=[]; ps=[]
    with torch.no_grad():
        for xb, yb in loader:
            p = model(xb.to(device)).cpu()
            ys.append(yb); ps.append(p)
    y = torch.cat(ys).numpy().ravel()
    p = torch.cat(ps).numpy().ravel()
    auc = roc_auc_score(y, p)
    acc = ((p>0.5) == (y>0.5)).mean()
    return acc, auc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--patience", type=int, default=6)
    ap.add_argument("--num_workers", type=int, default=0)   # keep 0 on Windows
    ap.add_argument("--verbose_eval", action="store_true")
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # --- load n build features ---
    pairs = np.load(EMB)                 # (N,2,D)
    labels = np.load(LAB).astype(np.float32)[:,None]  # (N,1)
    A, C = pairs[:,0,:], pairs[:,1,:]
    X = np.concatenate([A, C], axis=1).astype(np.float32)

    # standardize
    mu = X.mean(0, keepdims=True)
    sd = X.std(0, keepdims=True) + 1e-8
    X = (X - mu) / sd

    # split
    Xtr, Xte, ytr, yte = train_test_split(X, labels, test_size=0.2, random_state=42, stratify=labels)
    train_loader = DataLoader(
        TensorDataset(torch.from_numpy(Xtr), torch.from_numpy(ytr)),
        batch_size=args.batch, shuffle=True, num_workers=args.num_workers
    )
    val_loader = DataLoader(
        TensorDataset(torch.from_numpy(Xte), torch.from_numpy(yte)),
        batch_size=args.batch, shuffle=False, num_workers=args.num_workers
    )

    # model/opt
    model = MLP(X.shape[1]).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    bce = nn.BCELoss()

    best_auc, bad, best_state = 0.0, 0, None

    # --- training loop with progress bar ---
    for ep in range(1, args.epochs+1):
        model.train()
        running = 0.0
        pbar = tqdm(train_loader, desc=f"epoch {ep:02d}/{args.epochs}", leave=False)
        for xb, yb in pbar:
            xb, yb = xb.to(device), yb.to(device)
            pred = model(xb)
            loss = bce(pred, yb)
            opt.zero_grad(); loss.backward(); opt.step()
            running += loss.item() * xb.size(0)
            pbar.set_postfix(loss=f"{loss.item():.4f}")

        train_loss = running / len(train_loader.dataset)
        val_acc, val_auc = eval_metrics(model, val_loader, device)
        print(f"epoch {ep:02d} | loss {train_loss:.4f} | val acc {val_acc*100:4.1f}% | val AUC {val_auc:.3f}")

        if args.verbose_eval:
            # optional: print a few score samples (only for visuals really so dont add)
            pass

        # early stopping on AUC
        if val_auc > best_auc:
            best_auc = val_auc; bad = 0; best_state = model.state_dict()
        else:
            bad += 1
            if bad >= args.patience:
                print(f"Early stopping at epoch {ep} (best AUC {best_auc:.3f})")
                break

    # save best model state (just feed the data more antigens until its good asl)
    if best_state is not None:
        model.load_state_dict(best_state)
    torch.save({"state_dict": model.state_dict(), "mu": mu, "sd": sd}, "score_model.pt")
    print("saved score_model.pt")

if __name__ == "__main__":
    main()
