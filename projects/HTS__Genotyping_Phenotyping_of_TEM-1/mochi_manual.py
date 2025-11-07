
# mochi_manual_fixed.py
# Modified version with bugfixes and anti-overfitting measures for order2.
# Saves: model checkpoints and prints training logs.
# Usage: place your input table as "example.tsv" (tab-separated) in the same folder,
# or edit the DATA_PATH variable to point to your file.

import os
import sys
from collections import Counter, defaultdict
import math
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.amp import autocast, GradScaler
import scipy.sparse as sp

# ========== USER PARAMETERS (edit if desired) ==========
DATA_PATH = "/mnt/c/Users/20145/Desktop/kinetics/pymochi_catal_nostop_dataset_weighted.tsv"   # replace with your actual path if different
MIN_PAIR_SUPPORT = 3        # minimum co-occurrence of a pair to be kept
L1_PAIR = 1e-3              # L1 penalty on pairwise coefficients
WEIGHT_DECAY = 1e-6         # L2 (weight decay) for optimizer
LR_ORDER1 = 1e-3
LR_ORDER2 = 1e-3
N_EPOCHS_ORDER1 = 2000
N_EPOCHS_ORDER2 = 5000
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
R = 0.001987                # kcal/(mol*K) - correct physical constant
T = 303.0                   # temperature in K (default 30°C)
RT = R * T

# ========== Utility functions ==========
def find_mutations(seq, wt):
    """Return list of (pos, mut_aa) for positions where seq differs from wt."""
    muts = []
    for i, (a, b) in enumerate(zip(seq, wt)):
        if a != b:
            muts.append((i, a))
    return muts

def find_mutation_pairs_for_seq(seq, wt):
    muts = find_mutations(seq, wt)
    pairs = []
    for i in range(len(muts)):
        for j in range(i+1, len(muts)):
            # canonical order for pair key
            a = muts[i]; b = muts[j]
            if a <= b:
                pairs.append((a, b))
            else:
                pairs.append((b, a))
    return pairs

# ========== Read data ==========
if not os.path.exists(DATA_PATH):
    raise FileNotFoundError(f"Input data file not found: {DATA_PATH}\nPlease provide a tab-separated file with columns: aa_seq, fitness, sigma, WT (optional)")

df = pd.read_csv(DATA_PATH, sep=None, engine='python')  # auto-detect sep
# normalize column names
df.columns = [c.strip() for c in df.columns]

if 'aa_seq' not in df.columns or 'fitness' not in df.columns or 'sigma' not in df.columns:
    raise ValueError("Input file must contain columns: 'aa_seq', 'fitness', 'sigma'")

# choose WT sequence
if 'WT' in df.columns:
    wt_rows = df[df['WT'].astype(bool)]
    if len(wt_rows) == 0:
        wt_seq = df['aa_seq'].iloc[0].strip()
    else:
        wt_seq = wt_rows['aa_seq'].iloc[0].strip()
else:
    wt_seq = df['aa_seq'].iloc[0].strip()

mut_seqs = [s.strip() for s in df['aa_seq'].tolist()]
fitness_np = np.array(df['fitness'].astype(float).tolist(), dtype=float)
sigma_np = np.array(df['sigma'].astype(float).tolist(), dtype=float)

N = len(mut_seqs)
print(f"Loaded {N} sequences. WT length = {len(wt_seq)}")

# ========== Build one-hot feature mapping (position+AA -> index) ==========
# Collect all single-mutation types observed across dataset.
single_mut_set = set()
for s in mut_seqs:
    for pos, aa in find_mutations(s, wt_seq):
        single_mut_set.add((pos, aa))
single_mut_list = sorted(single_mut_set)
M = len(single_mut_list)
mut_to_idx = {m: i for i, m in enumerate(single_mut_list)}
print(f"Total distinct single-mutation features (M): {M}")

# Build dense X (N x M) as float32 tensor
X_np = np.zeros((N, M), dtype=np.float32)
for i, s in enumerate(mut_seqs):
    for pos, aa in find_mutations(s, wt_seq):
        idx = mut_to_idx[(pos, aa)]
        X_np[i, idx] = 1.0

# ========== Build pair list and filter by support ==========
pair_counter = Counter()
for s in mut_seqs:
    for p in find_mutation_pairs_for_seq(s, wt_seq):
        pair_counter[p] += 1

filtered_pairs = [p for p, c in pair_counter.items() if c >= MIN_PAIR_SUPPORT]
filtered_pairs = sorted(filtered_pairs)
P = len(filtered_pairs)
pair_to_idx = {p: i for i, p in enumerate(filtered_pairs)}
print(f"Kept {P} pairwise features after min_support={MIN_PAIR_SUPPORT} filter (from {len(pair_counter)} total pairs)")

# Construct sparse X_pair (COO) of shape (N, P)
rows, cols, vals = [], [], []
for i, s in enumerate(mut_seqs):
    for p in find_mutation_pairs_for_seq(s, wt_seq):
        if p in pair_to_idx:
            rows.append(i)
            cols.append(pair_to_idx[p])
            vals.append(1.0)
if P > 0:
    X_pair_sp = sp.coo_matrix((np.array(vals, dtype=np.float32), (np.array(rows), np.array(cols))), shape=(N, P))
else:
    X_pair_sp = sp.coo_matrix((N, 0))

# ========== Convert to PyTorch tensors and move to device ==========
X = torch.tensor(X_np, dtype=torch.float32, device=DEVICE)
fitness = torch.tensor(fitness_np, dtype=torch.float32, device=DEVICE)
sigma = torch.tensor(sigma_np, dtype=torch.float32, device=DEVICE)
weight = 1.0 / (sigma + 1e-12)  # weighted L1 as in MoCHI
weight = weight.to(DEVICE)

# For X_pair (sparse) create torch.sparse_coo_tensor on device
if P > 0:
    indices = torch.LongTensor(np.vstack((rows, cols))).to(DEVICE)
    values = torch.FloatTensor(np.array(vals, dtype=np.float32)).to(DEVICE)
    X_pair = torch.sparse_coo_tensor(indices, values, size=(N, P)).coalesce().to(DEVICE)
else:
    X_pair = None

# ========== Model definitions ==========
class MoCHI_TwoState_Order1(nn.Module):
    def __init__(self, M, R=RT):
        super().__init__()
        self.R = R
        self.theta = nn.Parameter(torch.zeros(M, dtype=torch.float32))
        self.phi0 = nn.Parameter(torch.tensor(0.0, dtype=torch.float32))
        self.a = nn.Parameter(torch.tensor(1.0, dtype=torch.float32))
        self.b = nn.Parameter(torch.tensor(0.0, dtype=torch.float32))

    def forward(self, X):
        phi = self.phi0 + X @ self.theta
        p = 1.0 / (1.0 + torch.exp(phi / self.R))
        yhat = self.a * p + self.b
        return yhat, phi

class MoCHI_TwoState_Order2(nn.Module):
    def __init__(self, M, P, R=RT):
        super().__init__()
        self.R = R
        self.theta = nn.Parameter(torch.zeros(M, dtype=torch.float32))
        self.phi_pair = nn.Parameter(torch.zeros(P, dtype=torch.float32))
        self.phi0 = nn.Parameter(torch.tensor(0.0, dtype=torch.float32))
        self.a = nn.Parameter(torch.tensor(1.0, dtype=torch.float32))
        self.b = nn.Parameter(torch.tensor(0.0, dtype=torch.float32))

    def forward(self, X, X_pair):
        # X: dense (N,M), X_pair: sparse (N,P) torch.sparse_coo_tensor
        if X_pair is None or X_pair.numel() == 0:
            pair_term = 0.0
        else:
            # torch.sparse.mm 不支持 float16 on CUDA.
            # 在这里临时禁用 autocast，以确保 sparse.mm 在 float32 上执行.
            # 使用 torch.amp.autocast(enabled=False) 以和训练时的 autocast 协同工作。
            with autocast(enabled=False):
                pair_term = torch.sparse.mm(X_pair, self.phi_pair.unsqueeze(1)).squeeze(1)
        phi = self.phi0 + X @ self.theta + pair_term
        p = 1.0 / (1.0 + torch.exp(phi / self.R))
        yhat = self.a * p + self.b
        return yhat, phi


# ========== Training order1 ==========
model1 = MoCHI_TwoState_Order1(M, R=RT).to(DEVICE)
optimizer1 = optim.Adam(model1.parameters(), lr=LR_ORDER1, weight_decay=WEIGHT_DECAY)
scaler1 = GradScaler()

print("Training order1 model...")
best_loss1 = float('inf')
for epoch in range(1, N_EPOCHS_ORDER1 + 1):
    optimizer1.zero_grad(set_to_none=True)
    with autocast():
        yhat1, phi1 = model1(X)
        data_loss1 = torch.mean(torch.abs((fitness - yhat1) * weight))
    scaler1.scale(data_loss1).backward()
    scaler1.step(optimizer1)
    scaler1.update()

    if epoch % 200 == 0 or epoch == 1:
        print(f"[Order1] epoch {epoch:5d} loss={data_loss1.item():.6f}")
        if data_loss1.item() < best_loss1:
            best_loss1 = data_loss1.item()

# ========== Training order2 ==========
model2 = MoCHI_TwoState_Order2(M, P, R=RT).to(DEVICE)
optimizer2 = optim.Adam(model2.parameters(), lr=LR_ORDER2, weight_decay=WEIGHT_DECAY)
scaler2 = GradScaler()
print("Training order2 model...")
best_loss2 = float('inf')
for epoch in range(1, N_EPOCHS_ORDER2 + 1):
    optimizer2.zero_grad(set_to_none=True)
    with autocast():
        yhat2, phi2 = model2(X, X_pair)
        data_loss2 = torch.mean(torch.abs((fitness - yhat2) * weight))
        l1_pen = L1_PAIR * torch.norm(model2.phi_pair, p=1) if P > 0 else 0.0
        loss2 = data_loss2 + l1_pen
    scaler2.scale(loss2).backward()
    scaler2.step(optimizer2)
    scaler2.update()

    if epoch % 500 == 0 or epoch == 1:
        nonzero_pairs = int((model2.phi_pair.abs() > 1e-3).sum().item()) if P > 0 else 0
        print(f"[Order2] epoch {epoch:5d} data_loss={data_loss2.item():.6f} l1={l1_pen.item() if P>0 else 0.0:.6f} nonzero_pairs={nonzero_pairs}/{P}")
        if data_loss2.item() < best_loss2:
            best_loss2 = data_loss2.item()

# ========== Final outputs ==========
with torch.no_grad():
    yhat_final, phi_final = model2(X, X_pair)
    df['predicted'] = yhat_final.cpu().numpy()
    df['phi'] = phi_final.cpu().numpy()

out_path = "mochi_manual_fixed_output.tsv"
df.to_csv(out_path, sep='\t', index=False)
print("Saved predictions to", out_path)
