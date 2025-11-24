import numpy as np
import scanpy as sc

# Load cell data
adata = sc.read_h5ad("sample_aorta_data_updated.h5ad")

adata.var_names = [g.upper() for g in adata.var_names]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Load embeddings
gene_list = np.load("gene_list.npy", allow_pickle=True)
gene_embeddings = np.load("gene_embeddings_ada002.npy")

gene_to_idx = {g: i for i, g in enumerate(gene_list)}

# --- BUILD MAPPING FROM ADATA GENE ORDER → EMBEDDING ROWS ---
idxs = []
for g in adata.var_names:
    idxs.append(gene_to_idx.get(g, -1))

idxs = np.array(idxs)

# Filter genes with valid embeddings
valid = idxs >= 0
idxs = idxs[valid]

# Subset adata.X and embeddings
X_counts = adata.X[:, valid]           # (cells × valid_genes)
E = gene_embeddings[idxs, :]           # (valid_genes × 1536)

print("Shapes:", X_counts.shape, E.shape)

# --- FAST MATRIX MULTIPLICATION ---
X_gp = X_counts @ E                    # (cells × 1536)

# Normalize rows
norms = np.linalg.norm(X_gp, axis=1, keepdims=True)
X_gp = np.divide(X_gp, norms, where=norms != 0)

np.save("aorta_genept_w.npy", X_gp)
print("Saved aorta_genept_w.npy, shape:", X_gp.shape)
