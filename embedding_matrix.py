import pickle
import numpy as np

with open("GenePT_gene_embedding_ada_text.pickle", "rb") as f:
    emb_dict = pickle.load(f)

gene_list = list(emb_dict.keys())
dim = len(next(iter(emb_dict.values())))

print("Number of genes:", len(gene_list))
print("Embedding dimension:", dim)

# Convert dict-of-lists → matrix
emb_matrix = np.zeros((len(gene_list), dim), dtype=float)

for i, g in enumerate(gene_list):
    emb_matrix[i] = np.array(emb_dict[g], dtype=float)

# Save for future use
np.save("gene_list.npy", np.array(gene_list, dtype=object))
np.save("gene_embeddings_ada002.npy", emb_matrix)

print("Saved gene_list.npy and gene_embeddings_ada002.npy")