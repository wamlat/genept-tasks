# GenePT Benchmarks – Aorta Dataset (Cell Type & Phenotype)

This repository contains a minimal, fully-working implementation of the **GenePT-w** downstream evaluation pipeline, reproduced from  
**Chen & Zou (2024), “GenePT: A Simple But Effective Foundation Model for Genes and Cells Built From ChatGPT.”**

We benchmark GenePT-w embeddings on the **Aorta** single-cell dataset for two tasks:
- **Cell type classification** (supervised kNN)
- **Cell type clustering** (unsupervised)
- (Optional) **Phenotype clustering** using a different column

All steps match the methodology described in the GenePT paper.

---
Embeddings (GenePT_gene_embedding_ada_text.pickle) can be downloaded here: https://zenodo.org/records/10833191

Aorta data (sample_aorta_data_updated.h5ad) can be downloaded here: https://drive.google.com/drive/folders/1LgFvJqWNq9BqHbuxB2tYf62kXs9KqL4t

## Repository Structure

GENEPT-TASKS/
│
├── sample_aorta_data_updated.h5ad # Aorta AnnData from GenePT authors
├── GenePT_gene_embedding_ada_text.pickle
│ # dict: gene → 1536-d embedding (ada-002)
├── gene_list.npy # extracted gene names
├── gene_embeddings_ada002.npy # extracted embedding matrix
│
├── build_genept.py # vectorized GenePT-w embedding builder
├── embedding_matrix.py # converts pickle → npy matrix
├── run_clustering.py # unsupervised ARI / AMI / ASW
├── run_knn.py # supervised accuracy / F1
│
└── aorta_genept_w.npy # generated GenePT-w cell embeddings



## Quickstart

1. Build GenePT-w embeddings
python build_genept.py
This loads the Aorta .h5ad, normalizes it, maps gene names to GenePT embeddings, and computes:
```python
aorta_genept_w.npy   # shape: (cells × 1536)
```

2. Run unsupervised clustering
```python
python run_clustering.py
```
Outputs ARI / AMI / ASW.

3. Run supervised kNN classification
```python
python run_knn.py
```
Outputs accuracy / weighted F1.

### Notes
The embedding file contains ≈93k gene entries; only overlapping genes are used.
Matrix multiplication is used for fast embedding computation.
The metrics align with the expected GenePT-w cell type performance from the paper.

### Citations
Chen, Yiqun & Zou, James (2024).
GenePT: A Simple But Effective Foundation Model for Genes and Cells Built From ChatGPT.

