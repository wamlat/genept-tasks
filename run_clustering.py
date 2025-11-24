import numpy as np
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score

adata = sc.read_h5ad("sample_aorta_data_updated.h5ad")
X = np.load("aorta_genept_w.npy")

labels = adata.obs["celltype"].astype(str).values  # OR whatever the column is called
k = len(np.unique(labels))

km = KMeans(n_clusters=k, random_state=0).fit(X)
pred = km.labels_

print("ARI:", adjusted_rand_score(labels, pred))
print("AMI:", adjusted_mutual_info_score(labels, pred))
print("ASW:", silhouette_score(X, labels))