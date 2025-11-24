import numpy as np
import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score

adata = sc.read_h5ad("sample_aorta_data_updated.h5ad")
X = np.load("aorta_genept_w.npy")
y = adata.obs["celltype"].astype(str).values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

knn = KNeighborsClassifier(n_neighbors=10, metric="cosine")
knn.fit(X_train, y_train)

y_pred = knn.predict(X_test)

print("Accuracy:", accuracy_score(y_test, y_pred))
print("F1:", f1_score(y_test, y_pred, average="weighted"))
