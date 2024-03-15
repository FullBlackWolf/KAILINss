import scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import accuracy_score

adata = sc.read_h5ad('F:\\CJP_cellfate\\data\\reprogramming.h5ad')
adata = adata[adata.obs.label.isin(['failed_prog','reprogrammed_prog'])]
mtx_use = adata.X.toarray()
label = adata.obs.label
label = label.replace('failed_prog',0)
label = label.replace('reprogrammed_prog',1)

# Kmeans
kmeans_model = KMeans(n_clusters=2)
kmeans_model.fit(mtx_use)
kmeans_labels = kmeans_model.labels_

# hclust
hclust_model = AgglomerativeClustering(n_clusters=2)
hclust_labels = hclust_model.fit_predict(mtx_use)

# acc
kmeans_accuracy = accuracy_score(label, kmeans_labels)
hclust_accuracy = accuracy_score(label, hclust_labels)

print('------reprogramming data------')
print("kmeans acc:", kmeans_accuracy)
print("hclust acc:", hclust_accuracy)


adata = sc.read_h5ad('E:\\public\larry\\larry_mono_neu.h5ad')
adata = adata[adata.obs.label.isin(['Neutrophil_prog','Monocyte_prog'])]
mtx_use = adata.X.toarray()
label = adata.obs.label
label = label.replace('Neutrophil_prog',0)
label = label.replace('Monocyte_prog',1)

# Kmeans
kmeans_model = KMeans(n_clusters=2)
kmeans_model.fit(mtx_use)
kmeans_labels = kmeans_model.labels_

# hclust
hclust_model = AgglomerativeClustering(n_clusters=2)
hclust_labels = hclust_model.fit_predict(mtx_use)

# acc
kmeans_accuracy = accuracy_score(label, kmeans_labels)
hclust_accuracy = accuracy_score(label, hclust_labels)

print('------hema data------')
print("kmeans acc:", kmeans_accuracy)
print("hclust acc:", hclust_accuracy)