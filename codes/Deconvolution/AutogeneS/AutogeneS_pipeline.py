import autogenes as ag
import scanpy as as
import numpy as np
import pandas as pd

import sys
import os

scrna_path = sys.argv[1] //cells*genes
spatial_path = sys.argv[2] //locations*genes
celltype_key = sys.argv[3]
output_path = sys.argv[4]

sc_adata = sc.read_h5ad(scrna_path)
st_adata = sc.read_h5ad(spatial_path)

sc_adata.obs["cell_type"]=pd.read_csv(celltype_key)

intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)

# let us filter some genes
G = 2000
sc.pp.filter_genes(sc_adata, min_counts=10)

sc_adata.layers["counts"] = sc_adata.X.copy()

sc.pp.highly_variable_genes(
    sc_adata,
    n_top_genes=G,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)

sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata
# filter genes to be the same on the spatial data
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)
st_adata=st_adata.transpose()   //genes*locations
data_bulk_raw= pd.DataFrame(st_adata.X.toarray())
data_bulk_raw.index= st_adata.obs_names.to_list()
data_bulk_raw=data_bulk_raw.astype('float64')

#calculating the centroids of cell types
clusters = np.unique(sc_adata.obs['cell_type'].to_list())
sc_mean = pd.DataFrame(index=sc_adata.var_names,columns=clusters)
for cluster in clusters:
    
    cells =[]
    for x in range(len(sc_adata.obs_names)):
      if sc_adata.obs['cell_type'].values[x]==cluster:
        cells.append(sc_adata.obs['cell_type'].index[x])
    sc_part= sc_adata[cells,:].X.T
    sc_mean[cluster] = pd.DataFrame(np.mean(sc_part,axis=1),index=sc_adata.var_names)
    
centroids_sc_hv = sc_mean
centroids_sc_hv.shape

ag.init(centroids_sc_hv.T)
ag.optimize(ngen=2000,seed=0,nfeatures=200,mode='fixed',offspring_size=100,verbose=False)
coef_nusvr = ag.deconvolve(data_bulk_raw.T, model='nusvr')
coef_nnls = ag.deconvolve(data_bulk_raw.T, model='nnls')

pd.DataFrame(data=coef_nusvr,columns=clusters,index=data_bulk_raw.columns).to_csv(output_path+'/nusvr_result.csv')
pd.DataFrame(data=coef_nnls,columns=clusters,index=data_bulk_raw.columns).to_csv(output_path+'/nnls_result.csv')
