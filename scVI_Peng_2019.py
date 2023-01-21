#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Load the libraries

import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scanpy.external as sce
import seaborn as sns
import anndata as ad
import squidpy as sq
import scvi
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import warnings
warnings.filterwarnings("ignore")
from sccoda.util import comp_ana as mod


# In[2]:


# Load PRJCA001063

adata = sc.read("D:Scanpy/peng_pdac_besca2.raw.h5ad")


# In[3]:


adata.obs.rename(columns = {'CONDITION':'Condition'}, inplace = True)
adata.obs


# In[4]:


del(adata.obs['CELL'])
del(adata.obs['Type'])
del(adata.obs['Cell_type'])
adata.obs


# In[5]:


adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[7]:


adata.obs


# In[8]:


adata.layers['counts'] = adata.X.copy()


# In[9]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata


# In[11]:


sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=3000, layer='counts', batch_key='Patient',subset=True)


# In[12]:


adata


# In[13]:


scvi.model.SCVI.setup_anndata(adata, layer = "counts",
                             categorical_covariate_keys=["Patient"],
                             continuous_covariate_keys=['total_counts'])


# In[14]:


model = scvi.model.SCVI(adata)


# In[15]:


model.train()


# In[16]:


adata.obsm['X_scVI'] = model.get_latent_representation()


# In[17]:


adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)


# In[18]:


sc.pp.neighbors(adata, use_rep = 'X_scVI')


# In[19]:


sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.3)


# In[20]:


sc.pl.umap(adata, color = ['leiden', 'Condition', 'Patient'], frameon = False, wspace=0.3)


# In[24]:


#sc.pl.umap(adata, color=['leiden','Condition','FXYD3','COL1A1','PTPRC','PLVAP'],ncols=2,cmap='viridis',frameon=False,wspace=0.2)


# In[25]:


sc.pl.umap(adata, color='leiden', legend_loc='on data', frameon=False)


# In[32]:


#Epithelial markers

sc.pl.umap(adata, color=['AMY2A','CFTR','FXYD3','INS'],cmap='viridis',frameon=False,ncols=2)


# In[33]:


#Immune markers

sc.pl.umap(adata, color=['CD3D','CD68','MS4A1','MZB1'],cmap='viridis',frameon=False,ncols=2)


# In[35]:


#Stromal markers

sc.pl.umap(adata, color=['COL1A1','RGS5','PLVAP','ACTA2'],cmap='viridis',ncols=2,frameon=False)


# In[38]:


sc.pl.umap(adata, color='Patient',frameon=False)


# In[63]:


cell_type = {"0":"Endothelial cells",
"1":"Ductal cells",
"2":"Fibroblasts",
"3":"Stellate cells",
"4":"Myeloid cells",
"5":"Malignant cells",
"6":"T-NK cells",
"7":"Malignant cells",
"8":"Ductal cells",
"9":"Acinar cells",
"10":"B cells",
"11":"Malignant cells",
"12":"Malignant cells",
"13":"Malignant cells",
"14":"Endocrine cells",
"15":"B cells",
"16":"Plasmablasts",
"17":"Endothelial cells"
}


# In[64]:


adata.obs['cell type'] = adata.obs.leiden.map(cell_type)


# In[65]:


sc.pl.umap(adata, color = ['cell type'], frameon = False)


# In[66]:


adata.obs.groupby(['Patient']).count()


# In[67]:


num_tot_cells = adata.obs.groupby(['Patient']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.leiden))
num_tot_cells


# In[68]:


cell_type_counts = adata.obs.groupby(['Patient', 'Condition', 'cell type']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]
cell_type_counts


# In[69]:


cell_type_counts['total_cells'] = cell_type_counts.Patient.map(num_tot_cells).astype(int)

cell_type_counts['frequency'] = cell_type_counts.n_genes_by_counts / cell_type_counts.total_cells

cell_type_counts


# In[70]:


import matplotlib.pyplot as plt

plt.figure(figsize = (10,4))

ax = sns.boxplot(data = cell_type_counts, x = 'cell type', y = 'frequency', hue = 'Condition')

plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')

plt.show()


# In[60]:


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')


# In[61]:


markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers


# In[62]:


markers_scvi = model.differential_expression(groupby = 'leiden')


# In[71]:


adata.uns['scvi_markers'] = markers_scvi
adata.uns['markers'] = markers


# In[72]:


adata.write_h5ad('D:/Scanpy/scvipeng.h5ad')


# In[73]:


model.save('D:/Scanpy/pengscvimodel.model')

