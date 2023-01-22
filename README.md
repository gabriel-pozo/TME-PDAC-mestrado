# TME-PDAC
Atlas level integration with scVI. To reproduce the notebook, first setup the data as it follows in seurat_prepare_atlas.
seurat_prepare_atlas --> Pre-process and basic QC filtering
scvi_pdac_integrate_ipynb --> Import data and setup anndata objects. Concatenate the objects and perform data integration with scVI.
