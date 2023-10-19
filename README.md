# masterthesis
Code for masterthesis "Unravelling regulatory interactions in stem cell differentiation in sponges", which explores stem cell differentiation using single-cell sequencing methods. 

## Aim 1: Generating an improved reference transcriptome for Suberites domuncula
* **transcriptome_rename_gtf.py**: script to rename reference transcriptome gtf file to more intuitive gene and transcript IDs

## Aim 2: Analysing single-cell gene expression data
* **sponge_scRNA.ipynb**: script to perform full scanpy analysis, including preprocessing, clustering and CytoTRACE pseudotime analysis
* **monocle3_conversion.R**: script to perform monocle3 analysis

## Aim 3: Investigating differentiation trajectories
* **monocle3_conversion.R**: the last code lines are used to select the differentiation trajectories, extract the cell barcodes and save them as an SVG file
* **sponge_diff_traj.ipynb**: reading in the SVG file of the cell barcodes, the adata object is subsetted to continue the analysis with only the current trajectory
* **full_trajectory_protocol**.ipynb: script which analyses the differentiation trajectory: Pearson correlation analysis, finding regulating motives, plotting motif abundance throughout a list of genes 


Packages:
Python                  3.10.1
R                       4.2.1
anndata                 0.7.5
biopython               1.80
conda                   23.3.1
h5py                    3.9.0
matplotlib              3.7.2
mudata                  0.2.3
muon                    0.1.5
numpy                   1.24.4
pandas                  2.0.3
pip                     23.0.1
scanpy                  1.9.3
scikit-learn            1.3.0
scipy                   1.11.1
seaborn                 0.12.2
umap-learn              0.5.3
wheel                   0.38.4
Seurat                  4.3.0
Monocle3                1.2.9
SeuratWrappers          0.3.0
SingleCellExperiment    1.20.1
zellkonverter           1.6.5
