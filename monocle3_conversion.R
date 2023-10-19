#script to convert anndata h5ad file from scanpy to a seurat and later cds file for monocle3
Sys.setenv(LANG = "en")
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(reticulate)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(zellkonverter)

#read in h5ad file from scanpy and create seurat object
ad <- readH5AD("/Users/carolinaatria/Desktop/ADesktop/Studium/Master/master_thesis/sponge/data_from_roger/python/PART1/isoseq_data/subdom.h5ad")
adata_Seurat <- as.Seurat(ad, counts = "X", data = NULL)

#due to version issues, fix UMAP annotation
DimPlot(adata_Seurat,reduction="X_umap")
adata_Seurat[["UMAP"]] <- adata_Seurat[["X_umap"]]
DimPlot(adata_Seurat,reduction="UMAP")

#convert to cds for monocle3
cds <- SeuratWrappers::as.cell_data_set(adata_Seurat)
#cds <- preprocess_cds(cds,method="PCA")
#cds <- reduce_dimension(cds,max_components=2,reduction_method = "UMAP")

cds <- cluster_cells(cds = cds, reduction_method = "UMAP" ,resolution = 0.2)
cds <- learn_graph(cds)

#choose root cells manually or by using the corresponding node
root_set <- choose_cells(cds)
root = (root_set@colData@rownames)

cds <- order_cells(cds,root_cells=root)
plot_cells(cds,
           color_cells_by = "ct_pseudotime", #"ct_pseudotime", CYTOtrace
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=T,
           label_roots = F,
           label_principal_points = FALSE,
           trajectory_graph_color = 'cyan',
           trajectory_graph_segment_size = 1,
           cell_size = 1,
           cell_stroke = 0,)

#FOR aim 3, to select a differentiation trajectory!
stem_cell = colnames(cds)[cds@colData[["ct_pseudotime"]] == 0.0]
cds_sub <- choose_graph_segments(cds,starting_pr_node = 34,ending_pr_nodes = 3)
cds_sub <- choose_graph_segments(cds)


#create new seurat object for subset
traj1_seurat <- subset(x=adata_Seurat,cells=cds_sub@colData@rownames)
list_of_cells <- traj1_seurat@assays[["originalexp"]]@counts@Dimnames[[2]]

# Export the list or dataframe to a CSV file
write.csv(list_of_cells, "new_cluster_trajectory.csv", row.names = FALSE)
