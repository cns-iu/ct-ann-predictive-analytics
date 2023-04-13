# Setting up FRMatch on Local system
if (!require("devtools",quietly=T))  install.packages("devtools")
if (!require("FRmatch",quietly=T))  devtools::install_github("JCVenterInstitute/FRmatch")


# Setting up SCDIOR on Local system: Required for h5ad -> SingleCellExperiment datatype conversion
if (!require("dior",quietly=T))  devtools::install_github('JiekaiLab/dior')
library(dior)
adata = dior::read_h5(file='../Datasets/TS_Lung/TS_Lung_predictions_anndata.h5ad', targ='anndata')



# Setting up ZellKonverter on Local system: Required for h5ad -> SingleCellExperiment datatype conversion
if (!require("BiocManager",quietly=T))  install.packages("BiocManager")
if (!require("zellkonverter",quietly=T))  BiocManager::install("zellkonverter")
if (!require("basilisk",quietly=T))  install.packages("basilisk")
if (!require("scRNAseq",quietly=T))  install.packages("scRNAseq")
library(zellkonverter)
library(basilisk)
library(scRNAseq)


adata <- zellkonverter::readH5AD(file='../Datasets/TS_Lung/TS_Lung_predictions_anndata.h5ad')
sce <- AnnData2SCE(adata)



# FR-Match
library(FRmatch)
library(SingleCellExperiment)
library(dplyr)
library(tibble)
library(data.table)



## read in pieces of input data - this may take a few minutes
DATASET_NAME <- "TS_Lung"
FOLDER <- paste0("Desktop/", DATASET_NAME, "/")
ANNOTATION_SOURCE_NAME <- "celltypist"
FILE_NAME_EXT <- "csv"
ANNOTATION_FILE_NAME <- paste0(FOLDER, ANNOTATION_SOURCE_NAME, "_preds.", FILE_NAME_EXT)

cell_cluster_labels <- fread(ANNOTATION_FILE_NAME,  sep=",") %>%
  # select(c(V1, predicted.ann_finest_level)) %>%
  mutate(V1=row_number()) %>%
  select(c(V1, majority_voting)) %>%
  # rename(c(Sample=V1, Cluster=predicted.ann_finest_level))
  rename(c(Sample=V1, Cluster=majority_voting))
cell_by_gene_expression <- fread(paste0(FOLDER, "/SingleCellExperiment_src_files/cell_by_gene_expression.csv")) %>%
  mutate(Sample=cell_cluster_labels$Sample) %>%
  select("Sample", everything())

NSForest_marker_genes <- fread(paste0(FOLDER, "NSForest_", ANNOTATION_SOURCE_NAME, "_markers.csv"))
# NSForest_fscores <- fread("NSForest_fscores.csv")
# MTG_taxonomy <- fread("MTG_taxonomy.csv")$x #need to be a vector


unique_markers <- unique(NSForest_marker_genes$markerGene)

# # Use the make_data_object() function to compile the pieces of input data into a data object.
# sce_dataset <- make_data_object(
#   dat = cell_by_gene_expression,
#   tab = cell_cluster_labels,
#   markers = unique_markers,
#   ## below are optional input data
#   cluster_marker_info = NSForest_marker_genes,
#   # f_score = NSForest_fscores,
#   # cluster_order = MTG_taxonomy
# )












###############################################################################
## rename key columns
names(cell_by_gene_expression)[1] <- "Sample"
names(cell_cluster_labels) <- c("Sample", "Cluster")

## replace special symbols by "_"
# cat("Replace any special symbol in sample and cluster names by '_'. \n")
# cell_by_gene_expression %<>% mutate(Sample=gsub("-| |\\.|/", "_", Sample))
# names(cell_by_gene_expression) <- gsub("-| |\\.|/", "_", names(cell_by_gene_expression))
# cell_cluster_labels %<>% mutate(Sample=gsub("-| |\\.|/", "_", Sample), Cluster=gsub("-| |\\.|/", "_", Cluster))
# unique_markers <- gsub("-| |\\.|/", "_", unique_markers)
# if(!is.null(NSForest_marker_genes)){
#   NSForest_marker_genes %<>% mutate(clusterName=gsub("-| |\\.|/", "_", clusterName), markerGene=gsub("-| |\\.|/", "_", markerGene))
# }
# if(!is.null(f_score)){
#   f_score %<>% mutate(clusterName=gsub("-| |\\.|/", "_", clusterName))
# }
# if(!is.null(cluster_order)){
#   cluster_order <- gsub("-| |\\.|/", "_", cluster_order)
# }

## cell_by_gene_expressiona cell_cluster_labelsle with "Sample", "Cluster", and gene columns for constructing the sce.object2
cell_by_gene_expression <- cell_by_gene_expression %>% inner_join(cell_cluster_labels, by="Sample") %>% #inner_join: make sure that cells are in the SAME order!!!
  select("Sample", "Cluster", everything())

##----------------------------------##
## make SingleCellExperiment object ## sce.object2
##----------------------------------##

## expression matrix
expr <- cell_by_gene_expression %>% column_to_rownames("Sample") %>% select(-Cluster) %>% t() #transpose: gene-by-cell
## row = genes
genenames <- rownames(expr) #gene names in order
## column = cells
sampnames <- colnames(expr)

## rowData: marker gene
df_marker_gene <- data.frame("marker_gene"=as.numeric(genenames %in% unique_markers), row.names=genenames)
## colData: cluster membership
df_cluster_membership <- cell_by_gene_expression %>% column_to_rownames("Sample") %>% select(cluster_membership=Cluster)

## make SingleCellExperiment object
sce.object2 <- SingleCellExperiment(assays = list(logcounts = expr), #logcounts!!!!!!!!!!!!!
                                   colData = df_cluster_membership,
                                   rowData = df_marker_gene)

## metacell_by_gene_expressiona:
if(!is.null(NSForest_marker_genes)){
  names(NSForest_marker_genes) <- c("clusterName", "markerGene", "scores")
  # NSForest_marker_genes %<>% mutate(clusterName=gsub("-| |\\.|/", "_", clusterName))
  # NSForest_marker_genes %<>% mutate(markerGene=gsub("-| |\\.|/", "_", markerGene))
}
if(!is.null(f_score)){
  names(f_score) <- c("clusterName", "score")
  # f_score %<>% mutate(clusterName=gsub("-| |\\.|/", "_", clusterName))
}
if(!is.null(cluster_order)){
  # cluster_order <- gsub("-| |\\.|/", "_", cluster_order)
}
metadata(sce.object2)$NSForest_marker_genes <- NSForest_marker_genes
metadata(sce.object2)$f_score <- f_score
metadata(sce.object2)$cluster_order <- cluster_order



rm(df_cluster_membership, df_marker_gene, cell_by_gene_expression, NSForest_marker_genes, genenames, sampnames, unique_markers, cell_cluster_labels, expr)
gc()


# AT2, CD4 T cells
plot_cluster_by_markers(sce.object2, cluster.name = "AT2", name.self = paste0(ANNOTATION_SOURCE_NAME," "))

# library(zellkonverter)
# library(basilisk)
# file <- tempfile(paste0(FOLDER, DATASET_NAME), fileext='.h5ad')
# writeH5AD( 
#   sce.object2, 
#   file,
#   X_name = NULL, 
#   skip_assays = FALSE, 
#   # compression = c("gzip") 
# )
