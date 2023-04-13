# if (!requireNamespace("satijalab/seurat-data", quietly = TRUE)) {
#   remove.packages(grep("spatstat", installed.packages(), value = T))
#   devtools::install_version("spatstat", version = "1.64-1")
#   devtools::install_github("satijalab/seurat-data")
# }

# if (!requireNamespace("satijalab/azimuth", quietly = TRUE)) {
#   devtools::install_github("satijalab/azimuth")
# }

# if (!requireNamespace("Seurat", quietly = TRUE)) {
#   install.packages("Seurat")
# }

# if (!requireNamespace("logr", quietly = TRUE)) {
#   install.packages("logr")
# }




rm(list = ls(all.names = TRUE))

library(Seurat)
print('loaded Seurat')
library(Azimuth)
print('loaded Azimuth')
library(SeuratData)
print('loaded SeuratData')

library(patchwork)
print('loaded patchwork')

library(logr)
print('loaded logr')


LOGS_DIRECTORY <- 'logs'
CURR_TIMESTAMP <- format(Sys.time(), '%Y%M%d%H')
LOG_FILENAME <- paste0('azimuth_preds_', CURR_TIMESTAMP, '.log')
print(LOG_FILENAME)
log_open(paste0(LOGS_DIRECTORY, '/', LOG_FILENAME))


print_logs <- function(msg){
  log_print(paste0(Sys.time(), ' : ', msg))
}

start <- Sys.time()




# METHOD 1: DIRECTLY USE THE RunAzimuth(...) function

# Remove all objects in current R Workspace
args <- commandArgs(trailingOnly = TRUE)
args <- as.list(strsplit(args, " ")[[1]])
print_logs(args)


# Initialize all args
print_logs('Initializing arguments')
ANNDATA_FOLDER <- args[1][[1]] # "Datasets"
QUERY_DATASET_NAME <- args[2][[1]] # "TS_Lung"
OUTPUT_PREDICTIONS_FILE <- args[3][[1]] # "azimuth_preds.tsv"
REFERENCE <- args[4][1] # "lungref"
ABS_DATASET_PATH <- paste0(ANNDATA_FOLDER, "/", QUERY_DATASET_NAME, "/", QUERY_DATASET_NAME,".h5ad")


print_logs(getwd())
print_logs(paste0("Loading the query dataset : ", ABS_DATASET_PATH))
print(paste0("****  Loading the query dataset : ", ABS_DATASET_PATH))
query_adata <- LoadFileInput(path = ABS_DATASET_PATH)
print_logs("Loaded the query dataset")
options(timeout=360)


lung_results <- RunAzimuth(query_adata, reference=REFERENCE)
print_logs(paste0("Running Azimuth using the reference [", REFERENCE, "] dataset."))
# Fails at: ??InstallData(reference)


azimuth_preds <- data.frame(
  predicted.ann_level_1 = lung_results$predicted.ann_level_1,
  predicted.ann_level_1.score = lung_results$predicted.ann_level_1.score,
  predicted.ann_level_2 = lung_results$predicted.ann_level_2,
  predicted.ann_level_2.score = lung_results$predicted.ann_level_2.score,
  predicted.ann_level_3 = lung_results$predicted.ann_level_3,
  predicted.ann_level_3.score = lung_results$predicted.ann_level_3.score,
  predicted.ann_level_4 = lung_results$predicted.ann_level_4,
  predicted.ann_level_4.score = lung_results$predicted.ann_level_4.score,
  predicted.ann_level_5 = lung_results$predicted.ann_level_5,
  predicted.ann_level_5.score = lung_results$predicted.ann_level_5.score,
  predicted.ann_finest_level = lung_results$predicted.ann_finest_level,
  predicted.ann_finest_level.score = lung_results$predicted.ann_finest_level.score
)

tgt_filename <- paste0(ANNDATA_FOLDER, "/", QUERY_DATASET_NAME, "/", OUTPUT_PREDICTIONS_FILE)
print_logs(paste0("Writing the azimuth predictions into ", tgt_filename))
write.table(x=azimuth_preds, file=tgt_filename)
print_logs(paste0("Successfully wrote the azimuth predictions into ", tgt_filename))

end=Sys.time()
cat("\n")
print_logs(end-start)







# METHOD 2: USE CUSTOM REFERENCE-FILES TO GET PREDICTIONS
#!/usr/bin/env Rscript

# # Ensure Seurat v4.0 or higher is installed
# if (packageVersion(pkg = "Seurat") < package_version(x = "4.0.0")) {
#   stop("Mapping datasets requires Seurat v4 or higher.", call. = FALSE)
# }

# # Ensure glmGamPoi is installed
# if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
#   if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     BiocManager::install("glmGamPoi")
#   }
# }

# # Ensure Azimuth is installed
# if (packageVersion(pkg = "Azimuth") < package_version(x = "0.4.0")) {
#   stop("Please install azimuth - remotes::install_github('satijalab/azimuth')", call. = FALSE)
# }

# library(Seurat)
# library(Azimuth)

# # Download the Azimuth reference and extract the archive

# # Load the reference
# # Change the file path based on where the reference is located on your system.
# reference <- LoadReference(path = "Reference_directory_hlca/")

# # Load the query object for mapping
# # Change the file path based on where the query file is located on your system.
# query <- LoadFileInput(path = ABS_DATASET_PATH)
# query <- ConvertGeneNames(
#   object = query,
#   reference.names = rownames(x = reference$map),
#   homolog.table = 'https://seurat.nygenome.org/azimuth/references/homologs.rds'
# )

# # Calculate nCount_RNA and nFeature_RNA if the query does not
# # contain them already
# if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
#     calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
#     colnames(x = calcn) <- paste(
#       colnames(x = calcn),
#       "RNA",
#       sep = '_'
#     )
#     query <- AddMetaData(
#       object = query,
#       metadata = calcn
#     )
#     rm(calcn)
# }

# # Calculate percent mitochondrial genes if the query contains genes
# # matching the regular expression "^MT-"
# if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
#   query <- PercentageFeatureSet(
#     object = query,
#     pattern = '^MT-',
#     col.name = 'percent.mt',
#     assay = "RNA"
#   )
# }

# # Filter cells based on the thresholds for nCount_RNA and nFeature_RNA
# # you set in the app
# cells.use <- query[["nCount_RNA", drop = TRUE]] <= 18063 &
#   query[["nCount_RNA", drop = TRUE]] >= 567 &
#   query[["nFeature_RNA", drop = TRUE]] <= 3992 &
#   query[["nFeature_RNA", drop = TRUE]] >= 497

# # If the query contains mitochondrial genes, filter cells based on the
# # thresholds for percent.mt you set in the app
# if ("percent.mt" %in% c(colnames(x = query[[]]))) {
#   cells.use <- cells.use & (query[["percent.mt", drop = TRUE]] <= 99 &
#     query[["percent.mt", drop = TRUE]] >= 0)
# }

# # Remove filtered cells from the query
# query <- query[, cells.use]

# # Preprocess with SCTransform
# query <- SCTransform(
#   object = query,
#   assay = "RNA",
#   new.assay.name = "refAssay",
#   residual.features = rownames(x = reference$map),
#   reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
#   method = 'glmGamPoi',
#   ncells = 2000,
#   n_genes = 2000,
#   do.correct.umi = FALSE,
#   do.scale = FALSE,
#   do.center = TRUE
# )

# # Find anchors between query and reference (##### TRAINING STEP #####)
# anchors <- FindTransferAnchors(
#   reference = reference$map,
#   query = query,
#   k.filter = NA,
#   reference.neighbors = "refdr.annoy.neighbors",
#   reference.assay = "refAssay",
#   query.assay = "refAssay",
#   reference.reduction = "refDR",
#   normalization.method = "SCT",
#   features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
#   dims = 1:50,
#   n.trees = 20,
#   mapping.score.k = 100
# )

# # Transfer cell type labels and impute protein expression (##### PREDICTION STEP #####)
# #
# # Transferred labels are in metadata columns named "predicted.*"
# # The maximum prediction score is in a metadata column named "predicted.*.score"
# # The prediction scores for each class are in an assay named "prediction.score.*"
# # The imputed assay is named "impADT" if computed

# refdata <- lapply(X = "ann_finest_level", function(x) {
#   reference$map[[x, drop = TRUE]]
# })

# names(x = refdata) <- "ann_finest_level"
# if (FALSE) {
#   refdata[["impADT"]] <- GetAssayData(
#     object = reference$map[['ADT']],
#     slot = 'data'
#   )
# }

# query <- TransferData(
#   reference = reference$map,
#   query = query,
#   dims = 1:50,
#   anchorset = anchors,
#   refdata = refdata,
#   n.trees = 20,
#   store.weights = TRUE
# )

# # Calculate the embeddings of the query data on the reference SPCA
# query <- IntegrateEmbeddings(
#   anchorset = anchors,
#   reference = reference$map,
#   query = query,
#   reductions = "pcaproject",
#   reuse.weights.matrix = TRUE
# )

# # Calculate the query neighbors in the reference
# # with respect to the integrated embeddings
# query[["query_ref.nn"]] <- FindNeighbors(
#   object = Embeddings(reference$map[["refDR"]]),
#   query = Embeddings(query[["integrated_dr"]]),
#   return.neighbor = TRUE,
#   l2.norm = TRUE
# )

# # The reference used in the app is downsampled compared to the reference on which
# # the UMAP model was computed. This step, using the helper function NNTransform,
# # corrects the Neighbors to account for the downsampling.
# query <- NNTransform(
#   object = query,
#   meta.data = reference$map[[]]
# )

# # Project the query to the reference UMAP.
# query[["proj.umap"]] <- RunUMAP(
#   object = query[["query_ref.nn"]],
#   reduction.model = reference$map[["refUMAP"]],
#   reduction.key = 'UMAP_'
# )


# # Calculate mapping score and add to metadata
# query <- AddMetaData(
#   object = query,
#   metadata = MappingScore(anchors = anchors),
#   col.name = "mapping.score"
# )

# # VISUALIZATIONS

# # First predicted metadata field, change to visualize other predicted metadata
# id <- "ann_finest_level"[1]
# predicted.id <- paste0("predicted.", id)

# # DimPlot of the reference
# DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()

# # DimPlot of the query, colored by predicted cell type
# DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id, label = TRUE) + NoLegend()

# # Plot the score for the predicted cell type of the query
# FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "proj.umap")
# VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()

# # Plot the mapping score
# FeaturePlot(object = query, features = "mapping.score", reduction = "proj.umap")
# VlnPlot(object = query, features = "mapping.score", group.by = predicted.id) + NoLegend()

# # Plot the prediction score for the class CD16 Mono
# FeaturePlot(object = query, features = "CD16 Mono", reduction = "proj.umap")
# VlnPlot(object = query, features = "CD16 Mono", group.by = predicted.id) + NoLegend()

# # Plot an RNA feature
# FeaturePlot(object = query, features = "EMP2", reduction = "proj.umap")
# VlnPlot(object = query, features = "EMP2", group.by = predicted.id) + NoLegend()

# # Plot an imputed protein feature
# if (FALSE) {
#   FeaturePlot(object = query, features = "CD3-1", reduction = "proj.umap")
#   VlnPlot(object = query, features = "CD3-1", group.by = predicted.id) + NoLegend()
# }



# Close the log file
log_close()