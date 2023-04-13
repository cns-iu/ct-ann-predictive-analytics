import scanpy as sc, pandas as pd, logging, time, popv, argparse, numpy as np, requests, os
from pprint import pprint
from anndata import AnnData





def init_logging(log_filename, log_dir='logs'):
    """Starts a logger object as `log_dir/log_filename_todays_date.log`.

    Args:
        log_filename (str): General log-filename, ex- 'celltypist_preds'__todays_date.log
        log_dir (str, optional): Target folder for logging. Defaults to 'logs'.

    Returns:
        Logger
    """
    today = time.strftime("%Y%m%d%H", time.localtime(time.time()))
    print(today)
    logging.basicConfig(
        filename=f'{log_dir}/{log_filename}_{today}.log',
        level=logging.DEBUG,
        format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        force=True
    )
    logger = logging.getLogger()
    return logger
 
 
def parse_arguments(logger):
    """Generic function to parse script-arguments.

    Args:
        logger (str): Logger object for writing to `log_dir/log_filename_todays_date.log`.

    Returns:
        bool, str, str, str, str: MOUNT_GOOGLE_DRIVE, DATASET_NAME, QUERY_DATA_FILE, OUTPUT_PREDICTIONS_FILE, MODEL_NAME
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--TISSUE', default='Lung', type=str, help="Tabula Sapiens TISSUE to annotate your data with: ['Bladder','Blood','Bone_Marrow','Fat','Heart','Kidney','Large_Intestine','Liver','Lung','Lymph_Node','Mammary','Muscle','Pancreas','Prostate','Salivary Gland','Skin','Small_Intestine','Spleen','Thymus','Trachea','Vasculature']")
    parser.add_argument('--query-batch-key', default='donor_method', type=str, help='Key in query_adata.obs for batch correction. Set to None for no batch correction.')
    parser.add_argument('--algorithms', default='None', type=str, help='8 algorithm options: ["knn_on_scvi","scanvi","knn_on_bbknn","svm","rf","onclass","knn_on_scanorama","celltypist"]. By default uses all 8.')
    parser.add_argument('--query-labels-key', default='None', type=str, help='scANVI has the option to use labeled cells in the query dataset during training. To use some prelabeled cells from the query dataset, set query_labels_key to the corresponding key in query_adata.obs.')
    parser.add_argument('--unknown-celltype-label', default='None', type=str, help='If query_labels_key is not None, will treat everything not labeled unknown_celltype_label as a labeled cell.')
    parser.add_argument('--folder-name', const='Datasets', default='Datasets', type=str, nargs='?', help='Folder where input AnnData object resides. ex- Datasets/TS_Lung/TS_Lung.h5ad')
    parser.add_argument('--dataset-name', const='TS_Lung', default='TS_Lung', type=str, nargs='?', help='What dataset do you want to read? This should be the folder-name as well. ex- TS_Lung/TS_Lung.h5ad')
    parser.add_argument('--output-predictions-file', const='CellTypist_annotations.csv', default='CellTypist_annotations.csv', type=str, nargs='?', help='Name of the output predictions CSV file. Will be stored inside the folder with dataset-name. ex- TS_Lung/CellTypist_annotations.csv')
    
    args = parser.parse_args()
    logger.info(f'{args}=')
    TISSUE = args.TISSUE.lower()
    QUERY_BATCH_KEY = args.query_batch_key.lower()
    ALGORITHMS = [alg.lower() for alg in args.algorithms]
    QUERY_LABELS_KEY = args.query_labels_key.lower()
    UNKNOWN_CELLTYPE_LABEL = args.unknown_celltype_label.lower()
    
    ANNDATA_FOLDER = args.folder_name
    DATASET_NAME = args.dataset_name
    OUTPUT_PREDICTIONS_FILE = args.output_predictions_file
    
    QUERY_DATA_FILE = f'{DATASET_NAME}/{DATASET_NAME}.h5ad'

    return TISSUE, QUERY_BATCH_KEY, ALGORITHMS, QUERY_LABELS_KEY, UNKNOWN_CELLTYPE_LABEL, ANNDATA_FOLDER, DATASET_NAME, QUERY_DATA_FILE, OUTPUT_PREDICTIONS_FILE




start = time.time()
logger = init_logging(log_dir='logs', log_filename='cellTypist_preds')

TISSUE, QUERY_BATCH_KEY, ALGORITHMS, QUERY_LABELS_KEY, UNKNOWN_CELLTYPE_LABEL, ANNDATA_FOLDER, DATASET_NAME, QUERY_DATA_FILE, OUTPUT_PREDICTIONS_FILE = parse_arguments(logger)

ANNDATA_FOLDER = 'Datasets'

# import sys
# if 1==1:
#     sys.exit()


query_adata = sc.read_h5ad(f'{ANNDATA_FOLDER}/{QUERY_DATA_FILE}')

logger.info(f'var: {query_adata.var.columns.values}')
logger.info(f'uns: {query_adata.uns.keys()}')
logger.info(f'obs: {query_adata.obs.columns.values}')
logger.info(f'obsm: {query_adata.obsm.keys()}')



# Access the TISSUE-metadata hosted on Zenodo
TISSUE_url_metadata = requests.get("https://zenodo.org/api/records/7587774")
TISSUE_download_path = {
    ind["key"][3:-14]: ind["links"]["self"] for ind in TISSUE_url_metadata.json()["files"]
}
logger.info(TISSUE_download_path)



# `model` argument defaults to `Immune_All_Low.pkl`.
model_url_metadata = requests.get("https://zenodo.org/api/records/7580707")
pretrained_models_download_path = {
    ind["key"][18:-10]: ind["links"]["self"] for ind in model_url_metadata.json()["files"]
}

output_folder = ANNDATA_FOLDER
refdata_url = TISSUE_download_path[TISSUE]
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

output_fn = f"{output_folder}/TS_{TISSUE}.h5ad"
if not os.path.exists(output_fn):
    !wget -O $output_fn $refdata_url
    
logger.info(f'Downloaded PopV reference-training file from {refdata_url=}, at {output_fn=}')



model_url = pretrained_models_download_path[TISSUE]
output_model_tar_fn = f"{output_folder}/pretrained_model_{TISSUE}.tar.gz"
output_model_fn = f"{output_folder}/pretrained_model_{TISSUE}"
if not os.path.exists(output_model_fn):
    os.mkdir(output_model_fn)
    
if not os.path.exists(output_model_tar_fn):
    !wget -O $output_model_tar_fn $model_url
    !tar -xzf $output_model_tar_fn -C $output_model_fn


ref_adata = sc.read_h5ad(output_fn)

# Following parameters are specific to Tabula Sapiens dataset and contain the annotated cell-type and the batch_key that are corrected for during model training.
ref_labels_key = "cell_ontology_class"
ref_batch_key = "donor_assay"

min_celltype_size = np.min(ref_adata.obs.groupby(ref_labels_key).size())
n_samples_per_label = np.max((min_celltype_size, 500))


from popv.annotation import annotate_data
from popv.preprocessing import Process_Query

predictions_adata = Process_Query(
    query_adata,
    ref_adata,
    query_labels_key=QUERY_LABELS_KEY,
    query_batch_key=None,
    ref_labels_key=ref_labels_key,
    ref_batch_key=ref_batch_key,
    unknown_celltype_label=UNKNOWN_CELLTYPE_LABEL,
    save_path_trained_models=output_model_fn,
    cl_obo_folder="./PopV/ontology/",
    prediction_mode="inference",  # 'fast' mode gives fast results (does not include BBKNN and Scanorama and makes more inaccurate errors)
    n_samples_per_label=n_samples_per_label,
    use_gpu=False,
    compute_embedding=True,
    hvg=None,
).adata


# save_path = f"{output_folder}/popv_output"
annotate_data(predictions_adata, save_path=None)


# Write the predictions to a csv file
ABS_TGT_CSV_PATH = f'{ANNDATA_FOLDER}/{DATASET_NAME}/{OUTPUT_PREDICTIONS_FILE}'
logger.info(f'Writing the predictions anndata to {ABS_TGT_CSV_PATH}')
predictions_adata.obs.to_csv(ABS_TGT_CSV_PATH, index=False)


end = time.time()
logger.info(f'Finished executing PopV prediction pipeline in {end-start}!')