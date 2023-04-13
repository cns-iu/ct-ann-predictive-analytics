import scanpy as sc, pandas as pd, logging, time, celltypist, argparse, numpy as np
from pprint import pprint
from anndata import AnnData
from celltypist import models




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
    parser.add_argument('--mount-google-drive', default='False', type=str, help='Do you want to link to Google Drive? (True/False)')
    parser.add_argument('--existing-annotations-column', const='', default='', type=str, nargs='?', help='Column in the obs-dataframe of AnnData that contains annotations already. Defaults to blank. ex- cell_ontology_type')
    parser.add_argument('--folder-name', const='Datasets', default='Datasets', type=str, nargs='?', help='Folder where input AnnData object resides. ex- Datasets/TS_Lung/TS_Lung.h5ad')
    parser.add_argument('--dataset-name', const='TS_Lung', default='TS_Lung', type=str, nargs='?', help='What dataset do you want to read? This should be the folder-name as well. ex- TS_Lung/TS_Lung.h5ad')
    parser.add_argument('--output-predictions-file', const='CellTypist_annotations.csv', default='CellTypist_annotations.csv', type=str, nargs='?', help='Name of the output predictions CSV file. Will be stored inside the folder with dataset-name. ex- TS_Lung/CellTypist_annotations.csv')
    parser.add_argument('--model-name', const='Human_Lung_Atlas.pkl', default='Human_Lung_Atlas.pkl', type=str, nargs='?', help='CellTypist pretrained-model to use.')
    
    args = parser.parse_args()
    logger.info(f'{args}=')
    MOUNT_GOOGLE_DRIVE = args.mount_google_drive.lower() in ['true', 1, 't', 'yes']
    EXISTING_ANNOTATIONS_COLUMN = args.existing_annotations_column if args.existing_annotations_column!='' else None
    ANNDATA_FOLDER = args.folder_name
    DATASET_NAME = args.dataset_name
    OUTPUT_PREDICTIONS_FILE = args.output_predictions_file
    MODEL_NAME = args.model_name
    
    QUERY_DATA_FILE = f'{DATASET_NAME}/{DATASET_NAME}.h5ad'

    return MOUNT_GOOGLE_DRIVE, EXISTING_ANNOTATIONS_COLUMN, ANNDATA_FOLDER, DATASET_NAME, QUERY_DATA_FILE, OUTPUT_PREDICTIONS_FILE, MODEL_NAME




start = time.time()
logger = init_logging(log_dir='logs', log_filename='cellTypist_preds')

MOUNT_GOOGLE_DRIVE, EXISTING_ANNOTATIONS_COLUMN, ANNDATA_FOLDER, DATASET_NAME, QUERY_DATA_FILE, OUTPUT_PREDICTIONS_FILE, MODEL_NAME = parse_arguments(logger)
print(f'{MOUNT_GOOGLE_DRIVE=}, {EXISTING_ANNOTATIONS_COLUMN=}, {ANNDATA_FOLDER=}, {DATASET_NAME=}, {QUERY_DATA_FILE=}, {OUTPUT_PREDICTIONS_FILE=}, {MODEL_NAME=}')

if MOUNT_GOOGLE_DRIVE:
    from google.colab import drive
    drive.mount('/content/drive')
    ANNDATA_FOLDER = 'drive/MyDrive/CT-Annotation-tools-dataset'

# import sys
# if 1==1:
#     sys.exit()


query_adata = sc.read_h5ad(f'{ANNDATA_FOLDER}/{QUERY_DATA_FILE}')

logger.info(f'var: {query_adata.var.columns.values}')
logger.info(f'uns: {query_adata.uns.keys()}')
logger.info(f'obs: {query_adata.obs.columns.values}')
logger.info(f'obsm: {query_adata.obsm.keys()}')


# Enabling `force_update = True` will overwrite existing (old) models.
models.download_models(force_update = True)
logger.info(models.models_description())



# `model` argument defaults to `Immune_All_Low.pkl`.
model = models.Model.load(model = MODEL_NAME)
logger.info(f'model.cell_types.shape={model.cell_types.shape}, model.cell_types={model.cell_types}')



'''
Standardization of unclean datasets:
    1. `AnnData.var` object is not in standard format: `pd.DataFrame` with only indexes
    2. `AnnData.obs` object is not in standard format: `pd.DataFrame` with only 1 column of existing pre-annotated cell-labels.
    3. `AnnData.X` object contains duplicates: `pd.DataFrame`/`sparse-dataframe` contains duplicate columns of gene-expressions.
    4. `AnnData.var` object contains duplicates: `pd.DataFrame` contains duplicate values of genes-expressed.
    5. `AnnData.X` object is not normalized: The `CxG` is required to be `log1p` normalised expression to `10,000` counts per cell.
'''
clean_cxg_matrix = False

if 'feature_name' not in query_adata.var.columns:
    ALTERNATE_VARNAME = 'gene_symbol'
    logger.warning(f'feature_name not in query_adata.var.columns. Trying to use {ALTERNATE_VARNAME}')
    query_adata.var = query_adata.var.rename({ALTERNATE_VARNAME:'feature_name'}, axis=1)
if 'feature_name' in query_adata.var.columns and clean_cxg_matrix:
    query_adata.var = pd.DataFrame(index=query_adata.var['feature_name'])
    logger.warning(f'Somehow the original dataset was created with duplicates in the INDEX column itself!\nReset index and use this dirty column to filter out the duplicates. ')
    clean_CxG_matrix = query_adata.to_df().T.reset_index().drop_duplicates(subset=['feature_name'], keep='first').set_index('feature_name').T
    logger.info(f'clean_CxG_matrix.shape={clean_CxG_matrix.shape}')

    # Get rid of the dirty duplicated values
    query_adata.var.index = pd.Categorical(values=clean_CxG_matrix.columns.values, categories=set(clean_CxG_matrix.columns.values))

    # Create a whole new AnnData since Scanpy doesn't allow rewriting a new-shaped CxG matrix into existing X object.
    query_adata = AnnData(X=clean_CxG_matrix, var=query_adata.var, obs=query_adata.obs, uns=query_adata.uns, obsm=query_adata.obsm, obsp=query_adata.obsp) #, layers=query_adata.layers)
    logger.info(f'Created a whole new AnnData object of shape: {query_adata.shape}')

# Normalize the CxG matrix for CellTypist constraints
sc.pp.normalize_total(query_adata, target_sum = 1e4)
sc.pp.log1p(query_adata)

# The expression matrix (`adata_2000.X`) is pre-processed (and required) as log1p normalised expression to 10,000 counts per cell (this matrix can be alternatively stashed in `.raw.X`).
logger.warning(f'The CxG df should be log1p normalized to 10,000 counts per cell:\n {np.expm1(query_adata.X).sum(axis = 1)}')


# Make the var names unique now. Hopefully it works
query_adata.var_names_make_unique()


# Ignore; predict cell identities using this loaded model.
# predictions = celltypist.annotate(adata_2000, model = model, majority_voting = True)

# Alternatively, just specify the model name (recommended as this ensures the model is intact every time it is loaded).
predictions = celltypist.annotate(query_adata, model = MODEL_NAME, majority_voting = True)

# Get an `AnnData` with predicted labels embedded into the cell metadata columns.
predictions_adata = predictions.to_adata()

# Write the predictions anndata-object to a csv file
if EXISTING_ANNOTATIONS_COLUMN:
    predictions_adata.obs['exp_vs_pred'] = predictions_adata.obs.apply(lambda row: row[EXISTING_ANNOTATIONS_COLUMN]==row['majority_voting'], axis=1)
    logger.info(f'{predictions_adata.obs["exp_vs_pred"].value_counts()=}')

ABS_TGT_CSV_PATH = f'{ANNDATA_FOLDER}/{DATASET_NAME}/{OUTPUT_PREDICTIONS_FILE}'
logger.info(f'Writing the predictions anndata to {ABS_TGT_CSV_PATH}')
predictions_adata.obs.to_csv(ABS_TGT_CSV_PATH, index=False)


end = time.time()
logger.info(f'Finished executing CellTypist prediction pipeline in {end-start}!')