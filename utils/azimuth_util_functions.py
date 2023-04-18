import pandas as pd

def get_specific_azimuth_marker_genes(filename='hlca_level1.csv', verbose=False):
    """Accesses the Github repo for Azimuth, and pulls the lung-information from specific URLs.
    
    Args:
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.
    
    Return:
        pd.DataFrame
    """
    USER = 'satijalab'
    REPOSITORY_NAME = 'azimuth_website'
    FILE_PATH = f'/static/csv/{filename}'
    URL = f'https://raw.githubusercontent.com/{USER}/{REPOSITORY_NAME}/master/{FILE_PATH}'
    df = pd.read_csv(URL, index_col=0)
    if verbose:  print(f'data = \n{df}')
    return df




def get_all_azimuth_marker_genes(organ='lung', verbose=False):
    """Accesses the Github repo for Azimuth, and pulls the lung-information from below URLs.
    `https://github.com/satijalab/azimuth_website/blob/master/static/csv/hlca_level1.csv`,
    `https://github.com/satijalab/azimuth_website/blob/master/static/csv/hlca_level2.csv`,
    `https://github.com/satijalab/azimuth_website/blob/master/static/csv/hlca_level3.csv`,
    `https://github.com/satijalab/azimuth_website/blob/master/static/csv/hlca_level4.csv`,
    `https://github.com/satijalab/azimuth_website/blob/master/static/csv/hlca_level5.csv`,
    `https://github.com/satijalab/azimuth_website/blob/master/static/csv/hlca_finest_level.csv`

    There are also some naming-convention issues to be fixed:
        The Azimuth HLCA_finest marker-genes file, contains "non classical monocytes" and "transitional club at2" while the 
        actual algorithm predictions are "non-classical monocytes" and "transitional club-at2" respectively.

    Args:
        organ (str, optional): Defaults to 'lung'.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        pd.DataFrame: Columns are [`'Label'`, `'Markers'`].
    """
    azimuth_marker_genes_df = pd.DataFrame()
    if organ=='lung':
        FILENAMES = ['hlca_level1', 'hlca_level2', 'hlca_level3', 'hlca_level4', 'hlca_level5', 'hlca_finest_level']
    for filename in FILENAMES:
        if verbose:  print(f'Trying to access {filename} in Azimuth Github-Repo.')
        curr_marker_genes_df = get_specific_azimuth_marker_genes(filename=f'{filename}.csv', verbose=verbose)
        curr_marker_genes_df = curr_marker_genes_df[['Label','Markers']]
        azimuth_marker_genes_df = pd.concat([azimuth_marker_genes_df, curr_marker_genes_df]).drop_duplicates(subset=['Label'])
    azimuth_marker_genes_df.columns = ['unique_cts', 'azimuth_markers']
    azimuth_marker_genes_df['unique_cts'] = azimuth_marker_genes_df['unique_cts'].str.lower().replace('φ','ï†')
    
    # Fixing some clearly naming-convention issues.
    azimuth_marker_genes_df.loc[azimuth_marker_genes_df['unique_cts']=='non classical monocytes' , 'unique_cts'] = 'non-classical monocytes'
    azimuth_marker_genes_df.loc[azimuth_marker_genes_df['unique_cts']=='transitional club at2' , 'unique_cts'] = 'transitional club-at2'
    azimuth_marker_genes_df = azimuth_marker_genes_df.reset_index(drop=True)
    return azimuth_marker_genes_df



def filter_valid_set_genes(input_genes, genes_available_in_anndata, verbose=False):
    """Keep only those marker-genes that exist within the input AnnData.
    After this standardization, we can utilize the scanpy.stacked_violin_plot() function.

    Args:
        input_genes (str): The data fetched from Azimuth's Github for Marker-genes, contains [`marker-genes`] column as a string.
        genes_available_in_anndata (iterable): You can access the unique gene-names present within the input AnnData as `predictions_adata.var_names`.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        list: List of valid marker-genes that exist within the input-AnnData.
    """
    input_genes = [gene.strip() for gene in input_genes.split(',')]
    valid_genes = set(genes_available_in_anndata) & set(input_genes)
    valid_genes = list(valid_genes)
    if verbose:  print(valid_genes)
    return valid_genes
