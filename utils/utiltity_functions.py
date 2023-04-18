import  pandas as pd
from subprocess import Popen, PIPE
from asctb_ct_label_mapper.utilities.nlp_preprocessing import execute_nlp_pipeline



def run_external_script(script_command, script_name, args):
    """Invokes the R-code from Python using 32-bit or Rscript 3.4.4 command.
        Uses the Python subprocess module to create a new Pipe.
        Reusable function to run any external Python/R script with arguments.

    Args:
        script_command (str): Rscript/python
        script_name (str): Script-to-invoke
        args (str): Args-for-script
    """
    try:
        cmd = [script_command, script_name, args]
        pipe = Popen( cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE )
        output, error = pipe.communicate()

        script_command = script_command.capitalize()
        if pipe.returncode == 0:
            print(f'{script_command} OUTPUT:\n')
            print(output.decode())
        else:
            print(f'{script_command} OUTPUT:\n')
            print(output.decode())
            print(f'{script_command} ERROR:\n')
            print(error.decode())
    except Exception as e:
        print(f'\nSomething went wrong in the {script_command}-Script {script_name}.')
        print(e)
    print('If the script fails due to usage issues, try the following command directly in a terminal:')
    print(script_command, script_name, args)






def get_unq_vals_as_str(x, verbose=False):
    """Fetches the unique values in the input iterable. Returns a string with unique values on each line. 

    Args:
        x (list/set/np.array): Iterable object containing values.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        str: String with unique values on each line.
    """
    res = '<br>'
    try:
        allvals = sorted(list(set(x)))
        if verbose:  print(allvals)
        res += '<br>'.join([str(v) for v in allvals])
    except Exception as e:
        if verbose:  print(f'Error: {e}')
        res += '## Error parsing values ##'
    return res




def clean_and_translate_annotation(input_label):
    """Returns the cleaned version of the annotation label after performing the following steps:
    
    ```python
    remove_whitespaces()
    expand_word_contractions()
    replace_special_chars()
    convert_number_to_word()
    make_lowercase()
    get_root_word()
    ```

    Args:
        input_label (str): Input annotation label text.

    Returns:
        str: Cleaned version of the annotation label text.
    """
    return ' '.join([execute_nlp_pipeline(word) for word in input_label.split()])




def create_summary_obs(adata, max_unique=50, verbose=False):
    """Create a summary dataframe containing all unique sorted values along with counts for the `obs` dataframe in the AnnData object.

    Args:
        adata (AnnData): Input AnnData object containing an `obs` dataframe that we can summarize and view.
        max_unique (int, optional): Threshold upto which we show all values on the interactive plot. Defaults to 50.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        pd.DataFrame: Summarized `obs` dataframe.
    """
    if verbose:  print(f'Getting unique values per columns in the obs-dataframe.')
    unique_vals_per_col_df = pd.DataFrame(adata.obs.nunique()).reset_index()
    unique_vals_per_col_df.columns = ['colname','nunique']

    if verbose:  print(f'Filtering out columns that have more than {max_unique} values.')
    cols_to_view_in_detail = unique_vals_per_col_df.loc[unique_vals_per_col_df['nunique'] < max_unique , 'colname'].tolist()

    if verbose:  print(f'Pulling all values for columns that have less than {max_unique} values.')
    details_subset_df = pd.DataFrame(adata.obs.loc[:, cols_to_view_in_detail].apply( get_unq_vals_as_str ).T).reset_index()
    details_subset_df.columns = ['colname', 'unique']

    if verbose:  print(f'Creating the final merged dataframe containing all unique sorted values along with counts.')
    summary_obs_df = pd.merge(left=unique_vals_per_col_df, right=details_subset_df, how='left', on='colname')

    return summary_obs_df





def get_aggregated_cellbygene_matrix(adata, celltype_labels_column='azimuth_preds', aggregation_method='mean', verbose=False):
    """Ensures a standardized CxG matrix in the input AnnData.
    Performs required CxG aggregations by grouping for each Cell-Type label.

    Args:
        adata (Scanpy.AnnData): Input AnnData containing a valid CxG matrix as 'X', and annotations in the 'obs' dataframe.
        celltype_labels_column (str, optional): Column containing annotations in the 'obs' dataframe. Defaults to 'azimuth_preds'.
        aggregation_method (str, optional): `mean/median/count`. Defaults to 'mean'.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.
    
    Returns:
        pd.DataFrame: Aggregated CxG dataframe containing each Cell-Type label mapped to `mean/median/count` values of all gene-expressions.
    """
    adata.X = adata.to_df()

    celltype_labels = adata.obs[celltype_labels_column]
    if verbose:  print(f'celltype_labels are: {celltype_labels}')

    cxg_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    if verbose:  print(f'new_cxg_df has columns: {cxg_df.columns}')

    cxg_annotated_df = pd.concat([cxg_df, celltype_labels], axis=1)
    if verbose:  print(f'Aggregation based on: {aggregation_method}')

    agg_cxg_df = cxg_annotated_df.groupby(by=[celltype_labels_column], as_index=False).agg(aggregation_method)
    if verbose:  print(agg_cxg_df.head())

    return agg_cxg_df



def get_top_and_bottom_K_genes(agg_cxg_df, celltype_labels_column='azimuth_preds', top_k=10, bottom_k=0, verbose=False):
    """Filters the aggregated CxG dataframe to retrieve the top-K Cell-Type labels showing highest `mean/median/count` gene-expressions.

    The `agg_cxg_df` dataframe should contain columns: [`celltype_labels_column`, `gene1_val`, `gene2_val`, ...] 
    ____ where `gene1_val`, `gene2_val`, ... indicate aggregated expression values for that celltype-label.

    Args:
        agg_cxg_df (pd.DataFrame): Aggregated CxG dataframe containing each Cell-Type label mapped to `mean/median/count` values of all gene-expressions.
        celltype_labels_column (str, optional): Column containing annotations in the 'obs' dataframe. Defaults to 'azimuth_preds'.
        top_k (int, optional): Top k genes to retrieve. Defaults to 10.
        bottom_k (int, optional): Bottom k genes to retrieve. Defaults to 0.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        dict: A dictionary containing `{'CT-label' : { 'index':[top-K genes] ,  'values':[top-K gene-expression values] }  }`
    """
    agg_cxg_T_df = agg_cxg_df.set_index(celltype_labels_column).T
    ct_labels = agg_cxg_T_df.columns.values.tolist()
    most_expressed_genes_dict = {}
    for ct_label in ct_labels:
        most_expressed_genes = agg_cxg_T_df.loc[agg_cxg_T_df[ct_label]>0, ct_label].sort_values(ascending=False)

        if verbose:  print(f'Fetching top-{top_k} genes expressed by `{ct_label}`...')
        top_k_expressed_genes = most_expressed_genes.head(n=top_k).reset_index()
        if verbose:  print(f'\ntop_most_expressed_genes={top_k_expressed_genes}')

        if verbose:  print(f'Fetching bottom-{bottom_k} genes expressed by `{ct_label}`...')
        bottom_k_expressed_genes = most_expressed_genes.tail(n=bottom_k).reset_index()
        if verbose:  print(f'\n\nbottom_most_expressed_genes={bottom_k_expressed_genes}')
        
        most_expressed_genes_dict[ct_label] = pd.concat([top_k_expressed_genes, bottom_k_expressed_genes], axis=0).to_dict(orient='list')
    return most_expressed_genes_dict








def create_agreeability_matrix(most_expressed_genes, canonical_markers, derived_markers, agreeability_map, verbose=True):
    """Creates a matrix of the MxK (Cell by Gene) subset with a character-value in each cell to indicate overlayed results of multiple sources.
    The character in each cell will be defined by the `agreeability_map` variable. For ex- 
    ```python
        {
            'both_matching' : '✖',
            'derived_marker' : 'Đ',
            'canonical_marker' : 'Ȼ'
        }
    ```

    Args:
        most_expressed_genes (dict): A dictionary containing `{'CT-label' : { 'index':[top-K genes] ,  'values':[top-K gene-expression values] }  }`
        canonical_markers (dict): A dictionary containing `{'CT-label' : [top-K genes]}` fetched from a reliable source. ex- Azimuth maintains top-10 canonical markers per CT.
        derived_markers (dict): A dictionary containing `{'CT-label' : [top-K genes]}` fetched algorithmically. ex- NSForest performs feature selection per CT-cluster.
        agreeability_map (dict): A dictionary to indicate what character to use for which scenario.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        list, list: List of gene-names, and List/Matrix of agreeability indicators.
    """
    gene_names = []
    agreeability_indicators = []
    
    for celltype, top_k_genes in most_expressed_genes.items():
        celltype = celltype.lower().replace('φ','ï†')
        gene_names.append(top_k_genes['index'])
        agreeability_indicators_celltype = ['' for _ in top_k_genes['index']]
        for i, gene in enumerate(top_k_genes['index']):
            if celltype not in derived_markers and celltype not in canonical_markers:
                if verbose:  print(f'{celltype} not found in derived_markers nor canonical_markers')
                break
            if gene in derived_markers.get(celltype,[]) and gene in canonical_markers.get(celltype,[]):
                agreeability_indicators_celltype[i] = agreeability_map['both_matching']
            elif gene in derived_markers.get(celltype,[]):
                agreeability_indicators_celltype[i] = agreeability_map['derived_marker']
            elif gene in canonical_markers.get(celltype,[]):
                agreeability_indicators_celltype[i] = agreeability_map['canonical_marker']
        agreeability_indicators.append(agreeability_indicators_celltype)
    return gene_names, agreeability_indicators







def get_matching_count_summaries(predictions1, predictions2, verbose=False):
    """Get stats for one-vs-one prediction comparisons.

    Args:
        predictions1 (pd.Series): Predicted labels from first algorithm (ex- Azimuth).
        predictions2 (pd.Series): Predicted labels from second algorithm (ex- CellTypist).
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        pd.DataFrame: Summary dataframe containing true/false counts for one-vs-one comparisons of both predicted labels.
    """
    predictions1 = predictions1.reset_index(drop=True)
    predictions2 = predictions2.reset_index(drop=True)
    comparisons = (predictions1==predictions2).value_counts().reset_index()
    if verbose:  print(comparisons)
    comparisons.columns = ['match', 'count']
    comparisons.loc[:, 'match'] = comparisons['match'].apply(lambda x : 'Yes' if x else 'No')
    return comparisons