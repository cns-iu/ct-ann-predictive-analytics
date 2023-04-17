import os, scanpy as sc, pandas as pd, numpy as np
from subprocess import Popen, PIPE



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