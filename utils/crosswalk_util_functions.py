import base64, json, sys, re, os, pandas as pd, json, http.client, numpy as np
from pprint import pprint
from requests import get

def fetch_and_parse_crosswalk_table(crosswalk_filename='Azimuth_CellTypist_PopV_Lung_ASCTB_Crosswalks.csv', raw_labels_colname='raw_input_column', asctb_crosswalk_colname='translation_column', tgt_crosswalk_colname='asctb_equivalent', verbose=False):
    """Processes the final-crosswalk data containing SME feedback for translating raw-labels into ASCTB naming convention.

    Args:
        crosswalk_filename (str, optional): Defaults to 'Azimuth_CellTypist_PopV_Lung_ASCTB_Crosswalks.csv'.
        raw_labels_colname (str, optional): Defaults to 'raw_input_column'.
        asctb_crosswalk_colname (str, optional): Defaults to 'translation_column'.
        tgt_crosswalk_colname (str, optional): Defaults to 'asctb_equivalent'.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        pd.DataFrame: Contains the final-crosswalk information.
    """
    # Merge the aggregated data with the translations file created using our ASCTB-Mapper package with finalized SME feedback
    crosswalk_df = pd.read_csv(crosswalk_filename)

    crosswalk_df['source'] = crosswalk_df['source'].replace('Azimuth-HLCAv2', 'azimuth').replace('PopV-Lung', 'popv').replace('CellTypist-Lung', 'celltypist')

    crosswalk_df[raw_labels_colname] = crosswalk_df[raw_labels_colname].str.lower().replace('φ','ï†')

    crosswalk_df[tgt_crosswalk_colname] = crosswalk_df[asctb_crosswalk_colname].replace('?', np.nan)
    crosswalk_df.loc[crosswalk_df[tgt_crosswalk_colname].isna(), tgt_crosswalk_colname] = crosswalk_df.loc[crosswalk_df[tgt_crosswalk_colname].isna(), 'best_matched_asctb_label']
    crosswalk_df[tgt_crosswalk_colname] = crosswalk_df[tgt_crosswalk_colname].str.lower().replace('φ','ï†')
    return crosswalk_df



def get_crosswalk_translation_hmap(crosswalk_df):
    """Create a hashmap of Gloria's crosswalk for cell-type labels from Azimuth/CellTypist/PopV -> ASCTB naming conventions.

    Args:
        crosswalk_df (pd.DataFrame): DataFrame containing final crosswalk information. Essential columns: `[unique_cts, asctb_equivalent]`.
    """
    translation_hmap = dict(
        zip(
            crosswalk_df['unique_cts'], crosswalk_df['asctb_equivalent']
        )
    )
    return translation_hmap

