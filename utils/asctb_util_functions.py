import base64, json, sys, re, os, pandas as pd, json, http.client
from pprint import pprint
from requests import get



def get_ccf_reporter_sheet_config(verbose=False):
    """Accesses the Github repo for CCF-reporter, pulls the sheet-config at below URL.
    `https://github.com/hubmapconsortium/ccf-asct-reporter/blob/main/projects/v2/src/assets/sheet-config.json`
    
    This sheet contains links for all ASCT+B organs --> dataset URL.

    Args:
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.
    
    Return:
        list(dict): Python version of the `sheet_config.json` file.
    """
    USER = 'hubmapconsortium'
    REPOSITORY_NAME = 'ccf-asct-reporter'
    FILE_PATH = '/projects/v2/src/assets/sheet-config.json'
    sheet_config = None
    response = get(url=f'https://api.github.com/repos/{USER}/{REPOSITORY_NAME}/contents/{FILE_PATH}')
    if response.status_code == 200:
        json_response = response.json()
        decoded_content = base64.b64decode(json_response['content']) # Github returns response in b64 encoding
        json_string = decoded_content.decode('utf-8')
        sheet_config = json.loads(json_string)
    else:
        print(f'Error {response.status_code}! Something went wrong while trying to read the CCF-Reporter Github repo sheet-config file...')
    if verbose:  print(f'sheet_config = \n{sheet_config}')
    return sheet_config



def get_asctb_data_url(asctb_organ='Lung', asctb_organ_version='v1.1', verbose=False):
    """Reads the sheet-config from the CCF-reporter Github repo, and parses it to fetch the dataset-link for the specific organ and version.
    
    Args:
        asctb_organ (str, optional): Defaults to 'Lung'.
        asctb_organ_version (str, optional): Defaults to 'v1.1'.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.
    
    Returns:
        str: `URL` for the specific ASCT+B organ and version.
        str: `sheetId` for the specific ASCT+B organ and version.
        str: `gid` for the required ASCT+B organ and version.
    """
    GOOGLE_SHEETS_BASE_URL = 'https://docs.google.com/spreadsheets/d/'
    SHEET_CONFIG = get_ccf_reporter_sheet_config(verbose=verbose)
    if not SHEET_CONFIG:
        sys.exit('Couldnt access the CCF-Reporter Github!')

    for asctb_dataset in SHEET_CONFIG:
        if asctb_dataset['name'].lower() == asctb_organ.lower():
            if verbose:  pprint(asctb_dataset)
            for version_metadata in asctb_dataset['version']:
                if version_metadata['viewValue'] == asctb_organ_version:
                    version_name = version_metadata['value']
                    google_sheets_url = GOOGLE_SHEETS_BASE_URL + version_metadata['sheetId']
                    if verbose:  print(version_name, google_sheets_url)
                    return google_sheets_url, version_metadata['sheetId'], version_metadata['gid']
    return None




def get_asctb_data(asctb_organ='Lung', asctb_organ_version='v1.2', verbose=False):
    """Fetches the URL for the specific ASCTB organ version, and reads the entire Google sheet data.
    Uses the CCF-reporter API documentation available here: mmpyikxkcp.us-east-2.awsapprunner.com

    Args:
        asctb_organ (str, optional): Defaults to 'Lung'.
        asctb_organ_version (str, optional): Defaults to 'v1.2'.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        dict: Json response containing rows of the ASCTB organ version dataset.
    """
    _, sheet_id, gid = get_asctb_data_url(asctb_organ=asctb_organ, asctb_organ_version=asctb_organ_version, verbose=verbose)
    conn = http.client.HTTPSConnection("mmpyikxkcp.us-east-2.awsapprunner.com")
    headers = { 'Content-Type': "application/json" }
    conn.request('GET', f'/v2/{sheet_id}/{gid}/', headers=headers)
    res = conn.getresponse()
    data = res.read().decode("utf-8")
    response = json.loads(data)
    return response


def parse_asctb_ct_vs_bg(response, verbose=False):
    """Parse each row of the ASCTB organ version json data and return a dictionary for:
        
        `{ Cell-Type name : Set of canonical biomarkers reported by ASCT+B authors }`

    Args:
        response (dict): Json response containing rows of the ASCTB organ version dataset.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        dict: Dictionary containing `{ Cell-Type name : Set of canonical biomarkers reported by ASCT+B authors }`
    """
    hmap = dict()
    for single_row in response['data']:
        biomarker_genes = set([biomarker['name'] for biomarker in single_row['biomarkers']]) # 'biomarker_genes'
        for row in single_row['cell_types']:
            ct_name = row['name'].lower().replace('φ','ï†')
            if ct_name not in hmap:
                hmap[ct_name] = []
            hmap[ct_name].extend(list(biomarker_genes))
    for ct, bgs in hmap.items():
        hmap[ct] = list(set(bgs))
    if verbose:  pprint(hmap)
    return hmap


def merge_asctb_crosswalk_into_derived_markers(translation_hmap, markers_dict, annotation_colnames=['azimuth_preds'], verbose=False):
    """Uses a translation hashmap of the final-crosswalk, and adds a copy of nsforest derived-markers but with ASCT+B naming convention.
    For ex- 
    ```python
        {
            {'azimuth_preds': {'at2' : [list of markers]}, ...},
            ...
        }
    ```
    copies over to
    ```python
        {
            {'azimuth_preds': {'at2' : [list of markers]}, ...},
            {'azimuth_preds_asctb_equivalent': {'proliferating at2 ' : [list of markers]}, ...},
            ...
        }
    ```

    Args:
        translation_hmap (dict): Contains the final-crosswalk translations: `[unique_cts, asctb_equivalent]`.
        markers_dict (dict): Json-structure of `{ source : { celltype : list(markers) } }`.
        annotation_colnames (list, optional): Original list of keys in the `algorithmically_derived_markers`. Defaults to `['azimuth_preds']`.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        dict: Json-structure of `{ source : { celltype : list(markers) } }`, similar to `algorithmically_derived_markers` including new asctb-equivalent-keys.
    """
    # Let's create a dictionary copy of all 3 predictions-and-markers, but with predictions translated to ASCTB naming conventions
    for colname in annotation_colnames:
        new_annotation_colname = f'{colname}_asctb_equivalent'
        if verbose:  print(f'Creating a new {new_annotation_colname} dictionary of nsforest-markers...')
        asctb_equivalent_markers_dict = {}
        for ctlabel, markers in markers_dict[colname].items():
            ctlabel = ctlabel.lower().replace('φ','ï†')
            ctlabel_asctb_equivalent = translation_hmap[ctlabel].lower().replace('φ','ï†')
            if verbose:  print(ctlabel, ctlabel_asctb_equivalent, markers)
            asctb_equivalent_markers_dict[ctlabel_asctb_equivalent] = markers
        markers_dict[new_annotation_colname] = asctb_equivalent_markers_dict

    if verbose:  print(markers_dict.keys())
    return markers_dict
