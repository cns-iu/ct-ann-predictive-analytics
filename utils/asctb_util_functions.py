import json, sys, json, http.client

from pprint import pprint
from asctb_ct_label_mapper.utilities.asctb_data_wrangling import get_ccf_reporter_sheet_config, get_asctb_data_url




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
