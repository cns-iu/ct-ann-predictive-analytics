# ct-ann-predictive-analytics

A companion code-repository of the CT-Annotations paper, for comparing the predictive behaviors of algorithms created using Atlas-like reference datasets.

TBD on how and when to open-source.

----------------------------------------------

## Azimuth vs CellTypist vs PopV for Lung datasets


----------------------------------------------


### Summary of all papers:

All condensed notes are in `“/papers/Vikrant-CxG-Notes-k.docx”`.

----------------------------------------------

### 3 main driver notebooks:

1. `1_Overview_and_Fetch_predictions.ipynb` - If modularized efficiently, can use this to fetch predictions from containerized microservices.
2. `2_Analyze_Labels_and_Feature_space.ipynb` - Explores the gene-expressions in the query dataset, and overlays the most expressive genes with canonically identified markers (ASCT+B or Azimuth), and algorithmically identified markers (NS-Forest). Uses an expertly curated "Crosswalk" for data harmonization.
3. `3_Analyze_scoring_metrics.ipynb` - Explores the distributions of scores (includes facets of metadata like assay-type), and naming conventions amongst 3 modalities of annotating single-cell data. Uses an expertly curated "Crosswalk" for data harmonization.



### Ad hoc visualizations and code:

4. `5_Sankey_Plot_and_PopV_compartments.ipynb` - WIP. The Sankey plot is tricky to implement. Still trying out multiple approaches.
5. `6_Sankey_plot.R` - Implemented in R for: Celltypist predicted annotations <-> Azimuth predicted annotations <-> PopV predicted annotations.


Exports of these notebooks and visualizations as html files (for presentability) are in the `"/html_presentations/"` folder.


----------------------------------------------

### Main folders with modularized code for fetching predictions and results:

A separate folder each for fetching outputs from CellTypist/Azimuth/PopV by using the driver notebook #a.

> PopV code is tricky to make modular, kept it aside for now. Use their tutorial notebook `4_Vikrant_PopV_tabula_sapiens_tutorial.ipynb` to fetch predictions.

A separate folder `/nsforest/` for fetching most discriminative features from NS-Forest algorithm by using the driver notebook #b.



The `/utilities/` folder contains reusable code modularized into functions:

1. `asctb_util_functions.py` - Data fetching and wrangling functionality using the [ASCT+B API](https://mmpyikxkcp.us-east-2.awsapprunner.com/#/).
2. `azimuth_util_functions.py` - Data fetching and wrangling functionality for the Azimuth-website's backend [repository](https://github.com/satijalab/azimuth_website/tree/master/static/csv).
3. `crosswalk_util_functions.py` - Data wrangling functionality to parse the local "crosswalk" file present in the "/asctb_mapper/" folder.
4. `utility_functions.py` - General reusable functions for summarizing the Cell-by-Gene query data, fetching top-k and bottom-k genes in the query data, analyzing pairwise predicted-annotations, etc.
5. `plotting.py` - Custom visualizations/deliverables that could be made reusable. Numerous plots in the 3 driver notebookes were not easy to functionalize.




----------------------------------------------

### Input and Outputs of AnnData objects will be in “Datasets” folder:

1. Naming convention of input file should be: `"Datasets/{data_name}/{data_name}.h5ad"`.
2. Output of predictions for each of the 3 algorithms in will be stored as `“Datasets/{data_name}/{algorithm_name}_preds.csv”`, or `“Datasets/{data_name}/{algorithm_name}_preds.tsv”`.
3. Output of NSForest computed for the 3 modalities will be stored as `“Datasets/{data_name}/NSForest_{algorithm_name}_preds.csv”`.
4. The locally stored outputs of #2 and #3 can be used as cache, while using the above notebooks.



----------------------------------------------

#### Query datsets explored so far:

<ol>
<li>LCA.h5ad - <a href="https://www.dropbox.com/s/mrf8y7emfupo4he/LCA.h5ad">Dropbox</a> (more details on Cell-by-Gene portal <a  href="https://cellxgene.cziscience.com/collections/5d445965-6f1a-4b68-ba3a-b8f765155d3a">link</a>)</li>

<li>TS_Lung.h5ad filtered and cleaned - <a href="https://www.dropbox.com/s/2kuzdamjevev2ci/Lung.h5ad?dl=1">Dropbox</a> (PopV Tutorial Notebook is on <a href="https://colab.research.google.com/drive/1Yw4ZDMoPgXNiC1ZQo2eS75Sw8Y_23rrb?usp=sharing#scrollTo=Zty7C8HAZwwr">Colab</a>)</li>

<li>TS_Lung.h5ad- <a href="https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219">FigShare</a></li>

<li>LungMap.h5ad - <a href="https://celltype.info/CAPinitialRelease/LungMAP-Human-data-from-a-broad-age-healthy-donor-group/3">Cell Annotation Platform</a></li>

<li>TS_Heart.h5ad - to be used as a Heart reference for creating CellTypist-Heart model - <a href="https://www.nature.com/articles/s41586-020-2797-4#data-availability"></a></li>

<li>Tucker_Heart_4_chambers.h5ad - <a href="https://nam12.safelinks.protection.outlook.com/?url=https%3A%2F%2Fpubmed.ncbi.nlm.nih.gov%2F32403949%2F&data=05%7C01%7Cvikdeshp%40iu.edu%7C5d9f524ad2ba4ef59ba208db2c04014f%7C1113be34aed14d00ab4bcdd02510be91%7C0%7C0%7C638152167935828448%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&sdata=OntW49vWEWw3iEpURAjEG5P63zE2zErgoqBz7mEuKgQ%3D&reserved=0">Paper</a> and <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7666104/#S6title">Dataset.</a></li>
</ol>

----------------------------------------------
