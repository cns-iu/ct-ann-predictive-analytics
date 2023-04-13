# ct-ann-predictive-analytics

A companion code-repository of the CT-Annotations paper, for comparing the predictive behaviors of algorithms created using Atlas-like reference datasets.


## Azimuth vs CellTypist vs PopV



### Summary notes of all papers is in “Vikrant-CxG-Notes-k.docx” (Link).

### 3 main driver notebooks:

<ol>
<li>1_Overview_and_Fetch_predictions.ipynb</li>
<li>2_Analyze_Labels_and_Feature_space.ipynb</li>
<li>3_Analyze_scoring_metrics.ipynb</li>
</ol>

I usually export these to html for presentations in the `"/html_presentations/"` folder.
I’m maintaining this as my own private Github Repo for now, but TBD on how and when to open-source.

### Main folders with modularized code for fetching predictions and results:

A separate folder each for fetching outputs from CellTypist/Azimuth/PopV by using the driver notebook #a.

> PopV code is tricky to make modular, kept it aside for now.


A separate folder for fetching most discriminative features from NS-Forest algorithm by using the driver notebook #b.


### Input and Outputs of AnnData objects will be in “Datasets” folder:

1. Naming convention of input file should be: Datasets/{data_name}/{data_name}.h5ad
2. Output of predictions for each of the 3 algorithms in step #2 will be stored as “Datasets/{data_name}/{algorithm_name}_preds.csv”, or “Datasets/{data_name}/{algorithm_name}_preds.tsv”.

#### Datsets explored so far:

LungMap.h5ad- Source data (https://celltype.info/CAPinitialRelease/LungMAP-Human-data-from-a-broad-age-healthy-donor-group/3)
PopV_Reference_Lung.h5ad- Training dataset for PopV model (https://www.dropbox.com/s/2kuzdamjevev2ci/Lung.h5ad?dl=1)


1. LungMap.h5ad - https://celltype.info/CAPinitialRelease/LungMAP-Human-data-from-a-broad-age-healthy-donor-group/3

2. TS_Lung.h5ad filtered and cleaned - https://www.dropbox.com/s/2kuzdamjevev2ci/Lung.h5ad?dl=1 (notebook is at https://colab.research.google.com/drive/1Yw4ZDMoPgXNiC1ZQo2eS75Sw8Y_23rrb?usp=sharing#scrollTo=Zty7C8HAZwwr)

3. TS_Lung.h5ad- https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219