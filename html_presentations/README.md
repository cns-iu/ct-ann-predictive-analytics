# HTML presentations

This folder contains the presentable format of the notebooks containing visualizations.



### 3 main driver notebooks:


a. `1_Overview_and_Fetch_predictions.ipynb` - If modularized efficiently, can use this to fetch predictions from containerized microservices.
b. `2_Analyze_Labels_and_Feature_space.ipynb` - Explores the gene-expressions in the query dataset, and overlays the most expressive genes with canonically identified markers (ASCT+B or Azimuth), and algorithmically identified markers (NS-Forest). Uses an expertly curated "Crosswalk" for data harmonization.
c. `3_Analyze_scoring_metrics.ipynb` - Explores the distributions of scores (includes facets of metadata like assay-type), and naming conventions amongst 3 modalities of annotating single-cell data. Uses an expertly curated "Crosswalk" for data harmonization.
d. `5_Sankey_Plot_and_PopV_compartments.ipynb` - Currently a WIP. The Sankey plot is tricky to implement. Still trying out multiple approaches.