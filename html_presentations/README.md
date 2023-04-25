# HTML presentations

This folder contains the presentable format of the notebooks containing visualizations.



### 3 main driver notebooks:

1. `1_Overview_and_Fetch_predictions.ipynb` - If modularized efficiently, can use this to fetch predictions from containerized microservices.
2. `2_Analyze_Labels_and_Feature_space.ipynb` - Explores the gene-expressions in the query dataset, and overlays the most expressive genes with canonically identified markers (ASCT+B or Azimuth), and algorithmically identified markers (NS-Forest). Uses an expertly curated "Crosswalk" for data harmonization.
3. `3_Analyze_scoring_metrics.ipynb` - Explores the distributions of scores (includes facets of metadata like assay-type), and naming conventions amongst 3 modalities of annotating single-cell data. Uses an expertly curated "Crosswalk" for data harmonization.


### Ad hoc visualizations and code:

4. `5_Sankey_Plot_and_PopV_compartments.ipynb` - WIP. The Sankey plot is tricky to implement. Still trying out multiple approaches.
5. `6_Sankey_plot.R` - Implemented in R for: Celltypist predicted annotations <-> Azimuth predicted annotations <-> PopV predicted annotations.