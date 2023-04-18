import pandas as pd
import ptitprince as pt, matplotlib.pyplot as plt
import plotly.graph_objects as go, plotly.express as px, numpy as np

from utils.utiltity_functions import create_agreeability_matrix, create_summary_obs, get_matching_count_summaries




def plot_obs_barchart(adata, max_unique=50, dataset_name='Dummy', verbose=False):
    """Create a summary dataframe and plot a barchart to show general distribution of metadata available in the AnnData object.

    Args:
        adata (AnnData): Input AnnData object containing an `obs` dataframe that we can summarize and view.
        max_unique (int, optional): Threshold upto which we show all values on the interactive plot. Defaults to 50.
        dataset_name (str, optional): Textual name for the input AnnData dataset. Defaults to 'Dummy'.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        plotly.Figure: A figure that we can plot later using `fig.show()`
    """
    summary_obs_df = create_summary_obs(adata, max_unique=max_unique, verbose=verbose)
    fig = go.Figure([
        go.Bar(
            x=summary_obs_df['colname'], 
            y=summary_obs_df['nunique'], 
            name='Counts for all metadata within AnnData.Obs',
            marker_color='rgb(39, 43, 102)', 
            marker_line_color='rgb(8,48,107)',
            hovertext=summary_obs_df['unique'].fillna('## Too many values ##'),
            hovertemplate=\
                    'Column: %{x}' + \
                        '<br>Unique-%{y}' + \
                            '<br>%{hovertext}',
            marker_line_width=.8,
            opacity=0.6
        )
    ])

    fig.update_layout(
        yaxis={'title':'Unique Counts'},
        xaxis={'title':'CellType label', 'categoryorder':'category ascending'}, 
        title=f'Summarized metadata from the "{dataset_name}" dataset',
        width=1000,
        height=1000,
        updatemenus=[
                dict(
                    buttons=list([
                        dict(
                            args=[{'yaxis.type': 'linear'}],
                            label='Linear scale',
                            method='relayout'
                        ),
                        dict(
                            args=[{'yaxis.type': 'log'}],
                            label='Log Scale',
                            method='relayout'
                        )
                    ]),
                    direction='down',
                    showactive=True,
                    x=1.,
                    y=1.1
                )
            ]
        )

    fig.update_yaxes(tick0=0, dtick=1)
    if verbose:  fig.show()
    return fig





def plot_obs_treemap(adata, max_unique=50, dataset_name='Dummy', verbose=False):
    """Create a summary dataframe and plot a treemap to show general distribution of metadata available in the AnnData object.

    Args:
        adata (AnnData): Input AnnData object containing an `obs` dataframe that we can summarize and view.
        max_unique (int, optional): Threshold upto which we show all values on the interactive plot. Defaults to 50.
        dataset_name (str, optional): Textual name for the input AnnData dataset. Defaults to 'Dummy'.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        plotly.Figure: A figure that we can plot later using `fig.show()`
    """
    summary_obs_df = create_summary_obs(adata, max_unique=max_unique, verbose=verbose)
    summary_obs_df.loc[:, 'nunique_log'] = np.log(summary_obs_df.loc[:, 'nunique'])
    
    
    fig = px.treemap(
        summary_obs_df, 
        path=['colname'], 
        values='nunique',
        color='colname',
        hover_data={'nunique': True, 'unique' : True},
        hover_name='colname'
    )
    fig.update_layout(
        margin = dict(t=50, l=25, r=25, b=25),
        yaxis={'title':'Unique Counts'},
        xaxis={'title':'CellType label', 'categoryorder':'category ascending'}, 
        title=f'Summarized metadata from the "{dataset_name}" dataset',
        width=800,
        height=800
    )

    fig.update_traces(marker={'colorscale':'viridis'})
    if verbose:  fig.show()
    return fig







def create_raincloudplot_top_k_bottom_k_genes(
        adata, most_expressed_genes, 
        celltype_labels_column='azimuth_preds', aggregation_method='mean', filter_sparse_cell_types=False, 
        top_k=100, bottom_k=10,
        verbose=False
    ):
    """Filters out the sparse Cell-Type clusters with <5 cells. These would not give much info in the RainCloud plot.
    
    For each unique and valid Cell-Type label cluster, make a Raincloud plot for the top-K gene-expressions in a row.

    Args:
        adata (Scanpy.AnnData): Input AnnData containing a valid CxG matrix as 'X', and annotations in the 'obs' dataframe.
        celltype_labels_column (str, optional): Column containing annotations in the 'obs' dataframe. Defaults to 'azimuth_preds'.
        aggregation_method (str, optional): `mean/median/mode/count`. Defaults to 'mean'.
        filter_sparse_cell_types (bool, optional): Flag to indicate if . Defaults to False.
        top_k (int, optional): Top-K genes to visualize. Defaults to 100.
        bottom_k (int, optional): Bottom-K genes to visualize. Defaults to 10.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.
    """
    
    if filter_sparse_cell_types:
        non_sparse_cell_types = adata.obs[celltype_labels_column].value_counts()>=5
        fig, ax = plt.subplots(figsize=(50,80), nrows=non_sparse_cell_types[non_sparse_cell_types==True].count())
        ordered_cell_types = reversed(sorted(non_sparse_cell_types[non_sparse_cell_types==True].index))
    else:
        fig, ax = plt.subplots(figsize=(50,80), nrows=len(adata.obs[celltype_labels_column].unique()))
        ordered_cell_types = reversed(sorted(adata.obs[celltype_labels_column].unique()))
        
    if verbose:  print(ordered_cell_types)

    for i, ct_label in enumerate(ordered_cell_types):
        if verbose:  print('Fetching the CxG subset-dataframe for current CT-annotation and top K gene-expressions: CTxK')
        ct_most_exp_gene_distributions = adata[adata.obs[celltype_labels_column]==ct_label, most_expressed_genes[ct_label]['index']].X
        ct_most_exp_gene_distributions = pd.DataFrame(ct_most_exp_gene_distributions, columns=most_expressed_genes[ct_label]['index'])
        if verbose:  print(f'ct_label[{i}]={ct_label} has shape: {ct_most_exp_gene_distributions.shape}...')

        if filter_sparse_cell_types and ct_most_exp_gene_distributions.shape[0]<5:
            if verbose:  print(f'ct_label={ct_label} observed in <5 cells. Skipping its Raincloud plot...')
            fig.delaxes(ax[i])
            continue

        order = pd.DataFrame(ct_most_exp_gene_distributions.agg(aggregation_method), columns=[f'{aggregation_method}_expression'])
        order = order.sort_values(by=[f'{aggregation_method}_expression'], ascending=[False]).index
        
        ax[i] = pt.RainCloud(
            ax=ax[i], 
            data=ct_most_exp_gene_distributions,
            palette='crest' if i%2==0 else 'flare',
            point_size=1,       # cloud point-size
            bw=0.05,            # bandwidth-param for Kernel-Density-Estimation plot
            cut=0,              # cut-param for Kernel-Density-Estimation plot
            orient='v',
            width_viol=1,
            pointplot=True,     # Line connecting the means across the RainCloud plots
            linecolor='grey' if i%2==0 else 'orange',
            order=list(order)
        )
        ax[i].set_ylabel(ct_label.replace(' ','\n'))
        labels = ax[i].set_xticklabels(labels=list(order), rotation=0, fontsize=8)
        for i, label in enumerate(labels):
            label.set_y(label.get_position()[1]*0.0001)
    
    fig_title = f'Top-{top_k}{f" and Bottom-{bottom_k}" if bottom_k else ""} {aggregation_method} gene-expressions per CT-annotation from {celltype_labels_column.replace("_preds","").capitalize()}'
    fig.suptitle(
        t=fig_title,
        x=.5, 
        y=.9, 
        verticalalignment='top', 
        size=50
    )
    return fig









def create_heatmap_top_k_bottom_k_genes(
        most_expressed_genes, canonical_markers_one, canonical_markers_two, derived_markers, 
        top_k=10, bottom_k=0, 
        aggregation_method='median', celltype_labels_colname='azimuth_preds', query_dataset_name='LCA',
        verbose=False
    ):
    """Creates a dataframe containing CT-Label, Aggregated-Gene-expression-value (CxG).
    Plots a heatmap for this aggregated view generated from the input dictionary.

    Args:
        most_expressed_genes (dict): A dictionary containing `{'CT-label' : { 'index':[top-K genes] ,  'values':[top-K gene-expression values] }  }`
        canonical_markers (dict): A dictionary containing `{'CT-label' : [top-K genes]}` fetched from a reliable source. ex- Azimuth maintains top-10 canonical markers per CT.
        derived_markers (dict): A dictionary containing `{'CT-label' : [top-K genes]}` fetched algorithmically. ex- NSForest performs feature selection per CT-cluster.
        top_k (int, optional): Top k genes to retrieve. Defaults to 10.
        bottom_k (int, optional): Bottom k genes to retrieve. Defaults to 0.
        aggregation_method (str, optional): Aggregation method used previously in the `get_aggregated_cellbygene_matrix()`. Defaults to 'median'.
        celltype_labels_colname (str, optional): Annotation column in the AnnData object for which we've identified the `most_expressed_genes`. Defaults to 'azimuth_preds'.
        query_dataset_name (str, optional): Name of the input dataset. Defaults to 'LCA'.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        plt.Figure
    """
    k_colnames = [f'top marker #{i}' for i in range(1,top_k+1)] + [f'bottom marker #{i}' for i in range(1,bottom_k+1)]
    df = pd.DataFrame(columns=['celltype_labels_column'] + k_colnames)
    df['celltype_labels_column'] = most_expressed_genes.keys()

    for ct_label, gene_exp_dict in most_expressed_genes.items():
        df.loc[df['celltype_labels_column']==ct_label, k_colnames] = gene_exp_dict[ct_label][:top_k+bottom_k]
    df = df.set_index('celltype_labels_column').astype(np.float64)

    # Populate the gene-names for top-k gene expressions (used on hover)
    # Populate the agreeability-indicator marks to indicate if current gene-name is also in canonical marker-list and also in algorithmically derived markers
    canonical_markers_one = {ct.lower().replace('φ','ï†') : marker_genes for ct, marker_genes in canonical_markers_one.items()}
    canonical_markers_two = {ct.lower().replace('φ','ï†') : marker_genes for ct, marker_genes in canonical_markers_two.items()}
    derived_markers = {ct.lower().replace('φ','ï†') : marker_genes for ct, marker_genes in derived_markers.items()}
    
    agreeability_map = {
        'both_matching' : '✖',
        'derived_marker' : 'Đ',
        'canonical_marker' : 'Ȼ'
    }
    if verbose:  print(f'Creating agreability matrix for first set of canonical markers')
    gene_names_one, agreeability_indicators_one = create_agreeability_matrix(most_expressed_genes, canonical_markers_one, derived_markers, agreeability_map, verbose)
    if verbose:  print(f'Creating agreability matrix for second set of canonical markers')
    gene_names_two, agreeability_indicators_two = create_agreeability_matrix(most_expressed_genes, canonical_markers_two, derived_markers, agreeability_map, verbose)
    
    
    fig = go.Figure()
    
    fig.add_trace(
        go.Heatmap(
            y=df.index.values,
            x=k_colnames,
            z=df.to_numpy(),
            name='Azimuth Canonical markers',
            text=agreeability_indicators_one,
            texttemplate="%{text}",
            textfont={"size":12, 'color':'white'},
            hovertext=gene_names_one,
            hovertemplate=\
                'CellType: %{y}' + \
                    '<br>Gene-%{x}: %{hovertext}' + \
                        f'<br>{aggregation_method.capitalize()}-val: '+'%{z:.2f}',
            colorscale=[[0, 'lightgreen'], [1, 'darkblue']]
        )
    )

    fig.add_trace(
        go.Heatmap(
            y=df.index.values,
            x=k_colnames,
            z=df.to_numpy(),
            name='ASCT+B Canonical markers',
            text=agreeability_indicators_two,
            texttemplate="%{text}",
            textfont={"size":12, 'color':'white'},
            hovertext=gene_names_two,
            hovertemplate=\
                'CellType: %{y}' + \
                    '<br>Gene-%{x}: %{hovertext}' + \
                        f'<br>{aggregation_method.capitalize()}-val: '+'%{z:.2f}',
            colorscale=[[0, 'lightgreen'], [1, 'darkblue']]
        )
    )

    base_title = f'Top-{top_k} and Bottom-{bottom_k} genes' if bottom_k else f'Top-{top_k} genes'
    fig.update_layout(
        width=1800,
        height=1800,
        xaxis={'title':f'{base_title} per CellType Cluster'},
        yaxis={'title':'CellType Label'},
        margin={'t': 200},
        title_text=f'{base_title}-{aggregation_method} values per CellType-cluster ({celltype_labels_colname.split("_")[0].capitalize()}) in {query_dataset_name}\
            <br>{agreeability_map["both_matching"]}: Canonical and algorithmically-derived marker\
                <br>{agreeability_map["canonical_marker"]}: Canonical marker\
                    <br>{agreeability_map["derived_marker"]}: Derived marker using an algorithm (ex- NSForest)'
    )

    fig.update_layout(
        updatemenus=[
            {
                'active': 0,
                'direction': 'down',
                'font': {'size' : 15},
                'pad': {'r': 10, 't': 10},
                'showactive': True,
                'xanchor': 'left',
                'x': .825,
                'y': 1.05,
                'yanchor': 'top',
                'buttons': [
                    {
                        'label' : 'ASCT+B Canonical markers',
                        'method' : 'restyle',
                        'args' : [{"visible": [False, True]}]
                    },
                    {
                        'label' : 'Azimuth Canonical markers',
                        'method' : 'restyle',
                        'args' : [{"visible": [True, False]}]
                    }
                ],
            }
        ]
    )

    fig.show()
    return fig








def add_custom_subplot_of_summaries(predictions1, predictions2, fig, i, j, verbose=False):
    """Adds a pie-chart subplot to the existing plotly-express figure.

    Args:
        predictions1 (pd.Series): Predicted labels from first algorithm (ex- Azimuth).
        predictions2 (pd.Series): Predicted labels from second algorithm (ex- CellTypist).
        fig (plotly.graph_objects.Figure): Figure containing multiple subplots.
        i (int): Row-Index for subplot position.
        j (int): Col-Index for subplot position.
        verbose (bool, optional): Flag to indicate logging in verbose mode. Defaults to False.

    Returns:
        pd.DataFrame: Summary dataframe containing true/false counts for one-vs-one comparisons of both predicted labels.
    """
    comparisons = get_matching_count_summaries(predictions1, predictions2, verbose)
    comparisons = comparisons.sort_values(by=['match'])
    fig.add_trace(
        go.Pie(
            labels=comparisons['match'], 
            values=comparisons['count'],
            marker=dict(colors=['orangered','green']),
            opacity=.7
        ),
        i, j
    )
    return comparisons