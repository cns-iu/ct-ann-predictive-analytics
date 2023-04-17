import plotly.graph_objects as go, plotly.express as px, numpy as np
from utils.utiltity_functions import create_summary_obs

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
