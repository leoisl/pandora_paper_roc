import dash
from dash.dependencies import Input, Output, State

from plot_helpers import *



# set up the app
external_stylesheets = ["css/style.css"]
dash_app = dash.Dash("ROC_pandora", assets_folder="assets", external_stylesheets=external_stylesheets)

# set up the layout
dash_app.layout = html.Div([
    html.Div(html.H1(f"Pandora ROC evaluation"), style={'textAlign': "center"}),
    html.H3("1. Choose analysis:"),
    dcc.Dropdown(
        id='plots-dropdown',
        options=get_available_analyses()
    ),

    html.H3("2. Apply filters and press go:"),
    get_empty_div_for_the_filters(filter_name=f"tool", filter_description="Tool"),
    get_empty_div_for_the_filters(filter_name=f"dataset", filter_description="Dataset coverage filter (pandora only)"),
    get_empty_div_for_the_filters(filter_name=f"coverage", filter_description="Coverage filter (pandora only)"),
    get_empty_div_for_the_filters(filter_name=f"strand_bias", filter_description="Strand bias filter (pandora only)"),
    get_empty_div_for_the_filters(filter_name=f"gaps", filter_description="Gaps filter (pandora only)"),
    html.Button('Go', id='button', className="row", style={"display": "block", "width": "60%", "margin-left": "auto",
                                     "margin-right": "auto"}),

    html.H3("3. See data (might take some seconds to plot):"),
    html.H3("Proportion ROC curve:"),
    html.Div(children=[
        dcc.Loading(id=f"loading_graph_proportion",
                    children=[html.Div(id="graph_proportion")],
                    type="graph")]),
    html.H3("Raw numbers ROC curve:"),
    html.Div(children=[
        dcc.Loading(id=f"loading_graph_raw",
                    children=[html.Div(id="graph_raw")],
                    type="graph")]),

    html.H3("Extra informations:"),
    html.H3("Data with no GT conf filtering:"),
    html.Div(children=[
        dcc.Loading(id=f"loading_data_with_no_gt_conf_filter",
                    children=[html.Div(id="data_with_no_gt_conf_filter")],
                    type="graph")])
    ])



# set up the callbacks
@dash_app.callback(
    [Output('tool_checklist', 'options'),
     Output('tool_checklist', 'value'),
     Output('dataset_checklist', 'options'),
     Output('dataset_checklist', 'value'),
     Output('coverage_checklist', 'options'),
     Output('coverage_checklist', 'value'),
     Output('strand_bias_checklist', 'options'),
     Output('strand_bias_checklist', 'value'),
     Output('gaps_checklist', 'options'),
     Output('gaps_checklist', 'value'),
     ],
    [Input('plots-dropdown', 'value')])
def update_all_checklists(plots_value):
    df = get_df_given_analysis_name(plots_value)
    list_of_output_tuples = update_tool_checklist(df), update_dataset_checklist(df), update_coverage_checklist(df),\
           update_strand_bias_checklist(df), update_gaps_checklist(df)
    flattened_output = list(itertools.chain(*list_of_output_tuples))
    return flattened_output



@dash_app.callback(
    [Output(f'graph_proportion', 'children'),
     Output(f'graph_raw', 'children'),
     Output('data_with_no_gt_conf_filter', 'children')],
    [Input('plots-dropdown', 'value'),
     Input('button', 'n_clicks')],
    [
     State(component_id='tool_checklist', component_property='value'),
     State(component_id='dataset_checklist', component_property='value'),
     State(component_id='coverage_checklist', component_property='value'),
     State(component_id='strand_bias_checklist', component_property='value'),
     State(component_id='gaps_checklist', component_property='value')
     ])
def update_all_graphs (*args):
    df = get_df_and_check_args_for_graph (*args)
    return get_graph_proportion(df, *args), get_graph_raw(df, *args), get_data_table_with_no_gt_conf_filter(df)


if __name__ == '__main__':
    dash_app.run_server(debug=True)