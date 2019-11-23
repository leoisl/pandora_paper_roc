import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output


def get_div_for_the_filters(filter_name, filter_values):
    return html.Div([html.H4(f"{filter_name} filters:"),
              dcc.Checklist(
                  id=f"{filter_name}_checklist",
                  options=[{'label': f'{filter_name} {threshold}', 'value': threshold} for
                           threshold in filter_values],
                  value=filter_values,
                  labelStyle={'display': 'inline-block', 'margin-right': '30px'}
              )],
             className="row", style={"display": "block", "width": "60%", "margin-left": "auto",
                                     "margin-right": "auto"})


def create_interactive_visualisation_app(ROC_data_path):
    # load df
    df = pd.read_csv(ROC_data_path, sep="\t")

    # load filters
    tools = df["tool"].unique()
    dataset_coverages = df["coverage"].unique()
    coverage_filters = df["coverage_threshold"].unique()
    strand_bias_filters = df["strand_bias_threshold"].unique()
    gaps_filters = df["gaps_threshold"].unique()

    # set up the app
    external_stylesheets = ["css/style.css"]
    dash_app = dash.Dash("ROC_pandora", assets_folder="assets", external_stylesheets=external_stylesheets)

    # set up the layout
    dash_app.layout = html.Div([html.Div([html.H1("Pandora ROC evaluation")], style={'textAlign': "center"}),
                          get_div_for_the_filters("Tool", tools),
                          get_div_for_the_filters("Dataset coverage", dataset_coverages),
                          get_div_for_the_filters("Coverage", coverage_filters),
                          get_div_for_the_filters("Strand bias", strand_bias_filters),
                          get_div_for_the_filters("Gaps", gaps_filters),
                          html.Div([dcc.Graph(id="my-graph", style={'height': '1000px', 'width': '1000px'})]),
                          ], className="container")

    # register the update callback
    @dash_app.callback(
        Output('my-graph', 'figure'),
        [Input(component_id='Tool_checklist', component_property='value'),
         Input(component_id='Dataset coverage_checklist', component_property='value'),
         Input(component_id='Coverage_checklist', component_property='value'),
         Input(component_id='Strand bias_checklist', component_property='value'),
         Input(component_id='Gaps_checklist', component_property='value')])
    def update_figure(tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values, strand_bias_checklist_values, gaps_checklist_values):
        trace = []
        for tool in tool_checklist_values:
            for dataset_coverage in dataset_coverage_checklist_values:
                for coverage_threshold in coverage_checklist_values:
                    for strand_bias_threshold in strand_bias_checklist_values:
                        for gaps_threshold in gaps_checklist_values:
                            df_for_label = df.query(f"tool==\"{tool}\" & coverage==\"{dataset_coverage}\" & coverage_threshold=={coverage_threshold} & strand_bias_threshold=={strand_bias_threshold} & gaps_threshold=={gaps_threshold}")
                            trace_name = f"Tool: {tool}, Coverage: {dataset_coverage}, Coverage threshold: {coverage_threshold}, Strand bias threshold: {strand_bias_threshold}, Gaps threshold: {gaps_threshold}"
                            trace.append(go.Scatter(x=df_for_label["error_rate"], y=df_for_label["recall"], name=trace_name, mode='lines',
                                                    marker={'size': 8, "opacity": 0.6, "line": {'width': 0.5}}, ))

        return {"data": trace,
                "layout": go.Layout(title="Pandora ROC (selected data)",
                                    yaxis={"title": "Recall"}, xaxis={"title": "Error rate"},
                                    legend_orientation="h")}

    return dash_app


if __name__ == '__main__':
    dash_app = create_interactive_visualisation_app(ROC_data_path ="ROC_data_gt_min_0.0_step_1.0_max_1000000000.0.tsv")
    dash_app.run_server(debug=True)