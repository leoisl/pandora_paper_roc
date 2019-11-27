import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output


def get_div_for_the_filters(filter_name, filter_description, filter_values, selected_values):
    return html.Div([html.H4(f"{filter_description}:"),
              dcc.Checklist(
                  id=f"{filter_name}_checklist",
                  options=[{'label': f'{threshold}', 'value': threshold} for
                           threshold in filter_values],
                  value=selected_values,
                  labelStyle={'display': 'inline-block', 'margin-right': '30px'}
              )],
             className="row", style={"display": "block", "width": "60%", "margin-left": "auto",
                                     "margin-right": "auto"})

def get_highest_error_rate_after_given_recall(df, recall):
    return df.query(f"recall >= {recall}")["error_rate"].max()

def compute_all_possible_traces(df, tools, dataset_coverages, coverage_filters, strand_bias_filters, gaps_filters):
    config_to_trace = {}
    for tool in tools:
        if tool == "snippy_all":
            for snippy_ref in [snippy_ref_tool for snippy_ref_tool in df["tool"].unique() if snippy_ref_tool.startswith("snippy")]:
                df_for_label = df.query(f"tool==\"{snippy_ref}\"")

                trace_name = f"Tool: {snippy_ref}"
                highest_error_rate_after_20_percent_recall = get_highest_error_rate_after_given_recall(df_for_label, 0.2)

                trace = go.Scatter(x=df_for_label["error_rate"], y=df_for_label["recall"], name=trace_name,
                                   mode='lines',
                                   marker={'size': 8, "opacity": 0.6, "line": {'width': 0.5}})

                config_to_trace[snippy_ref] = {
                    "trace": trace,
                    "highest_error_rate": highest_error_rate_after_20_percent_recall
                }
        elif tool.startswith("pandora"):
            for dataset_coverage in dataset_coverages:
                for coverage_threshold in coverage_filters:
                    for strand_bias_threshold in strand_bias_filters:
                        for gaps_threshold in gaps_filters:
                            df_for_label = df.query(
                                f"tool==\"{tool}\" & coverage==\"{dataset_coverage}\" & coverage_threshold==\"{coverage_threshold}\" & strand_bias_threshold==\"{strand_bias_threshold}\" & gaps_threshold==\"{gaps_threshold}\"")

                            trace_name = f"Tool: {tool}, Coverage: {dataset_coverage}, Coverage threshold: {coverage_threshold}, Strand bias threshold: {strand_bias_threshold}, Gaps threshold: {gaps_threshold}"
                            highest_error_rate_after_20_percent_recall = get_highest_error_rate_after_given_recall(df_for_label, 0.2)

                            trace = go.Scatter(x=df_for_label["error_rate"], y=df_for_label["recall"], name=trace_name,
                                                    mode='lines',
                                                    marker={'size': 8, "opacity": 0.6, "line": {'width': 0.5}})


                            config_to_trace[(tool, dataset_coverage, coverage_threshold, strand_bias_threshold, gaps_threshold)] = {
                                "trace": trace,
                                "highest_error_rate": highest_error_rate_after_20_percent_recall
                            }

    return config_to_trace

def remove_value_from_ndarray(array, value):
    return array[array != value]

def keep_pandora_only(tools):
    return [tool for tool in tools if "pandora" in tool]

def create_interactive_visualisation_app(ROC_data_path):
    # load df
    df = pd.read_csv(ROC_data_path, sep="\t")

    # load filters
    tools = df["tool"].unique()
    dataset_coverages = df["coverage"].unique()
    coverage_filters = df["coverage_threshold"].unique()
    strand_bias_filters = df["strand_bias_threshold"].unique()
    gaps_filters = df["gaps_threshold"].unique()

    tools = keep_pandora_only(tools)
    tools.append("snippy_all")
    dataset_coverages = remove_value_from_ndarray(dataset_coverages, "all")
    coverage_filters = remove_value_from_ndarray(coverage_filters, "Not_App")
    strand_bias_filters = remove_value_from_ndarray(strand_bias_filters, "Not_App")
    gaps_filters = remove_value_from_ndarray(gaps_filters, "Not_App")


    # set up the app
    external_stylesheets = ["css/style.css"]
    dash_app = dash.Dash("ROC_pandora", assets_folder="assets", external_stylesheets=external_stylesheets)

    # set up the layout
    dash_app.layout = html.Div([html.Div([html.H1("Pandora ROC evaluation")], style={'textAlign': "center"}),
                          get_div_for_the_filters(filter_name="tool", filter_description="Tool", filter_values=tools, selected_values=tools),
                          get_div_for_the_filters(filter_name="dataset", filter_description="Dataset coverage filter (pandora only)", filter_values=dataset_coverages, selected_values=dataset_coverages),
                          get_div_for_the_filters(filter_name="coverage", filter_description="Coverage filter (pandora only)", filter_values=coverage_filters, selected_values=[coverage_filters[0]]),
                          get_div_for_the_filters(filter_name="strand_bias", filter_description="Strand bias filter (pandora only)", filter_values=strand_bias_filters, selected_values=[strand_bias_filters[0]]),
                          get_div_for_the_filters(filter_name="gaps", filter_description="Gaps filter (pandora only)", filter_values=gaps_filters, selected_values=[gaps_filters[0]]),
                          html.Div([dcc.Graph(id="my-graph", style={'height': '1000px', 'width': '1000px'})]),
                          ], className="container")

    # pre-compute all traces
    config_to_all_traces = compute_all_possible_traces(df, tools, dataset_coverages, coverage_filters, strand_bias_filters, gaps_filters)


    # register the update callback
    @dash_app.callback(
        Output('my-graph', 'figure'),
        [Input(component_id='tool_checklist', component_property='value'),
         Input(component_id='dataset_checklist', component_property='value'),
         Input(component_id='coverage_checklist', component_property='value'),
         Input(component_id='strand_bias_checklist', component_property='value'),
         Input(component_id='gaps_checklist', component_property='value')])
    def update_figure(tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values, strand_bias_checklist_values, gaps_checklist_values):
        traces = []
        highest_error_rate = 0.01
        for tool in tool_checklist_values:
            if tool=="snippy_all":
                for config, trace in config_to_all_traces.items():
                    if isinstance(config, str) and config.startswith("snippy"):
                        trace = config_to_all_traces[config]
                        traces.append(trace["trace"])
                        highest_error_rate = max(highest_error_rate, trace["highest_error_rate"])
            elif tool.startswith("pandora"):
                for dataset_coverage in dataset_coverage_checklist_values:
                    for coverage_threshold in coverage_checklist_values:
                        for strand_bias_threshold in strand_bias_checklist_values:
                            for gaps_threshold in gaps_checklist_values:
                                trace = config_to_all_traces[(tool, dataset_coverage, coverage_threshold, strand_bias_threshold, gaps_threshold)]
                                traces.append(trace["trace"])
                                highest_error_rate = max(highest_error_rate, trace["highest_error_rate"])

        return {"data": traces,
                "layout": go.Layout(title="Pandora ROC (selected data)",
                                    yaxis={"title": "Recall", "range": [0.1, 1.0]},
                                    xaxis={"title": "Error rate", "range" : [0, highest_error_rate]},
                                    legend_orientation="h")}

    return dash_app


if __name__ == '__main__':
    dash_app = create_interactive_visualisation_app(ROC_data_path ="ROC_data.tsv")
    dash_app.run_server(debug=True)