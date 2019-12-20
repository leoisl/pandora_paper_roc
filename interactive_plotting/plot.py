import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
import itertools
import dash_table


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

def get_snippy_trace_data_for_given_ref(df, snippy_ref):
    trace_name = f"Tool: {snippy_ref}"
    config = (snippy_ref, )
    df_for_snippy = df.query(f"tool==\"{snippy_ref}\"")
    trace_data = {
        "trace_name": trace_name,
        "config": config,
        "df_for_snippy": df_for_snippy
    }
    return trace_data

def compute_all_possible_traces(df, tools, dataset_coverages, coverage_filters, strand_bias_filters, gaps_filters,
                                x_label, y_label):
    config_to_trace = {}
    for tool in tools:
        if tool == "snippy_all":
            color_vector = ["red", "green", "blue", "black", "yellow", "orange", "brown", "grey"]
            snippy_refs_tools = [snippy_ref_tool for snippy_ref_tool in df["tool"].unique() if snippy_ref_tool.startswith("snippy")]
            for snippy_index, snippy_ref in enumerate(snippy_refs_tools):
                snippy_trace_data = get_snippy_trace_data_for_given_ref(df, snippy_ref)

                trace_name = snippy_trace_data["trace_name"]
                config = snippy_trace_data["config"]
                df_for_label = snippy_trace_data["df_for_snippy"]

                highest_error_rate_after_20_percent_recall = get_highest_error_rate_after_given_recall(df_for_label, 0.2)

                trace = go.Scatter(x=df_for_label[x_label], y=df_for_label[y_label], name=trace_name,
                                   mode='lines',
                                   marker={'size': 8, "opacity": 0.6, "line": {'width': 0.5}},
                                   marker_color=color_vector[snippy_index % len(color_vector)])

                config_to_trace[config] = {
                    "trace": trace,
                    "highest_error_rate": highest_error_rate_after_20_percent_recall
                }
        else:
            for dataset_coverage, coverage_threshold, strand_bias_threshold, gaps_threshold in \
                itertools.product(dataset_coverages, coverage_filters, strand_bias_filters, gaps_filters):
                df_for_label = df.query(
                    f"tool==\"{tool}\" & coverage==\"{dataset_coverage}\" & coverage_threshold==\"{coverage_threshold}\" & strand_bias_threshold==\"{strand_bias_threshold}\" & gaps_threshold==\"{gaps_threshold}\"")

                trace_name = f"Tool: {tool}, Coverage: {dataset_coverage}, Coverage threshold: {coverage_threshold}, Strand bias threshold: {strand_bias_threshold}, Gaps threshold: {gaps_threshold}"
                highest_error_rate_after_20_percent_recall = get_highest_error_rate_after_given_recall(df_for_label, 0.2)

                trace = go.Scatter(x=df_for_label[x_label], y=df_for_label[y_label], name=trace_name,
                                        mode='lines',
                                        marker={'size': 8, "opacity": 0.6, "line": {'width': 0.5}})


                config_to_trace[(tool, dataset_coverage, coverage_threshold, strand_bias_threshold, gaps_threshold)] = {
                    "trace": trace,
                    "highest_error_rate": highest_error_rate_after_20_percent_recall
                }

    return config_to_trace

def remove_value_from_list(array, value):
    return [elem for elem in array if elem != value]

def remove_snippy(tools):
    return [tool for tool in tools if "snippy" not in tool]


def snippy_was_run(tools):
    for tool in tools:
        if "snippy" in tool:
            return True
    return False

def get_first_elem_of_list_as_list_itself (the_list):
    if len(the_list) >= 1:
        return [the_list[0]]
    else:
        return []


def update_figure_core(config_to_all_traces, tool_checklist_values, dataset_coverage_checklist_values,
                       coverage_checklist_values, strand_bias_checklist_values, gaps_checklist_values,
                       xaxis_label, yaxis_label, set_ranges):
    traces = []
    highest_error_rate = 0.01
    for tool in tool_checklist_values:
        if tool == "snippy_all":
            for config, trace in config_to_all_traces.items():
                if config[0].startswith("snippy"):
                    trace = config_to_all_traces[config]
                    traces.append(trace["trace"])
                    highest_error_rate = max(highest_error_rate, trace["highest_error_rate"])
        else:
            for dataset_coverage, coverage_threshold, strand_bias_threshold, gaps_threshold in \
                    itertools.product(dataset_coverage_checklist_values, coverage_checklist_values,
                                      strand_bias_checklist_values, gaps_checklist_values):
                trace = config_to_all_traces[
                    (tool, dataset_coverage, coverage_threshold, strand_bias_threshold, gaps_threshold)]
                traces.append(trace["trace"])
                highest_error_rate = max(highest_error_rate, trace["highest_error_rate"])

    if set_ranges:
        layout = go.Layout(title="Pandora ROC (selected data)",
                                yaxis={"title": yaxis_label, "range": [0.1, 1.0]},
                                xaxis={"title": xaxis_label, "range": [0, highest_error_rate]},
                                legend_orientation="h")
    else:
        layout = go.Layout(title="Pandora ROC (selected data)",
                           yaxis={"title": yaxis_label, 'autorange': True},
                           xaxis={"title": xaxis_label, 'autorange': True},
                           legend_orientation="h")
    return {"data": traces,
            "layout": layout}


def add_visualization_page_to_dash_app(dash_app, page_name, ROC_data_path):
    # load df
    df = pd.read_csv(ROC_data_path, sep="\t")

    # load filters
    tools = list(df["tool"].unique())
    dataset_coverages = list(df["coverage"].unique())
    coverage_filters = list(df["coverage_threshold"].unique())
    strand_bias_filters = list(df["strand_bias_threshold"].unique())
    gaps_filters = list(df["gaps_threshold"].unique())

    # adjust the filters
    if snippy_was_run(tools):
        tools = remove_snippy(tools)
        tools.append("snippy_all")
    dataset_coverages = remove_value_from_list(dataset_coverages, "all")
    coverage_filters = remove_value_from_list(coverage_filters, "Not_App")
    strand_bias_filters = remove_value_from_list(strand_bias_filters, "Not_App")
    gaps_filters = remove_value_from_list(gaps_filters, "Not_App")

    # set up the layout
    df_with_no_gt_conf_filters = df.query("step_GT == 0")
    dash_app.layouts[page_name] = html.Div([
        html.Div([html.H1(f"Pandora ROC evaluation - {page_name}")], style={'textAlign': "center"}),
        get_div_for_the_filters(filter_name=f"{page_name}_tool", filter_description="Tool", filter_values=tools, selected_values=tools),
        get_div_for_the_filters(filter_name=f"{page_name}_dataset", filter_description="Dataset coverage filter (pandora only)", filter_values=dataset_coverages, selected_values=dataset_coverages),
        get_div_for_the_filters(filter_name=f"{page_name}_coverage", filter_description="Coverage filter (pandora only)", filter_values=coverage_filters, selected_values=get_first_elem_of_list_as_list_itself(coverage_filters)),
        get_div_for_the_filters(filter_name=f"{page_name}_strand_bias", filter_description="Strand bias filter (pandora only)", filter_values=strand_bias_filters, selected_values=get_first_elem_of_list_as_list_itself(strand_bias_filters)),
        get_div_for_the_filters(filter_name=f"{page_name}_gaps", filter_description="Gaps filter (pandora only)", filter_values=gaps_filters, selected_values=get_first_elem_of_list_as_list_itself(gaps_filters)),
        html.H1("Proportion ROC curve"),
        html.Div([dcc.Graph(id=f"{page_name}_graph_proportion", style={'height': '1000px', 'width': '1000px'})]),
        html.H1("Raw numbers ROC curve"),
        html.Div([dcc.Graph(id=f"{page_name}_graph_raw", style={'height': '1000px', 'width': '1000px'})]),
        html.H1("Extra informations:"),
        html.H3("Data with no GT conf filtering:"),
        html.Div(dash_table.DataTable(id="data_with_no_gt_conf_filter",
                              columns=[{"name": column, "id": column} for column in df_with_no_gt_conf_filters.columns if "Unnamed" not in column],
                             data=df_with_no_gt_conf_filters.to_dict('records'))),
        # html.H3("Whole data:"),
        # html.Div(dash_table.DataTable(id="whole_data",
        #                               columns=[{"name": column, "id": column} for column in
        #                                        df.columns if "Unnamed" not in column],
        #                               data=df.to_dict('records')))
    ], className="container")

    # pre-compute all traces
    config_to_all_traces_proportion = compute_all_possible_traces(df, tools, dataset_coverages, coverage_filters,
                                                       strand_bias_filters, gaps_filters,
                                                       x_label="error_rate", y_label="recall")
    config_to_all_traces_raw = compute_all_possible_traces(df, tools, dataset_coverages, coverage_filters,
                                                                  strand_bias_filters, gaps_filters,
                                                                  x_label="nb_of_correct_calls", y_label="nb_of_truth_probes_found")


    # register the update callback
    @dash_app.callback(
        Output(f'{page_name}_graph_proportion', 'figure'),
        [Input(component_id=f'{page_name}_tool_checklist', component_property='value'),
         Input(component_id=f'{page_name}_dataset_checklist', component_property='value'),
         Input(component_id=f'{page_name}_coverage_checklist', component_property='value'),
         Input(component_id=f'{page_name}_strand_bias_checklist', component_property='value'),
         Input(component_id=f'{page_name}_gaps_checklist', component_property='value')])
    def update_figure_proportion(tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values, strand_bias_checklist_values, gaps_checklist_values):
        return update_figure_core(config_to_all_traces_proportion, tool_checklist_values, dataset_coverage_checklist_values,
                                  coverage_checklist_values, strand_bias_checklist_values, gaps_checklist_values,
                                  xaxis_label="Error rate", yaxis_label="Recall", set_ranges=True)

    @dash_app.callback(
        Output(f'{page_name}_graph_raw', 'figure'),
        [Input(component_id=f'{page_name}_tool_checklist', component_property='value'),
         Input(component_id=f'{page_name}_dataset_checklist', component_property='value'),
         Input(component_id=f'{page_name}_coverage_checklist', component_property='value'),
         Input(component_id=f'{page_name}_strand_bias_checklist', component_property='value'),
         Input(component_id=f'{page_name}_gaps_checklist', component_property='value')])
    def update_figure_raw(tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values, strand_bias_checklist_values, gaps_checklist_values):
        return update_figure_core(config_to_all_traces_raw, tool_checklist_values, dataset_coverage_checklist_values,
                                  coverage_checklist_values, strand_bias_checklist_values, gaps_checklist_values,
                                  xaxis_label="Nb of correct calls", yaxis_label="Nb of truth variants found", set_ranges=False)