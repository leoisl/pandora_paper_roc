import itertools

import dash_core_components as dcc
import dash_html_components as html
import dash_table
from dash.exceptions import PreventUpdate
import plotly.graph_objs as go

import pandas as pd
from config import names_to_ROC_paths


###################################################################################
# update checklist functions
###################################################################################
def update_tool_checklist(df):
    if df is None:
        return [], None
    else:
        tools = list(df["tool"].unique())
        if snippy_was_run(tools):
            tools = remove_snippy(tools)
            tools.append("snippy_all")
        return get_options_for_filters(tools), tools

def update_dataset_checklist(df):
    if df is None:
        return [], None
    else:
        dataset_coverages = list(df["coverage"].unique())
        dataset_coverages = remove_value_from_list(dataset_coverages, "all")
        return get_options_for_filters(dataset_coverages), dataset_coverages

def update_filters_checklist(df, filter):
    if df is None:
        return [], None
    else:
        values_filter = list(df[f"{filter}_threshold"].unique())
        values_filter = remove_value_from_list(values_filter, "Not_App")
        return get_options_for_filters(values_filter), get_first_elem_of_list_as_list_itself(values_filter)

def update_coverage_checklist(df):
    return update_filters_checklist(df, "coverage")

def update_strand_bias_checklist(df):
    return update_filters_checklist(df, "strand_bias")

def update_gaps_checklist(df):
    return update_filters_checklist(df, "gaps")
###################################################################################
###################################################################################
###################################################################################






###################################################################################
# update graph functions
###################################################################################
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

def compute_traces_in_product_of_args(df, tools, dataset_coverages, coverage_filters, strand_bias_filters, gaps_filters,
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



def get_figure_for_graph(config_to_all_traces, tool_checklist_values, dataset_coverage_checklist_values,
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


def get_df_and_check_args_for_graph(plots_value, tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values,
                             strand_bias_checklist_values, gaps_checklist_values):
    if None in [plots_value, tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values,
                             strand_bias_checklist_values, gaps_checklist_values]:
        raise PreventUpdate

    df = get_df_given_analysis_name(plots_value)
    if df is None:
        raise PreventUpdate

    return df



def get_graph_proportion (plots_value, tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values,
                             strand_bias_checklist_values, gaps_checklist_values):
    df = get_df_and_check_args_for_graph (plots_value, tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values,
                             strand_bias_checklist_values, gaps_checklist_values)

    config_to_all_traces_proportion = compute_traces_in_product_of_args(df,
                                                                        tool_checklist_values,
                                                                        dataset_coverage_checklist_values,
                                                                        coverage_checklist_values,
                                                                        strand_bias_checklist_values,
                                                                        gaps_checklist_values,
                                                                        x_label="error_rate", y_label="recall")
    return dcc.Graph(figure=get_figure_for_graph(config_to_all_traces_proportion, tool_checklist_values,
                                                 dataset_coverage_checklist_values,
                                                 coverage_checklist_values, strand_bias_checklist_values,
                                                 gaps_checklist_values,
                                                 xaxis_label="Error rate", yaxis_label="Recall", set_ranges=True),
                     style={'height': '1000px', 'width': '1000px'})



def get_graph_raw (plots_value, tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values,
                             strand_bias_checklist_values, gaps_checklist_values):
    df = get_df_and_check_args_for_graph (plots_value, tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values,
                             strand_bias_checklist_values, gaps_checklist_values)

    config_to_all_traces_raw = compute_traces_in_product_of_args(df,
                                                                        tool_checklist_values,
                                                                        dataset_coverage_checklist_values,
                                                                        coverage_checklist_values,
                                                                        strand_bias_checklist_values,
                                                                        gaps_checklist_values,
                                                                        x_label="nb_of_correct_calls", y_label="nb_of_truth_probes_found")
    return dcc.Graph(figure=get_figure_for_graph(config_to_all_traces_raw, tool_checklist_values,
                                                 dataset_coverage_checklist_values,
                                                 coverage_checklist_values, strand_bias_checklist_values,
                                                 gaps_checklist_values,
                                                 xaxis_label="Nb of correct calls", yaxis_label="Nb of truth variants found", set_ranges=False),
                     style={'height': '1000px', 'width': '1000px'})


def get_data_table_with_no_gt_conf_filter(plots_value, tool_checklist_values, dataset_coverage_checklist_values, coverage_checklist_values,
                             strand_bias_checklist_values, gaps_checklist_values):
    df = get_df_and_check_args_for_graph(plots_value, tool_checklist_values, dataset_coverage_checklist_values,
                                         coverage_checklist_values, strand_bias_checklist_values, gaps_checklist_values)
    df_with_no_gt_conf_filters = df.query("step_GT == 0")
    return dash_table.DataTable(
        columns=[{"name": column, "id": column} for column in df_with_no_gt_conf_filters.columns if "Unnamed" not in column],
        data=df_with_no_gt_conf_filters.to_dict('records'))
###################################################################################
###################################################################################
###################################################################################







###################################################################################
# misc functions
###################################################################################
def get_available_analyses():
    return [ {'label': name, 'value': name} for name in names_to_ROC_paths.keys() ]

def get_options_for_filters(filter_values):
    return [{'label': f'{threshold}', 'value': threshold} for
                           threshold in filter_values]

def get_empty_div_for_the_filters(filter_name, filter_description):
    return html.Div([html.H4(f"{filter_description}:"),
              dcc.Checklist(
                  id=f"{filter_name}_checklist",
                  options=[],
                  value=None,
                  labelStyle={'display': 'inline-block', 'margin-right': '30px'}
              )],
             className="row", style={"display": "block", "width": "60%", "margin-left": "auto",
                                     "margin-right": "auto"})

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


def get_df_given_analysis_name(page_name):
    if page_name is None:
        return None

    # load df
    ROC_data_path = names_to_ROC_paths[page_name]
    df = pd.read_csv(ROC_data_path, sep="\t")
    return df


###################################################################################
###################################################################################
###################################################################################