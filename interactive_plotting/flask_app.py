import plot
dash_app = plot.create_interactive_visualisation_app(ROC_data_path = "/home/leoisl/mysite/ROC_data_gt_min_0.0_step_1.0_max_1000000000.0.tsv")
app = dash_app.server
