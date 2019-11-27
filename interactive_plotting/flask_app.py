import plot
dash_app = plot.create_interactive_visualisation_app(ROC_data_path = "/home/leoisl/mysite/ROC_data.tsv")
app = dash_app.server
