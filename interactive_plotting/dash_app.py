import dash
from plot import add_visualization_page_to_dash_app
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output
from config import pages
import base64


image_404 = base64.b64encode(open("assets/images/404.png", 'rb').read())

# set up the analyses links
analyses_links = []
for config_dict in pages:
    analyses_links.append(dcc.Link(config_dict["name"], href=config_dict["name"]))
    analyses_links.append(html.Br())

def add_all_visualization_pages_to_dash_app(dash_app):
    for config_dict in pages:
        add_visualization_page_to_dash_app(dash_app, config_dict["name"], config_dict["ROC_path"])

# set up the app
external_stylesheets = ["css/style.css"]
dash_app = dash.Dash("ROC_pandora", assets_folder="assets", external_stylesheets=external_stylesheets)
dash_app.config.suppress_callback_exceptions = True
dash_app.layouts = {}
add_all_visualization_pages_to_dash_app(dash_app)


dash_app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])


@dash_app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == "/":
        return html.Div([html.H3("Available analysis:")] + analyses_links)

    for config_dict in pages:
        name = config_dict["name"]
        if pathname == f'/{name}' or pathname == f'/{name}/':
            return dash_app.layouts[name]

    return html.Div([html.Img(src='data:image/png;base64,{}'.format(image_404.decode()))])

if __name__ == '__main__':
    dash_app.run_server(debug=True)