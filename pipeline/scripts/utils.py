import pandas as pd
def get_concatenated_df(files, separator, fields_to_keep):
    dfs = [pd.read_csv(file, sep=separator) for file in files]
    concatenated_df = pd.concat(dfs, ignore_index=True)[fields_to_keep]
    return concatenated_df