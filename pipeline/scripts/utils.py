import pandas as pd

def get_concatenated_df(files, separator, fields_to_keep = None):
    dfs = [pd.read_csv(file, sep=separator) for file in files]
    concatenated_df = pd.concat(dfs, ignore_index=True)
    if fields_to_keep is not None:
        concatenated_df = concatenated_df[fields_to_keep]
    return concatenated_df


def get_sample_pairs_containing_given_sample(sample_pairs, sample):
    sample_pairs_with_sample = [pair for pair in sample_pairs if sample in pair]
    sample_pairs_with_sample_as_str = [f"{sample1}_and_{sample2}" for sample1, sample2 in sample_pairs_with_sample]
    return sample_pairs_with_sample_as_str