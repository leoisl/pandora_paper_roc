import argparse
import pandas as pd

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Filter .tsv files.')
parser.add_argument('files', nargs='+', help='Paths to the .tsv files')
parser.add_argument('--output', help='Path to the output .tsv file')
parser.add_argument('--coverage', nargs='+', help='Filter by coverage')
parser.add_argument('--coverage_threshold', nargs='+', help='Filter by coverage threshold')
parser.add_argument('--strand_bias_threshold', nargs='+', help='Filter by strand bias threshold')
parser.add_argument('--gaps_threshold', nargs='+', help='Filter by gaps threshold')
parser.add_argument('--append_tool_name', nargs='+', help='Append the given value to each tool name to identify them later')
args = parser.parse_args()

# Read the .tsv files
data = [pd.read_csv(file, delimiter='\t') for file in args.files]
for i, tsv in enumerate(data):
    append = args.append_tool_name[i]
    tsv['tool'] = tsv['tool'] + append

# concatenate the .tsv files
data = pd.concat(data)

# Apply filters
filters = {
    'coverage': args.coverage,
    'coverage_threshold': args.coverage_threshold,
    'strand_bias_threshold': args.strand_bias_threshold,
    'gaps_threshold': args.gaps_threshold
}

filtered_data = data.copy()
for column, values in filters.items():
    if values:
        filtered_data = filtered_data[filtered_data[column].isin(values)]

# Read and concatenate the .tsv files
data = pd.concat([pd.read_csv(file, delimiter='\t') for file in args.files])


# Remove the second column
filtered_data.drop(columns=filtered_data.columns[0], inplace=True)

# Reset the index
filtered_data.reset_index(drop=True, inplace=True)

# Save the filtered data as a .tsv file
output_filename = args.output if args.output else 'filtered_data.tsv'
filtered_data.to_csv(output_filename, sep='\t', index=True)

print(f"Filtered data saved to {output_filename}")
