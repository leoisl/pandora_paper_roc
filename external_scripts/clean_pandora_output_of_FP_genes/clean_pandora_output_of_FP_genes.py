import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description='Clean pandora_multisample_genotyped_global.vcf of FPs.')
    parser.add_argument('--vcf', type=str, help='Path to pandora_multisample_genotyped_global.vcf', required=True)
    parser.add_argument('--classification_df_renamed', type=str, help='Gene classification csv, obtained with https://github.com/leoisl/snippy_calls_gene_distance/blob/FP_genes/notebooks/gene_presence_matrix/FP_genes.ipynb',
                        required=True)
    args = parser.parse_args()
    return args


def infer_if_gene_is_present_in_sample(gene_classification_series, sample):
    try:
        gene_classification_in_sample = gene_classification_series[sample].to_list()[0]
        return gene_classification_in_sample == "TP"
    except IndexError:
        return False



def nullify_call(line_split, sample, sample_index):
    original_call = line_split[sample_index + 9]
    nullified_call = "." + original_call[original_call.index(":"):]

    if original_call != nullified_call:
        print(f"Call nullified for:")
        print(f"line_split = {line_split}")
        print(f"sample = {sample}")
        print(f"original_call = {original_call}")
        print(f"nullified_call = {nullified_call}")
        line_split[sample_index + 9] = nullified_call


def get_nullified_calls_of_non_TP_genes(vcf_lines, gene_classification_csv, samples):
    nullified_calls = [None] * len(vcf_lines)

    for line_index, line in enumerate(vcf_lines):
        if line_index % 10000 == 0:
            print(f"{line_index} lines processed")

        line = line.strip()
        if line.startswith("#"):
            nullified_calls[line_index] = line
            continue

        line_split = line.split()
        gene = line_split[0]
        gene_classification_series = gene_classification_csv[gene_classification_csv.gene_name == gene]

        for sample_index, sample in enumerate(samples):
            gene_is_present_in_sample = infer_if_gene_is_present_in_sample(gene_classification_series, sample)
            if not gene_is_present_in_sample:
                nullify_call(line_split, sample, sample_index)

        nullified_call = "\t".join(line_split)
        nullified_calls[line_index] = nullified_call

    return nullified_calls

def get_headers(vcf_filename):
    headers=[]
    with open(vcf_filename) as vcf_filehandler:
        for line in vcf_filehandler:
            line = line.strip()
            if line.startswith("#"):
                headers.append(line)
            else:
                break
    return headers

def get_samples(header_line):
    header_line_split = header_line.split()
    return header_line_split[9:]

def process(vcf_filename, gene_classification_csv_filename):
    print("Loading...")
    gene_classification_csv = pd.read_csv(gene_classification_csv_filename)
    headers = get_headers(vcf_filename)
    samples = get_samples(headers[-1])
    # this takes a lot of RAM, but it is faster I guess
    with open(vcf_filename) as vcf_filehandler:
        vcf_lines = vcf_filehandler.readlines()

    print("Processing...")
    nullified_calls = get_nullified_calls_of_non_TP_genes(vcf_lines, gene_classification_csv, samples)

    print("Outputting...")
    with open(vcf_filename+".no_FPs.vcf", "w") as vcf_no_fps_filehandler:
        print("\n".join(nullified_calls), file=vcf_no_fps_filehandler)

    print("Done!")

def main():
    args = get_args()
    process(args.vcf, args.classification_df_renamed)


if __name__=="__main__":
    main()