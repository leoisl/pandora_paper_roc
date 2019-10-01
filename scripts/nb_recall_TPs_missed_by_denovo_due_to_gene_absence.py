# with open("genes_exclusively_in_no_denovo") as fin:
with open("/hps/nobackup/research/zi/projects/pandora_paper/analysis/100x/filter/compare_no_denovo/genes_exclusively_in_no_denovo") as fin:
    genes_exclusively_in_no_denovo = set([line.strip() for line in fin])


all_recall_files_100x_no_denovo = [
    "analysis/recall/reports/063_STEC/100x/pandora_nanopore_filter_no_denovo/063_STEC_and_CFT073.report.tsv",
    "analysis/recall/reports/063_STEC/100x/pandora_nanopore_filter_no_denovo/063_STEC_and_H131800734.report.tsv",
    "analysis/recall/reports/063_STEC/100x/pandora_nanopore_filter_no_denovo/063_STEC_and_ST38.report.tsv",
    "analysis/recall/reports/CFT073/100x/pandora_nanopore_filter_no_denovo/063_STEC_and_CFT073.report.tsv",
    "analysis/recall/reports/CFT073/100x/pandora_nanopore_filter_no_denovo/CFT073_and_H131800734.report.tsv",
    "analysis/recall/reports/CFT073/100x/pandora_nanopore_filter_no_denovo/CFT073_and_ST38.report.tsv",
    "analysis/recall/reports/H131800734/100x/pandora_nanopore_filter_no_denovo/063_STEC_and_H131800734.report.tsv",
    "analysis/recall/reports/H131800734/100x/pandora_nanopore_filter_no_denovo/CFT073_and_H131800734.report.tsv",
    "analysis/recall/reports/H131800734/100x/pandora_nanopore_filter_no_denovo/H131800734_and_ST38.report.tsv",
    "analysis/recall/reports/ST38/100x/pandora_nanopore_filter_no_denovo/063_STEC_and_ST38.report.tsv",
    "analysis/recall/reports/ST38/100x/pandora_nanopore_filter_no_denovo/CFT073_and_ST38.report.tsv",
    "analysis/recall/reports/ST38/100x/pandora_nanopore_filter_no_denovo/H131800734_and_ST38.report.tsv"
]

# all_recall_files_100x_no_denovo = ["063_STEC_and_CFT073.report.tsv"]


count_of_TP_records_in_genes_is_only_in_no_denovo = 0
count_of_TP_records = 0
for recall_filepath in all_recall_files_100x_no_denovo:
    with open(recall_filepath) as recall_file:
        for line in recall_file:
            line_split = line.strip().split("\t")
            if line_split[3] in ["primary_correct", "secondary_correct", "supplementary_correct"]:
                count_of_TP_records+=1
                ref_probe_header = line_split[2]
                gene = ref_probe_header.split(";")[0].split("=")[-1]
                if gene in genes_exclusively_in_no_denovo:
                    count_of_TP_records_in_genes_is_only_in_no_denovo += 1

print(f"count_of_TP_records_in_genes_is_only_in_no_denovo = {count_of_TP_records_in_genes_is_only_in_no_denovo}")
print(f"count_of_TP_records = {count_of_TP_records}")