samples_names=["063_STEC", "CFT073", "Escherichia_coli_MINF_1A", "Escherichia_coli_MINF_1D", "Escherichia_coli_MINF_2E", "Escherichia_coli_MINF_7C", "Escherichia_coli_MINF_8D", "Escherichia_coli_MINF_9A", "Escherichia_coli_MSB1_1A", "Escherichia_coli_MSB1_3B", "Escherichia_coli_MSB1_3C", "Escherichia_coli_MSB1_3I", "Escherichia_coli_MSB1_4E", "Escherichia_coli_MSB1_4G", "Escherichia_coli_MSB1_4I", "Escherichia_coli_MSB1_6C", "Escherichia_coli_MSB1_6J", "Escherichia_coli_MSB1_7A", "Escherichia_coli_MSB1_7C", "Escherichia_coli_MSB1_8B", "Escherichia_coli_MSB1_8G", "Escherichia_coli_MSB1_9D", "Escherichia_coli_MSB1_9I", "Escherichia_coli_MSB2_1A", "H131800734", "ST38"]
base_path="/hps/nobackup/research/zi/projects/pandora_analysis_pipeline/data_26_way/data_for_pipeline/samples"

print("sample_id,reference,mask")
for sample_name in samples_names:
    print(f"{sample_name},{base_path}/{sample_name}/{sample_name}.ref.fa,{base_path}/{sample_name}/{sample_name}.mask.bed")