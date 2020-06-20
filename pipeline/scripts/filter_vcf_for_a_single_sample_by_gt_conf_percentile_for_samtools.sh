gzipped_singlesample_vcf_file=$1
gt_conf_percentile=$2
singlesample_vcf_file_gt_conf_percentile_filtered=$3

bcftools view $gzipped_singlesample_vcf_file -i "QUAL>=${gt_conf_percentile}" > $singlesample_vcf_file_gt_conf_percentile_filtered