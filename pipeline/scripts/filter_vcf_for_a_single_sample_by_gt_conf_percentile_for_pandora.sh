gzipped_singlesample_vcf_file=$1
gt_conf_percentile=$2
singlesample_vcf_file_gt_conf_percentile_filtered=$3

min_gt=`bcftools query -f '[ %GT_CONF]\n' $gzipped_singlesample_vcf_file | sort -n | head -n 1`
echo "min_gt = $min_gt"

max_gt=`bcftools query -f '[ %GT_CONF]\n' $gzipped_singlesample_vcf_file | sort -nr | head -n 1`
echo "max_gt = $max_gt"

expression=`echo "( ( $max_gt - $min_gt ) * $gt_conf_percentile * 0.01) + $min_gt"`
echo "expression = $expression"

gt_conf=`echo "$expression" | bc`
echo "gt_conf = $gt_conf"

bcftools view $gzipped_singlesample_vcf_file -i "FORMAT/GT_CONF>=${gt_conf}" > $singlesample_vcf_file_gt_conf_percentile_filtered