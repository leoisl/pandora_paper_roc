samples=(CFT073 H131800734 063_STEC ST38)
ref="mutated_ref_snippy_NZ_CP016497.1"

for sample in "${samples[@]}"
do
  nucmer --prefix=best_snippy_vs_${sample} ${sample}.ref_chrom.fa ${ref}_${sample}_chrom.fa
  delta-filter -q -r best_snippy_vs_${sample}.delta > best_snippy_vs_${sample}.filter
  mummerplot --png -f --large best_snippy_vs_${sample}.filter -R ${sample}.ref_chrom.fa -Q ${ref}_${sample}_chrom.fa
  grep -v "set mouse clipboardformat" out.gp > out2.gp
  gnuplot out2.gp
  dnadiff -d best_snippy_vs_${sample}.delta -p best_snippy_vs_${sample}_dnadiff
  mv out.png best_snippy_vs_${sample}.png
  rm out*
done


for sample_1 in "${samples[@]}"
do
  for sample_2 in "${samples[@]}"
  do
    echo "Running for ${sample_1} and ${sample_2}"
    nucmer --prefix=${sample_1}_${sample_2} ${sample_1}.ref_chrom.fa ${sample_2}.ref_chrom.fa
    delta-filter -q -r ${sample_1}_${sample_2}.delta > ${sample_1}_${sample_2}.filter
    mummerplot --png -f --large ${sample_1}_${sample_2}.filter -R ${sample_1}.ref_chrom.fa -Q ${sample_2}.ref_chrom.fa
    grep -v "set mouse clipboardformat" out.gp > out2.gp
    gnuplot out2.gp
    dnadiff -d ${sample_1}_${sample_2}.delta -p ${sample_1}_${sample_2}_dnadiff
    mv out.png ${sample_1}_${sample_2}.png
    rm out*
  done
done