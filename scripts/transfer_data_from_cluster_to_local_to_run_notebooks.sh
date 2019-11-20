PANDORA_FILES_PATH="noah-login-ext:/hps/nobackup/research/zi/projects/pandora_paper_leandro/old_analysis/100x/filter"
mkdir -p cluster/analysis
rsync --relative -zcvh ${PANDORA_FILES_PATH}/compare_with_denovo/pandora_multisample.vcf_ref.fa cluster/pandora
rsync --relative -zcvh ${PANDORA_FILES_PATH}/compare_with_denovo/pandora_multisample.matrix cluster/pandora
rsync --relative -zcvh ${PANDORA_FILES_PATH}/compare_with_denovo/pandora_multisample_genotyped.vcf cluster/pandora
rsync --relative -zcvh ${PANDORA_FILES_PATH}/compare_no_denovo/pandora_multisample_genotyped.vcf cluster/pandora
rsync --relative -zcvh ${PANDORA_FILES_PATH}/compare_no_denovo/pandora_multisample.vcf_ref.fa cluster/pandora
rsync --relative -zcvh ${PANDORA_FILES_PATH}/compare_no_denovo/pandora_multisample.matrix cluster/pandora

DATA_FILES_PATH="noah-login-ext:/hps/nobackup/research/zi/projects/pandora_paper/data"
mkdir -p cluster/data
rsync --relative -zcvh ${DATA_FILES_PATH}/H131800734/H131800734.ref.pilon.fa cluster/data
rsync --relative -zcvh ${DATA_FILES_PATH}/H131800734/mask/H131800734.mask.bed cluster/data
rsync --relative -zcvh ${DATA_FILES_PATH}/ST38/ST38.ref.pilon.fa cluster/data
rsync --relative -zcvh ${DATA_FILES_PATH}/ST38/mask/ST38.mask.bed cluster/data
rsync --relative -zcvh ${DATA_FILES_PATH}/CFT073/CFT073.ref.fa cluster/data
rsync --relative -zcvh ${DATA_FILES_PATH}/CFT073/mask/CFT073.mask.bed cluster/data
rsync --relative -zcvh ${DATA_FILES_PATH}/063_STEC/063_STEC.ref.fa cluster/data
rsync --relative -zcvh ${DATA_FILES_PATH}/063_STEC/mask/063_STEC.mask.bed cluster/data

EVALUATION_FILES_PATH="noah-login-ext:/hps/nobackup/research/zi/leandro/pandora1_paper/analysis/recall/reports"
rsync --relative -zcvhr ${EVALUATION_FILES_PATH}/ cluster/analysis



