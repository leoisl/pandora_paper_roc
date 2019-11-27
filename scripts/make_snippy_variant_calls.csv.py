from collections import Counter
from pathlib import Path

config = {
    "testing": False,
    "input_folder": "/hps/nobackup/research/zi/rmcolq/paper_4_way/snippy_pick_refs/snippy",
    "all_refs": ["063_STEC", "CFT073", "H131800734", "ST38"]
}

def get_dataset_name_from_vcf(refs, vcf):
    for ref in refs:
        prefix = f"snippy_{ref}"
        if prefix in vcf:
            dataset_name = vcf[len(prefix)+1 : -4]
            return dataset_name

    raise RuntimeError(f"Dataset not found in vcf: {vcf}")


def get_all_datasets_in_all_refs(all_vcf_files, refs):
    dataset_counter = Counter([get_dataset_name_from_vcf(refs, vcf) for vcf in all_vcf_files])
    datasets_in_all_refs = [dataset for dataset, count in dataset_counter.items() if count == len(refs)]
    return datasets_in_all_refs


def get_all_vcf_filenames_in_a_folder(folder):
    all_vcf_filepaths = Path(folder).glob("*.vcf")
    all_vcf_filenames = [filename.name for filename in all_vcf_filepaths]
    return sorted(all_vcf_filenames)


def main():
    input_folder = config["input_folder"]
    all_refs = config["all_refs"]

    all_vcfs = get_all_vcf_filenames_in_a_folder(input_folder)
    all_datasets_in_all_refs = get_all_datasets_in_all_refs(all_vcfs, all_refs)

    print("sample_id,tool,coverage,reference,vcf")
    for dataset in all_datasets_in_all_refs:
        for ref in all_refs:
            print(f"{ref},snippy,all,{input_folder}/snippy_{ref}_{dataset}.ref.fa,{input_folder}/snippy_{ref}_{dataset}.vcf")




########################################################################################################################
# TESTING
########################################################################################################################
# TODO : move this to a proper testing
def main_test():
    test_get_dataset_name_from_vcf()
    test_get_all_datasets_in_all_refs()
    test_get_all_vcf_filenames_in_a_folder()

def test_get_dataset_name_from_vcf():
    all_refs = ["063_STEC", "CFT073", "H131800734", "ST38"]
    assert "NZ_CP013483.1" == get_dataset_name_from_vcf(all_refs, "snippy_063_STEC_NZ_CP013483.1.vcf")
    assert "NZ_CP013483.1" == get_dataset_name_from_vcf(all_refs, "snippy_CFT073_NZ_CP013483.1.vcf")
    assert "NZ_CP013483.1" == get_dataset_name_from_vcf(all_refs, "snippy_H131800734_NZ_CP013483.1.vcf")
    assert "NZ_CP013483.1" == get_dataset_name_from_vcf(all_refs, "snippy_ST38_NZ_CP013483.1.vcf")

    try:
        get_dataset_name_from_vcf(all_refs, "snippy_ref_does_not_exist_NZ_CP013483.1.vcf")
        assert False
    except RuntimeError:
        assert True

    print("test_get_dataset_name_from_vcf(): ok")


def test_get_all_datasets_in_all_refs():
    all_refs = ["063_STEC", "CFT073", "H131800734", "ST38"]
    all_vcf_files = ["snippy_063_STEC_NZ_CP013483.1.vcf", "snippy_CFT073_NZ_CP013483.1.vcf",
                     "snippy_H131800734_NZ_CP013483.1.vcf", "snippy_ST38_NZ_CP013483.1.vcf",
                     "snippy_063_STEC_NZ_CP013483.2.vcf", "snippy_CFT073_NZ_CP013483.2.vcf",
                     "snippy_H131800734_NZ_CP013483.2.vcf"]
    assert ["NZ_CP013483.1"] == get_all_datasets_in_all_refs(all_vcf_files, refs=all_refs)
    print("test_get_all_datasets_in_all_refs(): ok")


def test_get_all_vcf_filenames_in_a_folder():
    assert ["a.vcf", "e.vcf"] == get_all_vcf_filenames_in_a_folder("test_vcf_folder")
    print("test_get_all_vcf_filenames_in_a_folder(): ok")








if __name__ == "__main__":
    if config["testing"]:
        main_test()
    else:
        main()
