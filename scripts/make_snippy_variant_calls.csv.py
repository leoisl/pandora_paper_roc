from collections import Counter
from pathlib import Path
import os
from shutil import copyfile
import subprocess


config = {
    "testing": False,
    "input_folder": "/hps/nobackup/research/zi/leandro/pandora1_paper/data/variant_calls/snippy/uncorrected_snippy_files",
    "output_folder": "/hps/nobackup/research/zi/leandro/pandora1_paper/data/variant_calls/snippy/corrected_snippy_files",
    "all_refs": ["063_STEC", "CFT073", "H131800734", "ST38"],
    "output_csv": "snippy_variant_calls.csv"
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


def vcf_contains_at_least_one_record(filepath):
    with open(filepath) as file:
        for line in file:
            line = line.strip()
            if len(line) > 0 and line[0] != "#":
                return True

    return False

def get_all_non_empty_vcf_filenames_in_a_folder(folder):
    all_vcf_filepaths = Path(folder).glob("*.vcf")

    all_vcf_filenames = []
    for filepath in all_vcf_filepaths:
        if vcf_contains_at_least_one_record(filepath):
            all_vcf_filenames.append(filepath.name)

    return sorted(all_vcf_filenames)

def correct_snippy_sample_name(ref, vcf_ref_uncorrected, vcf_ref_corrected, vcf_uncorrected, vcf_corrected):
    print(f"Correcting: {vcf_ref_uncorrected} -> {vcf_ref_corrected}")
    print(f"Correcting: {vcf_uncorrected} -> {vcf_corrected}")

    copyfile(vcf_ref_uncorrected, vcf_ref_corrected)

    with open("sample_file", "w") as sample_file:
        print(ref, file=sample_file)
    subprocess.run(f"bcftools reheader -s sample_file -o {vcf_corrected} {vcf_uncorrected}", check=True, shell=True)



def main():
    input_folder = config["input_folder"]
    all_refs = config["all_refs"]
    output_csv_filepath = config["output_csv"]
    output_folder = config["output_folder"]
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print(f"Getting all datasets in {all_refs} ...")
    all_vcfs = get_all_non_empty_vcf_filenames_in_a_folder(input_folder)
    all_datasets_in_all_refs = get_all_datasets_in_all_refs(all_vcfs, all_refs)
    print(f"All datasets in {all_refs}: {all_datasets_in_all_refs}")

    print(f"Creating {output_csv_filepath} ...")
    ref_and_vcfrefuncorrected_and_vcfrefcorrected_and_vcfuncorrected_and_vcfcorrected_list = []
    with open(output_csv_filepath, "w") as output_csv_file:
        for dataset in all_datasets_in_all_refs:
            for ref in all_refs:
                vcf_ref_uncorrected = f"{input_folder}/snippy_{ref}_{dataset}.ref.fa"
                vcf_ref_corrected = f"{output_folder}/snippy_{ref}_{dataset}.ref.fa"
                vcf_uncorrected = f"{input_folder}/snippy_{ref}_{dataset}.vcf"
                vcf_corrected = f"{output_folder}/snippy_{ref}_{dataset}.vcf"
                ref_and_vcfrefuncorrected_and_vcfrefcorrected_and_vcfuncorrected_and_vcfcorrected_list.append(
                    (ref, vcf_ref_uncorrected, vcf_ref_corrected, vcf_uncorrected, vcf_corrected))
                print(f"{ref},snippy_{dataset},all,{vcf_ref_corrected},{vcf_corrected}", file=output_csv_file)
    print(f"Created {output_csv_filepath} !")

    for ref, vcf_ref_uncorrected, vcf_ref_corrected, vcf_uncorrected, vcf_corrected in ref_and_vcfrefuncorrected_and_vcfrefcorrected_and_vcfuncorrected_and_vcfcorrected_list:
        correct_snippy_sample_name(ref, vcf_ref_uncorrected, vcf_ref_corrected, vcf_uncorrected, vcf_corrected)



########################################################################################################################
# TESTING
########################################################################################################################
# TODO : move this to a proper testing
def main_test():
    test_get_dataset_name_from_vcf()
    test_get_all_datasets_in_all_refs()
    test_get_all_non_empty_vcf_filenames_in_a_folder()
    test_correct_snippy_sample_name()

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


def test_get_all_non_empty_vcf_filenames_in_a_folder():
    assert ["a.vcf", "e.vcf"] == get_all_non_empty_vcf_filenames_in_a_folder("test_vcf_folder")
    print("test_get_all_non_empty_vcf_filenames_in_a_folder(): ok")


def test_correct_snippy_sample_name():
    correct_snippy_sample_name("CFT073",
                               "test_correct_sample_name_in_vcf/snippy_CFT073_CP_OK_HEADER.ref.fa",
                               "test_correct_sample_name_in_vcf/snippy_CFT073_CP_OK_HEADER.ref.corrected.fa",
                               "test_correct_sample_name_in_vcf/snippy_CFT073_CP_OK_HEADER.vcf",
                               "test_correct_sample_name_in_vcf/snippy_CFT073_CP_OK_HEADER.corrected.vcf")
    correct_snippy_sample_name("CFT073",
                               "test_correct_sample_name_in_vcf/snippy_CFT073_CP_WRONG_HEADER.ref.fa",
                               "test_correct_sample_name_in_vcf/snippy_CFT073_CP_WRONG_HEADER.ref.corrected.fa",
                               "test_correct_sample_name_in_vcf/snippy_CFT073_CP_WRONG_HEADER.vcf",
                               "test_correct_sample_name_in_vcf/snippy_CFT073_CP_WRONG_HEADER.corrected.vcf")
    print("test_correct_snippy_sample_name(): ok")




if __name__ == "__main__":
    if config["testing"]:
        main_test()
    else:
        main()
