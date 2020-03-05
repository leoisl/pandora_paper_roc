from pathlib import Path
import logging
log_level = "INFO"
logging.basicConfig(
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)
from typing import Dict
from evaluate.query import Query
from evaluate.filtered_vcf_file import FilteredVCFFile
from evaluate.vcf_filters import VCF_Filters
from evaluate.vcf import VCFFactory
from evaluate.bwa import BWA
import pysam
from evaluate.classifier import PrecisionClassifier
from evaluate.reporter import PrecisionReporter
from evaluate.calculator import PrecisionCalculator
from evaluate.report import PrecisionReport
import argparse
import subprocess

def get_args():
    parser = argparse.ArgumentParser(description='Generate precision report given a vcf, a vcf ref and a truth ref.')
    parser.add_argument('--vcf', type=str, help='Path to the vcf file.', required=True)
    parser.add_argument('--caller', type=str, help='Which caller was used (only supported: Pandora, Clockwork or Snippy).', required=True)
    parser.add_argument('--sample_id', type=str, help='The name of the sample id in the vcf (TODO: this could be inferred).', required=True)
    parser.add_argument('--vcf_ref', type=str, help='Path to the vcf ref file.', required=True)
    parser.add_argument('--truth_ref', type=str, help='Path to the truth ref file.', required=True)
    parser.add_argument('--prefix', type=str, help='output_pref', required=True)
    parser.add_argument('--flank_width', type=int, help='The length of the probe flanks.', default=150)
    parser.add_argument('--threads', type=int, help='Number of threads', default=1)


    args = parser.parse_args()
    return args

def main():
    args = get_args()
    vcf=args.vcf
    caller=args.caller
    vcf_ref=args.vcf_ref
    truth_ref=args.truth_ref
    sample_id=args.sample_id
    prefix = args.prefix
    flank_width=args.flank_width
    threads=args.threads

    assert caller in ["Pandora", "Clockwork", "Snippy"]

    if caller=="Pandora":
        VCF_creator_method = VCFFactory.create_Pandora_VCF_from_VariantRecord_and_Sample
    elif caller=="Clockwork":
        VCF_creator_method = VCFFactory.create_Clockwork_VCF_from_VariantRecord_and_Sample
    else:
        VCF_creator_method = VCFFactory.create_Snippy_VCF_from_VariantRecord_and_Sample()

    logging.info(f"Reading VCF")
    vcf_filters = VCF_Filters()
    with pysam.VariantFile(vcf) as pysam_variant_file:
        filtered_vcf_file = FilteredVCFFile(pysam_variant_file=pysam_variant_file, filters=vcf_filters, VCF_creator_method=VCF_creator_method)

    logging.info(f"Making probes")
    query_vcf = Query(
        filtered_vcf_file,
        Path(vcf_ref),
        samples=[sample_id],
        flank_width=flank_width,
    )
    vcf_probes: Dict[str, str] = query_vcf.make_probes()
    sample_vcf_probes = vcf_probes[sample_id]

    logging.info(f"Writing probes")
    Path(f"{prefix}_probes.fa").write_text(sample_vcf_probes)


    logging.info(f"Mapping probes to truth...")
    subprocess.check_call(f"bwa index {truth_ref}", shell=True)
    BWA.map_query_to_ref(
        query=Path(f"{prefix}_probes.fa"),
        ref=Path(truth_ref),
        output=Path(f"{prefix}_mapping.sam"),
        threads=threads,
    )



    logging.info("Creating classifier")
    with pysam.AlignmentFile(f"{prefix}_mapping.sam") as sam:
        records = [record for record in sam]
    classifier = PrecisionClassifier(sam=records, name=sample_id)

    logging.info("Creating reporter")
    reporter = PrecisionReporter(classifiers=[classifier])

    logging.info("Generating report")
    report = reporter.generate_report()

    logging.info("Saving report")
    with open(f"{prefix}_mapping_report.tsv", "w") as output:
        reporter.save_report(report, output)



    logging.info(f"Loading report")
    precision_report = PrecisionReport.from_files([Path(f"{prefix}_mapping_report.tsv")])

    logging.info(f"Creating calculator")
    precision_calculator = PrecisionCalculator(precision_report)

    logging.info(f"Calculating precision")
    precision_df = precision_calculator.get_precision_report(100)

    logging.info(f"Outputting precision file")
    precision_df.to_csv(f"{prefix}_precision_report.tsv", sep="\t")

if __name__=="__main__":
    main()