import logging
from io import StringIO
from pathlib import Path
from typing import Tuple, Dict, List

import pysam

from .bwa import BWA
from .cli import cli
from .mummer import Nucmer, DeltaFilter, ShowSnps
from .probe import ProbeHeader, Probe
from .query import Query
from .utils import strip_extensions


def generate_mummer_snps(
    reference: Path,
    query: Path,
    prefix: Path = Path("out"),
    flank_width: int = 0,
    indels: bool = False,
    print_header: bool = True,
) -> StringIO:
    logging.info("Generating MUMmer SNPs file.")

    nucmer_params = "--maxmatch"
    nucmer = Nucmer(reference, query, str(prefix), extra_params=nucmer_params)
    nucmer_result = nucmer.run()
    nucmer_result.check_returncode()

    deltafile = Path(str(prefix) + ".delta")
    deltafilter_params = "-1"
    deltafilter = DeltaFilter(deltafile, extra_params=deltafilter_params)
    deltafilter_result = deltafilter.run()
    deltafilter_result.check_returncode()

    filtered_deltafile = prefix.with_suffix(".delta1")
    _ = filtered_deltafile.write_text(deltafilter_result.stdout.decode())

    showsnps_params = "-rlTC"
    showsnps = ShowSnps(
        filtered_deltafile,
        context=flank_width,
        extra_params=showsnps_params,
        indels=indels,
        print_header=print_header,
    )
    showsnps_result = showsnps.run()
    showsnps_result.check_returncode()
    showsnps_content = showsnps_result.stdout.decode()

    snpsfile = prefix.with_suffix(".snps")
    _ = snpsfile.write_text(showsnps_content)

    logging.info("Finished generating MUMmer SNPs file.")

    return StringIO(showsnps_content)


def write_vcf_probes_to_file(
    vcf_probes: Dict[str, str], query_name: str, tempdir: Path
) -> Path:
    query_vcf_probes = vcf_probes[query_name]
    query_vcf_probes_path: Path = tempdir / f"{query_name}.query_probes.fa"
    query_vcf_probes_path.write_text(query_vcf_probes)
    logging.info(f"VCF probes written to file: {query_vcf_probes_path}")
    return query_vcf_probes_path


def map_panel_to_probes(
    panel: Path, probes: Path, output: Path = Path(), threads: int = 1
) -> Tuple[pysam.VariantHeader, List[pysam.AlignedSegment]]:
    bwa = BWA(threads)
    bwa.index(str(probes))
    stdout, stderr = bwa.align(panel.read_text())

    # write sam to file if output path given
    if output.name:
        output.write_text(stdout)

    return bwa.parse_sam_string(stdout)


def is_mapping_invalid(record: pysam.AlignedSegment) -> bool:
    return any([record.is_unmapped, record.is_secondary, record.is_supplementary])


def assess_sam_record(record: pysam.AlignedSegment) -> str:
    assessment = ""

    if record.is_unmapped:
        assessment = "unmapped"
    elif not whole_probe_maps(record):
        assessment = "partially_mapped"
    elif record.is_secondary:
        is_correct = do_probes_match(record)
        assessment = "seconday_correct" if is_correct else "secondary_incorrect"

    return assessment


def whole_probe_maps(record: pysam.AlignedSegment) -> bool:
    if record.is_unmapped:
        return False

    truth_interval = ProbeHeader.from_string(record.query_name).interval
    truth_starts_before_alignment = truth_interval.start < record.query_alignment_start
    truth_ends_after_alignment = truth_interval.end > record.query_alignment_end

    if truth_starts_before_alignment or truth_ends_after_alignment:
        return False

    return True


def do_probes_match(record: pysam.AlignedSegment) -> bool:
    truth_probe_header = ProbeHeader.from_string(record.query_name)
    truth_probe = Probe(header=truth_probe_header, full_sequence=record.query_sequence)
    truth = truth_probe.core_sequence

    vcf_probe_header = ProbeHeader.from_string(record.reference_name)
    ref_len = record.header.get_reference_length(record.reference_name)

    # todo: code up written down solution

    for query_pos, ref_pos, ref_base in record.get_aligned_pairs(with_seq=True):
        if query_pos == REF_PANEL_FLANK_WIDTH:
            return truth == ref_base

    return False


def main():
    args = cli()

    query1: Path = args.query1
    query1_name: str = strip_extensions(query1).name
    query2: Path = args.query2
    query2_name: str = strip_extensions(query2).name
    prefix: Path = args.temp / f"{query1_name}_{query2_name}"

    mummer_snps: StringIO = generate_mummer_snps(
        query1, query2, prefix, args.truth_flank, indels=args.indels
    )
    snps_df = ShowSnps.to_dataframe(mummer_snps)
    logging.info("Making truth probesets.")
    query1_truth_probes, query2_truth_probes = snps_df.get_probes()

    query1_truth_probes_path: Path = args.temp / f"{query1_name}.truth_probes.fa"
    query2_truth_probes_path: Path = args.temp / f"{query2_name}.truth_probes.fa"
    query1_truth_probes_path.write_text(query1_truth_probes)
    logging.info(
        f"{query1_name} truth probes written to: {str(query1_truth_probes_path)}"
    )
    query2_truth_probes_path.write_text(query2_truth_probes)
    logging.info(
        f"{query2_name} truth probes written to: {str(query2_truth_probes_path)}"
    )

    logging.info("Making probes for VCF")
    samples = [query1_name, query2_name]
    query_vcf = Query(
        args.vcf, args.vcf_ref, samples=samples, flank_width=args.query_flank
    )
    vcf_probes: Dict[str, str] = query_vcf.make_probes()
    query1_vcf_probes_path: Path = write_vcf_probes_to_file(
        vcf_probes, query1_name, args.temp
    )
    query2_vcf_probes_path: Path = write_vcf_probes_to_file(
        vcf_probes, query2_name, args.temp
    )
    logging.info(f"Mapping probes for {query1_name}")
    query1_sam_file = args.temp / (query1_name + ".panel_to_probes.sam")
    query1_header, query1_sam = map_panel_to_probes(
        query1_truth_probes_path,
        query1_vcf_probes_path,
        output=query1_sam_file,
        threads=args.threads,
    )
    logging.info(f"Mapping probes for {query2_name}")
    query2_sam_file = args.temp / (query2_name + ".panel_to_probes.sam")
    query2_header, query2_sam = map_panel_to_probes(
        query2_truth_probes_path,
        query2_vcf_probes_path,
        output=query2_sam_file,
        threads=args.threads,
    )
    # todo: Assess each valid SAM record

    # todo: because of the multi-mapping, need to be careful when assessing records that I only assess mappings that cover the pandora call. Reason being is that if the truth probe maps to the flank of a pandora call, but not the actuall call part of the probe, then we will just get whatever the vcf ref is, which is an unfair comparison.
    # todo: when assessing deletions I think it makes sense to asses the bases either side of the deletion site
    # todo: Write results for each SAM record


if __name__ == "__main__":
    main()
