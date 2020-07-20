from .report import Report, PrecisionReport, RecallReport
from collections import Counter
from enum import Enum
from .classification import AlignmentAssessment
import pandas as pd
from typing import Tuple
import math


class EmptyReportError(Exception):
    pass

class PrecisionInfo:
    def __init__(self, true_positives: float, number_of_calls: float):
        self.true_positives = float(true_positives)
        self.total = float(number_of_calls)
        try:
            self.precision = self.true_positives / self.total
        except ZeroDivisionError:
            raise EmptyReportError(
                "There are not classifications to compute precision on."
            )

class RecallInfo:
    def __init__(self, truth_probes_true_positives: float, truth_probes_total: float,
                       nb_variants_where_all_allele_seqs_were_found: float,
                       nb_variants_found_wrt_alleles: float, variants_total: float):
        self.truth_probes_true_positives = truth_probes_true_positives
        self.truth_probes_total = truth_probes_total
        self.nb_variants_where_all_allele_seqs_were_found = nb_variants_where_all_allele_seqs_were_found
        self.nb_variants_found_wrt_alleles = nb_variants_found_wrt_alleles
        self.variants_total = variants_total

        try:
            self.recall_wrt_truth_probes = self.truth_probes_true_positives / self.truth_probes_total
            self.recall_wrt_variants_where_all_allele_seqs_were_found = \
                self.nb_variants_where_all_allele_seqs_were_found / self.variants_total
            self.recall_wrt_variants_found_wrt_alleles = self.nb_variants_found_wrt_alleles / self.variants_total
        except ZeroDivisionError:
            raise EmptyReportError(
                "There are not classifications to compute recall on."
            )


class Calculator:
    def __init__(self, report: Report):
        self.report = report


class PrecisionCalculator(Calculator):
    def __init__(self, report: PrecisionReport):
        super().__init__(report)


    # Note: not tested, this is just data gathering
    def get_precision_report(self, all_gts) -> pd.DataFrame:
        gts = []
        precisions = []
        error_rates = []
        nb_of_correct_calls = []
        nb_of_total_calls = []
        for gt in all_gts:
            try:
                precision_info = self._calculate_precision_for_a_given_confidence(gt)
                gts.append(gt)
                precisions.append(precision_info.precision)
                error_rates.append(1 - precision_info.precision)
                nb_of_correct_calls.append(precision_info.true_positives)
                nb_of_total_calls.append(precision_info.total)
            except EmptyReportError:
                pass

        precision_df = pd.DataFrame(
            data={
                "GT": gts,
                "step_GT": list(range(len(gts))),
                "precision": precisions,
                "error_rate": error_rates,
                "nb_of_correct_calls": nb_of_correct_calls,
                "nb_of_total_calls": nb_of_total_calls
            }
        )

        return precision_df

    def _calculate_precision_for_a_given_confidence(self, conf_threshold: float = 0.0) -> PrecisionInfo:
        report_satisfying_confidence_threshold = self.report.get_report_satisfying_confidence_threshold(conf_threshold)
        confident_classifications = report_satisfying_confidence_threshold.get_classifications_as_list()
        true_positives = sum(confident_classifications)
        number_of_calls = len(confident_classifications)
        return PrecisionInfo(true_positives, number_of_calls)


class RecallCalculator(Calculator):
    def __init__(self, report: RecallReport):
        super().__init__(report)


    # Note: not tested, this is just data gathering
    def get_recall_report(self, all_gts) -> pd.DataFrame:
        gts = []
        recalls_wrt_truth_probes = []
        nbs_of_truth_probes_found = []
        nbs_of_truth_probes_in_total = []
        recalls_wrt_variants_where_all_allele_seqs_were_found = []
        recalls_wrt_variants_found_wrt_alleles = []
        nbs_variants_where_all_allele_seqs_were_found = []
        nbs_variants_found_wrt_alleles = []
        nbs_variants_total = []
        for gt in all_gts:
            try:
                recall_info = self._calculate_recall_for_a_given_confidence(gt)
                gts.append(gt)

                recalls_wrt_truth_probes.append(recall_info.recall_wrt_truth_probes)
                nbs_of_truth_probes_found.append(recall_info.truth_probes_true_positives)
                nbs_of_truth_probes_in_total.append(recall_info.truth_probes_total)

                recalls_wrt_variants_where_all_allele_seqs_were_found.append(
                    recall_info.recall_wrt_variants_where_all_allele_seqs_were_found)
                nbs_variants_where_all_allele_seqs_were_found.append(
                    recall_info.nb_variants_where_all_allele_seqs_were_found)
                recalls_wrt_variants_found_wrt_alleles.append(recall_info.recall_wrt_variants_found_wrt_alleles)
                nbs_variants_found_wrt_alleles.append(recall_info.nb_variants_found_wrt_alleles)
                nbs_variants_total.append(recall_info.variants_total)

            except EmptyReportError:
                pass

        recall_df = pd.DataFrame(
            data={
                "GT": gts,
                "step_GT": list(range(len(gts))),
                "recalls_wrt_truth_probes": recalls_wrt_truth_probes,
                "nbs_of_truth_probes_found": nbs_of_truth_probes_found,
                "nbs_of_truth_probes_in_total": nbs_of_truth_probes_in_total,
                "recalls_wrt_variants_where_all_allele_seqs_were_found": recalls_wrt_variants_where_all_allele_seqs_were_found,
                "recalls_wrt_variants_found_wrt_alleles": recalls_wrt_variants_found_wrt_alleles,
                "nbs_variants_where_all_allele_seqs_were_found": nbs_variants_where_all_allele_seqs_were_found,
                "nbs_variants_found_wrt_alleles": nbs_variants_found_wrt_alleles,
                "nbs_variants_total": nbs_variants_total,
            }
        )

        return recall_df

    # Note: not tested, this is a method called to build a specific plot
    def get_recall_report_wrt_truth_probes_for_those_present_in_a_given_nb_of_samples(self, list_of_nb_of_samples) -> pd.DataFrame:
        recalls_wrt_truth_probes = []
        nbs_of_truth_probes_found = []
        nbs_of_truth_probes_in_total = []
        for nb_of_samples in list_of_nb_of_samples:
            try:
                report_for_the_given_nb_of_samples = self.report.get_report_with_a_given_nb_of_samples(nb_of_samples)

                truth_probes_true_positives, truth_probes_total = self._calculate_info_wrt_truth_probes(
                    report_for_the_given_nb_of_samples)

                recall_wrt_truth_probes = truth_probes_true_positives / truth_probes_total

                recalls_wrt_truth_probes.append(recall_wrt_truth_probes)
                nbs_of_truth_probes_found.append(truth_probes_true_positives)
                nbs_of_truth_probes_in_total.append(truth_probes_total)

            except ZeroDivisionError:
                pass

        recall_df = pd.DataFrame(
            data={
                "GT": [0.0] * len(list_of_nb_of_samples),
                "step_GT": [0] * len(list_of_nb_of_samples),
                "nb_of_samples": list_of_nb_of_samples,
                "recalls_wrt_truth_probes": recalls_wrt_truth_probes,
                "nbs_of_truth_probes_found": nbs_of_truth_probes_found,
                "nbs_of_truth_probes_in_total": nbs_of_truth_probes_in_total,

            }
        )

        return recall_df

    def get_recall_vs_nb_of_samples_report(self, list_with_nb_of_samples) -> pd.DataFrame:
        df_with_all_nb_of_samples = self.report.get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples(
            binary=True
        )
        df_with_all_nb_of_samples = df_with_all_nb_of_samples.groupby(by="NB_OF_SAMPLES", as_index=False).mean() \
            .rename(columns={"proportion_of_allele_seqs_found_binary": "recall"})

        # get a subset of df_with_all_nb_of_samples
        recalls = []
        for nb_of_samples in list_with_nb_of_samples:
            if nb_of_samples in df_with_all_nb_of_samples["NB_OF_SAMPLES"].to_list():
                recall = df_with_all_nb_of_samples[df_with_all_nb_of_samples.NB_OF_SAMPLES == nb_of_samples]["recall"].to_list()[0]
                recalls.append(recall)
            else:
                recalls.append(0.0)

        df = pd.DataFrame({"NB_OF_SAMPLES": list_with_nb_of_samples, "recall": recalls})
        df = df.astype({"NB_OF_SAMPLES": int, "recall": float})
        return df


    @staticmethod
    def _calculate_info_wrt_truth_probes(report: RecallReport) -> Tuple[float, float]:
        classifications = report.get_classifications_as_list()
        counter = Counter(
            [
                RecallCalculator._statistical_classification(classification)
                for classification in classifications
            ]
        )
        true_positives = counter[RecallCalculator.StatisticalClassification.TRUE_POSITIVE]
        return true_positives, report.get_number_of_truth_probes()

    @staticmethod
    def _calculate_info_wrt_variants(report: RecallReport) -> Tuple[float, float, float]:
        proportions_of_allele_seqs_found = report.get_proportion_of_allele_seqs_found_for_each_variant(binary=True)
        nb_variants_where_all_allele_seqs_were_found = sum(proportions_of_allele_seqs_found)

        proportions_of_alleles_found = report.get_proportion_of_alleles_found_for_each_variant()
        nb_variants_found_wrt_alleles = sum(proportions_of_alleles_found)

        variants_total = report.get_number_of_variants()

        return nb_variants_where_all_allele_seqs_were_found, nb_variants_found_wrt_alleles, variants_total


    def _calculate_recall_for_a_given_confidence(self, conf_threshold: float = 0) -> RecallInfo:
        report_satisfying_confidence_threshold = self.report.get_report_satisfying_confidence_threshold(conf_threshold)

        truth_probes_true_positives, truth_probes_total = self._calculate_info_wrt_truth_probes(
            report_satisfying_confidence_threshold)

        nb_variants_where_all_allele_seqs_were_found, \
        nb_variants_found_wrt_alleles, \
        variants_total = self._calculate_info_wrt_variants(report_satisfying_confidence_threshold)

        return RecallInfo(truth_probes_true_positives, truth_probes_total,
                nb_variants_where_all_allele_seqs_were_found, nb_variants_found_wrt_alleles, variants_total)

    class StatisticalClassification(Enum):
        FALSE_NEGATIVE = "fn"
        FALSE_POSITIVE = "fp"
        TRUE_POSITIVE = "tp"
        TRUE_NEGATIVE = "tn"

    @staticmethod
    def _statistical_classification(classification: str) -> StatisticalClassification:
        return {
            AlignmentAssessment.UNMAPPED: RecallCalculator.StatisticalClassification.FALSE_NEGATIVE,
            AlignmentAssessment.PARTIALLY_MAPPED: RecallCalculator.StatisticalClassification.FALSE_NEGATIVE,
            AlignmentAssessment.PRIMARY_CORRECT: RecallCalculator.StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.PRIMARY_INCORRECT: RecallCalculator.StatisticalClassification.FALSE_POSITIVE,
            AlignmentAssessment.SECONDARY_CORRECT: RecallCalculator.StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.SECONDARY_INCORRECT: RecallCalculator.StatisticalClassification.FALSE_POSITIVE,
            AlignmentAssessment.SUPPLEMENTARY_CORRECT: RecallCalculator.StatisticalClassification.TRUE_POSITIVE,
            AlignmentAssessment.SUPPLEMENTARY_INCORRECT: RecallCalculator.StatisticalClassification.FALSE_POSITIVE,
        }[AlignmentAssessment(classification)]
