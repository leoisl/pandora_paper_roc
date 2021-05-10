from pathlib import Path
import pandas as pd
from typing import Iterable, List
import logging
import math

class DelimNotFoundError(Exception):
    pass
class ReturnTypeDoesNotMatchError(Exception):
    pass

class Report:
    def __init__(self, dfs: Iterable[pd.DataFrame]):
        self.report = pd.concat(dfs, ignore_index=True)

    def get_report_satisfying_confidence_threshold(
        self, conf_threshold: float
    ) -> ["RecallReport", "PrecisionReport"]:
        report_satisfying_confidence_threshold = self.report.query("GT_CONF >= @conf_threshold")
        return self.__class__([report_satisfying_confidence_threshold])

    # Note: trivial method, not tested
    def get_classifications_as_list(self) -> List[str]:
        return self.report.classification.to_list()

    def get_maximum_gt_conf(self) -> float:
        return self.report["GT_CONF"].max()

    def get_minimum_gt_conf(self) -> float:
        return self.report["GT_CONF"].min()

    @classmethod
    def from_files(cls, paths: List[Path]) -> ["PrecisionReport", "RecallReport"]:
        reports = (pd.read_csv(path, sep="\t", keep_default_na=False) for path in paths)
        return cls(reports)

    def __eq__(self, other: "Report"):
        return self.report.equals(other.report)

    @staticmethod
    def get_value_from_header_fast(header: str, field: str, return_type, value_to_return_if_not_found, delim: str = ";"):
        """
        The slower, but better version, would be using regex, and can be found at probe.py
        """
        try:
            string_with_value = header[header.index(field) + len(field) + 1:]  # "+ 1" is for the "="
        except ValueError:
            return value_to_return_if_not_found

        try:
            string_with_value = string_with_value[:string_with_value.index(delim)]
        except ValueError:
            raise DelimNotFoundError()

        try:
            return return_type(string_with_value)
        except ValueError:
            raise ReturnTypeDoesNotMatchError

    def _create_field_from_header(self, field:str, probe_header: str, field_type, default_value) -> None:
        self.report[field] = self.report[probe_header].apply(
            lambda header: Report.get_value_from_header_fast(header, field, field_type, default_value)
        )

    def _create_gt_conf_column_from(self, probe_header: str) -> None:
        self._create_field_from_header("GT_CONF", probe_header, float, 0.0)
    def _create_pangenome_variation_id_column_from(self, probe_header: str) -> None:
        self._create_field_from_header("PVID", probe_header, int, None)
    def _create_number_of_alleles_column_from(self, probe_header: str) -> None:
        self._create_field_from_header("NB_ALL", probe_header, int, None)
    def _create_allele_id_column_from(self, probe_header: str) -> None:
        self._create_field_from_header("ALL_ID", probe_header, int, None)
    def _create_number_of_different_allele_sequences_column_from(self, probe_header: str) -> None:
        self._create_field_from_header("NB_DIFF_ALL_SEQ", probe_header, int, None)
    def _create_allele_sequence_id_column_from(self, probe_header: str) -> None:
        self._create_field_from_header("ALL_SEQ_ID", probe_header, int, None)
    def _create_nb_of_samples_column_from(self, probe_header: str) -> None:
        self._create_field_from_header("NB_OF_SAMPLES", probe_header, int, None)


    def _create_good_eval_column(self) -> None:
        self.report["good_eval"] = self.report["classification"].apply(
            lambda classification: classification in ["primary_correct", "secondary_correct", "supplementary_correct"]
        )

    # Note: not tested (trivial method)
    def save_report(self, file_handle):
        self.report.to_csv(file_handle, sep="\t", header=True, index=False)


class PrecisionReport(Report):
    def __init__(self, dfs: Iterable[pd.DataFrame]):
        self.report = pd.concat(dfs)
        self._create_gt_conf_column_from("query_probe_header")


class RecallReport(Report):
    def __init__(self, dfs: Iterable[pd.DataFrame],
                 concatenate_dfs_one_by_one_keeping_only_best_mappings: bool = True):
        if concatenate_dfs_one_by_one_keeping_only_best_mappings:
            self._concatenate_dfs_one_by_one_keeping_only_best_mappings(dfs)
        else:
            #  simple concatenation
            super().__init__(dfs)
        self._create_helper_columns()
        self.assure_there_are_no_duplicated_evaluation()


    def get_report_with_a_given_nb_of_samples(
        self, nb_of_samples: int
    ) -> "RecallReport":
        report_with_a_given_nb_of_samples = self.report.query("NB_OF_SAMPLES == @nb_of_samples")
        return self.__class__([report_with_a_given_nb_of_samples])


    # TODO: should this be added to the PrecisionReport also?
    def assure_there_are_no_duplicated_evaluation(self):
        no_duplicated_evaluation = len(self.report) == self.get_number_of_truth_probes()
        assert no_duplicated_evaluation,\
            f"[FATAL ERROR @ RecallReport.assure_there_are_no_duplicated_evaluation] Expected one evaluation line per truth probe. " \
            f"Nb of evaluation lines: {len(self.report)}. " \
            f"Nb of truth probes: {self.get_number_of_truth_probes()}."


    def _create_helper_columns(self):
        self._create_gt_conf_column_from("ref_probe_header")
        self._create_pangenome_variation_id_column_from("query_probe_header")
        self._create_number_of_alleles_column_from("query_probe_header")
        self._create_allele_id_column_from("query_probe_header")
        self._create_number_of_different_allele_sequences_column_from("query_probe_header")
        self._create_allele_sequence_id_column_from("query_probe_header")
        self._create_nb_of_samples_column_from("query_probe_header")
        self._create_good_eval_column()


    @classmethod
    def from_files(cls, paths: List[Path], concatenate_dfs_one_by_one_keeping_only_best_mappings: bool = True) -> "RecallReport":
        reports = (pd.read_csv(path, sep="\t", keep_default_na=False) for path in paths)
        return cls(reports, concatenate_dfs_one_by_one_keeping_only_best_mappings)


    # Note: trivial method, not tested
    def get_number_of_truth_probes(self):
        return len(self.report["query_probe_header"].unique())
    def get_number_of_variants(self) -> int:
        return len(self.report["PVID"].unique())


    def _concatenate_dfs_one_by_one_keeping_only_best_mappings(self, dfs: Iterable[pd.DataFrame]) -> None:
        self.report = None
        original_columns = []
        for index, df in enumerate(dfs):
            logging.info(f"RecallReport._concatenate_dfs_one_by_one_keeping_only_best_mappings: processing df {index+1}...")
            if self.report is None:
                self.report = df
                original_columns = self.report.columns
                continue
            self.report = pd.concat([self.report, df], ignore_index=True)
            self._keep_only_best_mapping_for_all_truth_probes()
        self._keep_only_best_mapping_for_all_truth_probes()
        self.report = self.report[original_columns].copy(deep=True)


    def _keep_only_best_mapping_for_all_truth_probes(self) -> None:
        self._create_gt_conf_column_from("ref_probe_header")
        self._create_good_eval_column()
        self.report = self.report.sort_values(["good_eval", "GT_CONF"], ascending=False).groupby("query_probe_header", as_index=False).first()



    def _get_id_to_nb_of_found_objects(self, object_to_find: str) -> pd.DataFrame:
        id_to_objects_found = self.report[["PVID", object_to_find, "good_eval"]].copy(deep=True)
        id_to_objects_found.drop_duplicates(inplace=True, ignore_index=True)
        nb_objects_found_column_name = f"NB_OF_{object_to_find}_FOUND"
        id_to_objects_found = id_to_objects_found.groupby(by="PVID")\
            .sum()\
            .rename(columns={"good_eval": nb_objects_found_column_name})
        id_to_objects_found = id_to_objects_found[[nb_objects_found_column_name]].copy(deep=True)
        id_to_objects_found = id_to_objects_found.astype({nb_objects_found_column_name: 'int'})
        assert len(id_to_objects_found) == self.get_number_of_variants()
        return id_to_objects_found
    def _get_id_to_nb_of_allele_sequences_found(self) -> pd.DataFrame:
        return self._get_id_to_nb_of_found_objects(object_to_find="ALL_SEQ_ID")
    def _get_id_to_nb_of_alleles_found(self) -> pd.DataFrame:
        return self._get_id_to_nb_of_found_objects(object_to_find="ALL_ID")



    def _get_id_to_total_nb_of_objects(self, field_containing_total_nb_of_objects: str) -> pd.DataFrame:
        id_to_total = self.report[["PVID", field_containing_total_nb_of_objects]].copy(deep=True)
        id_to_total.drop_duplicates(inplace=True, ignore_index=True)
        id_to_total.sort_values(by="PVID", inplace=True)
        id_to_total.set_index("PVID", inplace=True)
        assert len(id_to_total) == self.get_number_of_variants()
        return id_to_total
    def _get_id_to_nb_of_different_allele_sequences(self) -> pd.DataFrame:
        return self._get_id_to_total_nb_of_objects(
            field_containing_total_nb_of_objects="NB_DIFF_ALL_SEQ")
    def _get_id_to_nb_of_alleles(self) -> pd.DataFrame:
        return self._get_id_to_total_nb_of_objects(
            field_containing_total_nb_of_objects="NB_ALL")
    def _get_id_to_nb_of_samples(self) -> pd.DataFrame:
        return self._get_id_to_total_nb_of_objects(
            field_containing_total_nb_of_objects="NB_OF_SAMPLES")


    def _get_df_with_recall_proportions(self) -> pd.DataFrame:
        id_to_nb_of_allele_sequences_found = self._get_id_to_nb_of_allele_sequences_found()
        id_to_nb_of_different_allele_sequences = self._get_id_to_nb_of_different_allele_sequences()
        allele_seqs_merged_data = id_to_nb_of_allele_sequences_found.merge(
            id_to_nb_of_different_allele_sequences, left_index=True, right_index=True)
        allele_seqs_merged_data["proportion_of_allele_seqs_found"] = \
            allele_seqs_merged_data["NB_OF_ALL_SEQ_ID_FOUND"] / allele_seqs_merged_data["NB_DIFF_ALL_SEQ"]
        allele_seqs_merged_data["proportion_of_allele_seqs_found_binary"] = \
            allele_seqs_merged_data["proportion_of_allele_seqs_found"].apply(
                lambda proportion: int(math.isclose(proportion, 1.0)))
        id_to_nb_of_alleles_found = self._get_id_to_nb_of_alleles_found()
        id_to_nb_of_alleles = self._get_id_to_nb_of_alleles()
        alleles_merged_data = id_to_nb_of_alleles_found.merge(
            id_to_nb_of_alleles, left_index=True, right_index=True)
        alleles_merged_data["proportion_of_alleles_found"] = \
            alleles_merged_data["NB_OF_ALL_ID_FOUND"] / alleles_merged_data["NB_ALL"]

        merged_data = allele_seqs_merged_data.merge(alleles_merged_data, left_index=True, right_index=True)
        return merged_data


    def _get_df_with_recall_proportions_and_with_nb_of_samples(self) -> pd.DataFrame:
        merged_data = self._get_df_with_recall_proportions()
        id_to_nb_of_samples = self._get_id_to_nb_of_samples()
        merged_data = merged_data.merge(id_to_nb_of_samples, left_index=True, right_index=True)
        return merged_data

    def get_proportion_of_allele_seqs_found_for_each_variant(self, binary) -> List[float]:
        """
        :param binary: binarise the proportion_of_allele_seqs_found: everything < 1 becomes 0, and 1s keep being 1s.
        """
        merged_data = self._get_df_with_recall_proportions()
        column = "proportion_of_allele_seqs_found" if not binary else "proportion_of_allele_seqs_found_binary"
        return merged_data[column].to_list()

    def get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples(self, binary) -> pd.DataFrame:
        """
        :param binary: binarise the proportion_of_allele_seqs_found: everything < 1 becomes 0, and 1s keep being 1s.
        """
        merged_data = self._get_df_with_recall_proportions_and_with_nb_of_samples()
        column = "proportion_of_allele_seqs_found" if not binary else "proportion_of_allele_seqs_found_binary"
        return merged_data[[column, "NB_OF_SAMPLES"]]


    def get_proportion_of_alleles_found_for_each_variant(self) -> List[float]:
        merged_data = self._get_df_with_recall_proportions()
        return merged_data["proportion_of_alleles_found"].to_list()


    def get_proportion_of_alleles_found_for_each_variant_with_nb_of_samples(self) -> pd.DataFrame:
        merged_data = self._get_df_with_recall_proportions_and_with_nb_of_samples()
        return merged_data[["proportion_of_alleles_found", "NB_OF_SAMPLES"]]

