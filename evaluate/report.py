from pathlib import Path
import pandas as pd
from typing import Iterable, Type, List
import logging

class Report:
    def __init__(self, dfs: Iterable[pd.DataFrame]):
        self.report = pd.concat(dfs)

    def get_confident_classifications(
        self, conf_threshold: float
    ) -> List[str or float]:
        confident_classifications = self.report.query(
            "gt_conf >= @conf_threshold"
        ).classification.to_list()
        return confident_classifications

    def get_maximum_gt_conf(self) -> float:
        return self.report["gt_conf"].max()

    def get_minimum_gt_conf(self) -> float:
        return self.report["gt_conf"].min()

    @classmethod
    def from_files(cls, paths: List[Path]) -> Type["Report"]:
        reports = (pd.read_csv(path, sep="\t", keep_default_na=False) for path in paths)
        return cls(reports)

    def __eq__(self, other: "Report"):
        return self.report.equals(other.report)

    @staticmethod
    def get_GT_conf_fast(header: str):
        if header=="":
            return 0.0

        string_with_GT_conf = header[header.index("GT_CONF=") + 8:]
        string_with_GT_conf = string_with_GT_conf[:string_with_GT_conf.index(";")]
        gt_conf = float(string_with_GT_conf)
        return gt_conf

    def _create_gt_conf_column_from(self, probe_header: str) -> None:
        self.report["gt_conf"] = self.report[probe_header].apply(
            lambda column_name: Report.get_GT_conf_fast(column_name)
        )

    def _create_good_eval_column_from(self, probe_header: str) -> None:
        self.report["good_eval"] = self.report[probe_header].apply(
            lambda column_name: column_name in ["primary_correct", "secondary_correct", "supplementary_correct"]
        )

    def save_report(self, file_handle):
        self.report.to_csv(file_handle, header=True, index=False)


class PrecisionReport(Report):
    def __init__(self, dfs: Iterable[pd.DataFrame]):
        self.report = pd.concat(dfs)
        self._create_gt_conf_column_from("query_probe_header")


class RecallReport(Report):
    def __init__(self, dfs: Iterable[pd.DataFrame], concatenate_dfs_one_by_one_keeping_only_best_mappings: bool):
        if concatenate_dfs_one_by_one_keeping_only_best_mappings:
            self._concatenate_dfs_one_by_one_keeping_only_best_mappings(dfs)
            self._create_gt_conf_column_from("ref_probe_header")
        else:
            #  normal concatenation
            super().__init__(dfs)

    @classmethod
    def from_files(cls, paths: List[Path], concatenate_dfs_one_by_one_keeping_only_best_mappings: bool) -> "RecallReport":
        reports = (pd.read_csv(path, sep="\t", keep_default_na=False) for path in paths)
        return cls(reports, concatenate_dfs_one_by_one_keeping_only_best_mappings)

    def get_number_of_truth_probes(self):
        return len(self.report)

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
        self.report = self.report[original_columns]


    def _keep_only_best_mapping_for_all_truth_probes(self) -> None:
        self._create_gt_conf_column_from("ref_probe_header")
        self._create_good_eval_column_from("classification")
        self.report = self.report.sort_values(["good_eval", "gt_conf"], ascending=False).groupby("query_probe_header", as_index=False).first()

