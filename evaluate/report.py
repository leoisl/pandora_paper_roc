from .classification import *
from pathlib import Path
import pandas as pd
from typing import Iterable, Type, List
from collections import defaultdict


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

    def _create_gt_conf_column_from(self, probe_header: str) -> None:
        self.report["gt_conf"] = self.report[probe_header].apply(
            lambda column_name: ProbeHeader.from_string(column_name).gt_conf
        )


class PrecisionReport(Report):
    def __init__(self, dfs: Iterable[pd.DataFrame]):
        super().__init__(dfs)
        self._create_gt_conf_column_from("query_probe_header")


class RecallReport(Report):
    def __init__(self, dfs: Iterable[pd.DataFrame]):
        self._concatenate_dfs_one_by_one_keeping_only_best_mappings(dfs)
        self._create_gt_conf_column_from("ref_probe_header")

    def get_number_of_truth_probes(self):
        return len(self.report)

    def _concatenate_dfs_one_by_one_keeping_only_best_mappings(self, dfs: Iterable[pd.DataFrame]) -> None:
        self.report = None
        for df in dfs:
            if self.report is None:
                self.report = df
                continue
            self.report = pd.concat([self.report, df])
            self._keep_only_best_mapping_for_all_truth_probes()


    def _keep_only_best_mapping_for_all_truth_probes(self) -> None:
        truth_probe_to_best_mapping = self._get_best_mapping_for_all_truth_probes()
        self.report = pd.DataFrame(columns=self.report.columns, data=truth_probe_to_best_mapping.values())


    def _get_best_mapping_for_all_truth_probes(self):
        all_truth_probes = self._get_all_truth_probes()
        truth_probe_to_all_mappings_dfs = self._get_truth_probe_to_all_mappings_dfs()
        truth_probe_to_best_mapping = {}
        for truth_probe in all_truth_probes:
            truth_probe_to_best_mapping[truth_probe] = self._get_best_mapping_for_truth_probe(truth_probe_to_all_mappings_dfs, truth_probe)
        return truth_probe_to_best_mapping


    def _get_all_truth_probes(self):
        return self.report.query_probe_header.unique()


    def _get_truth_probe_to_all_mappings_series(self):
        truth_probe_to_all_mappings_series = defaultdict(list)
        for index, series in self.report.iterrows():
            truth_probe_to_all_mappings_series[series.query_probe_header].append(series)
        return truth_probe_to_all_mappings_series


    def _get_truth_probe_to_all_mappings_dfs(self):
        truth_probe_to_all_mappings_series = self._get_truth_probe_to_all_mappings_series()
        truth_probe_to_all_mappings_dfs = {}
        for truth_probe, series_list in truth_probe_to_all_mappings_series.items():
            truth_probe_to_all_mappings_dfs[truth_probe] = pd.DataFrame (
                columns=self.report.columns,
                data = series_list
            ).sort_values(by=["gt_conf"], ascending=False) # this sort is necessary to select the highest gt_conf later easier
        return truth_probe_to_all_mappings_dfs


    def _get_best_mapping_for_truth_probe(self, truth_probe_to_all_mappings_dfs, truth_probe):
        all_mappings_for_the_truth_probe = truth_probe_to_all_mappings_dfs[truth_probe]

        # the truth probe has to be found in the df
        truth_probe_found_in_the_df = len(all_mappings_for_the_truth_probe) > 0
        assert truth_probe_found_in_the_df

        correct_classification_query_str = f"classification == @AlignmentAssessment.PRIMARY_CORRECT.value or " \
            f"classification == @AlignmentAssessment.SECONDARY_CORRECT.value or " \
            f"classification == @AlignmentAssessment.SUPPLEMENTARY_CORRECT.value"
        mappings_with_correct_classifications = all_mappings_for_the_truth_probe.query(correct_classification_query_str)
        if len(mappings_with_correct_classifications) > 0:
            # selects the highest gt_conf correct mapping, which is a TP
            mapping_to_return = mappings_with_correct_classifications.iloc[0].copy()
        else:
            # selects the highest gt_conf incorrect mapping, which is a FN
            mapping_to_return = all_mappings_for_the_truth_probe.iloc[0].copy()

        return mapping_to_return


