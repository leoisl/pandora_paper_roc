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
        self.report = pd.concat(dfs)
        self._create_gt_conf_column_from("query_probe_header")


class RecallReport(Report):
    def __init__(self, dfs: Iterable[pd.DataFrame]):
        self._concatenate_dfs_one_by_one_keeping_only_best_mappings(dfs)
        self._create_gt_conf_column_from("ref_probe_header")

    def get_number_of_truth_probes(self):
        return len(self.report)

    def _concatenate_dfs_one_by_one_keeping_only_best_mappings(self, dfs: Iterable[pd.DataFrame]) -> None:
        self.report = None
        for index, df in enumerate(dfs):
            print(f"RecallReport._concatenate_dfs_one_by_one_keeping_only_best_mappings: processing df {index+1}...")
            if self.report is None:
                self.report = df
                continue
            self.report = pd.concat([self.report, df], ignore_index=True)
            self._keep_only_best_mapping_for_all_truth_probes()
        self._keep_only_best_mapping_for_all_truth_probes()


    def _keep_only_best_mapping_for_all_truth_probes(self) -> None:
        self._create_gt_conf_column_from("ref_probe_header")
        truth_probe_to_best_mapping = self._get_best_mapping_for_all_truth_probes()
        self.report = pd.DataFrame(columns=self.report.columns, data=truth_probe_to_best_mapping.values())


    def _get_best_mapping_for_all_truth_probes(self):
        all_truth_probes = self._get_all_truth_probes()
        truth_probe_to_all_mappings_dfs = self._get_truth_probe_to_all_mappings_sorted_by_gt_conf()
        truth_probe_to_best_mapping = {}
        for truth_probe in all_truth_probes:
            truth_probe_to_best_mapping[truth_probe] = self._get_best_mapping_for_truth_probe(truth_probe_to_all_mappings_dfs, truth_probe)
        return truth_probe_to_best_mapping


    def _get_all_truth_probes(self):
        return self.report.query_probe_header.unique()


    def _get_truth_probe_to_all_mappings_list(self):
        truth_probe_to_all_mappings_list = defaultdict(list)
        for index, series in self.report.iterrows():
            truth_probe_to_all_mappings_list[series.query_probe_header].append(list(series))
        return truth_probe_to_all_mappings_list


    def _get_truth_probe_to_all_mappings_sorted_by_gt_conf(self):
        truth_probe_to_all_mappings_list = self._get_truth_probe_to_all_mappings_list()
        for truth_probe, info_list in truth_probe_to_all_mappings_list.items():
            truth_probe_to_all_mappings_list[truth_probe] = sorted(info_list, key=lambda row: row[-1], reverse=True)
        return truth_probe_to_all_mappings_list


    def _get_best_classification(self, all_mappings_for_the_truth_probe):
        best_classification = None

        # finds the best correct one
        for mapping in all_mappings_for_the_truth_probe:
            if mapping[3] in ["primary_correct", "secondary_correct", "supplementary_correct"]:
                best_classification = mapping # this is the best one as the list is sorted by gt_conf
                break

        if best_classification is None:
            best_classification = all_mappings_for_the_truth_probe[0] # There is no correct one, get the wrong one with highest gt_conf

        return best_classification


    def _get_best_mapping_for_truth_probe(self, truth_probe_to_all_mappings_dfs, truth_probe):
        all_mappings_for_the_truth_probe = truth_probe_to_all_mappings_dfs[truth_probe]

        # the truth probe has to be found in the df
        truth_probe_found_in_the_df = len(all_mappings_for_the_truth_probe) > 0
        assert truth_probe_found_in_the_df

        best_classification = self._get_best_classification(all_mappings_for_the_truth_probe)
        return best_classification


