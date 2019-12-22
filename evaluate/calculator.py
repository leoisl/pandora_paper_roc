from typing import Iterable, Type, List

import pandas as pd

from .classification import *

from pathlib import Path

from collections import defaultdict

class StatisticalClassification(Enum):
    FALSE_NEGATIVE = "fn"
    FALSE_POSITIVE = "fp"
    TRUE_POSITIVE = "tp"
    TRUE_NEGATIVE = "tn"


class EmptyReportError(Exception):
    pass


class CalculatorInfo:
    def __init__(self, true_positives: float, total: float):
        self.true_positives = float(true_positives)
        self.total = float(total)

        try:
            self.ratio = true_positives / total
        except ZeroDivisionError:
            raise EmptyReportError(
                "There are not classifications to compute recall/precision on."
            )


class Calculator:
    def get_confident_classifications(
        self, conf_threshold: float
    ) -> List[str or float]:
        confident_classifications = self.report.query(
            "gt_conf >= @conf_threshold"
        ).classification.to_list()
        return confident_classifications

    def create_gt_conf_column_from(self, probe_header: str) -> None:
        self.report["gt_conf"] = self.report[probe_header].apply(
            lambda column_name: ProbeHeader.from_string(column_name).gt_conf
        )

    def get_maximum_gt_conf(self) -> float:
        return self.report["gt_conf"].max()

    def get_minimum_gt_conf(self) -> float:
        return self.report["gt_conf"].min()

    def __init__(self, reports: Iterable[pd.DataFrame]):
        self.report = pd.concat(reports)


    @classmethod
    def from_files(cls, paths: List[Path]) -> Type["Calculator"]:
        reports = [pd.read_csv(path, sep="\t", keep_default_na=False) for path in paths]
        return cls(reports)

    def __eq__(self, other: "Calculator"):
        return self.report.equals(other.report)


class RecallInfo(CalculatorInfo):
    def __init__(self, true_positives: float, number_of_truth_probes: float):
        super().__init__(true_positives, number_of_truth_probes)
        self.recall = self.ratio


class RecallCalculator(Calculator):
    def __init__(self, reports: Iterable[pd.DataFrame]):
        super().__init__(reports)
        self.report.drop_duplicates(subset=["sample", "query_probe_header", "classification"], inplace=True)
        self.create_gt_conf_column_from("ref_probe_header")
        self.report = self.get_df_with_best_mapping_for_all_truth_probes()
        self.number_of_truth_probes = len(self.report)

    @staticmethod
    def _get_all_truth_probes(df):
        return df.query_probe_header.unique()

    @staticmethod
    def _get_truth_probe_to_all_mappings_series(df):
        truth_probe_to_all_mappings_series = defaultdict(list)
        for index, series in df.iterrows():
            truth_probe_to_all_mappings_series[series.query_probe_header].append(series)
        return truth_probe_to_all_mappings_series

    @staticmethod
    def _get_truth_probe_to_all_mappings_dfs(df):
        truth_probe_to_all_mappings_series = RecallCalculator._get_truth_probe_to_all_mappings_series(df)
        truth_probe_to_all_mappings_dfs = {}
        for truth_probe, series_list in truth_probe_to_all_mappings_series.items():
            truth_probe_to_all_mappings_dfs[truth_probe] = pd.DataFrame (
                columns=df.columns,
                data = series_list
            ).sort_values(by=["gt_conf"], ascending=False) # this sort is necessary to select the highest gt_conf later easier
        return truth_probe_to_all_mappings_dfs

    @staticmethod
    def _get_best_mapping_for_truth_probe(truth_probe_to_all_mappings_dfs, truth_probe):
        all_mappings_for_the_truth_probe = truth_probe_to_all_mappings_dfs[truth_probe]

        # the truth probe has to be found in the df
        assert len(all_mappings_for_the_truth_probe) > 0

        correct_classification_query_str = f"classification == @AlignmentAssessment.PRIMARY_CORRECT.value or " \
            f"classification == @AlignmentAssessment.SECONDARY_CORRECT.value or " \
            f"classification == @AlignmentAssessment.SUPPLEMENTARY_CORRECT.value"
        mappings_with_correct_classifications = all_mappings_for_the_truth_probe.query(correct_classification_query_str)
        if len(mappings_with_correct_classifications) > 0:
            # selects the highest gt_conf correct mapping, which is a TP
            mapping_to_return = mappings_with_correct_classifications.iloc[0]
            mapping_to_return["classification"] = StatisticalClassification.TRUE_POSITIVE
        else:
            # selects the highest gt_conf incorrect mapping, which is a FN
            mapping_to_return = all_mappings_for_the_truth_probe.iloc[0]
            mapping_to_return["classification"] = StatisticalClassification.FALSE_NEGATIVE

        return mapping_to_return

    def _get_best_mapping_for_all_truth_probes(self):
        all_truth_probes = self._get_all_truth_probes(self.report)
        truth_probe_to_all_mappings_dfs = self._get_truth_probe_to_all_mappings_dfs(self.report)
        truth_probe_to_best_mapping = {}
        for truth_probe in all_truth_probes:
            truth_probe_to_best_mapping[truth_probe] = self._get_best_mapping_for_truth_probe(truth_probe_to_all_mappings_dfs, truth_probe)
        return truth_probe_to_best_mapping

    def get_df_with_best_mapping_for_all_truth_probes(self):
        truth_probe_to_best_mapping = self._get_best_mapping_for_all_truth_probes()
        return pd.DataFrame(columns=self.report.columns, data=truth_probe_to_best_mapping.values())

    def calculate_recall(self, conf_threshold: float = 0) -> RecallInfo:
        confident_classifications = self.get_confident_classifications(conf_threshold)
        counter = Counter(confident_classifications)
        true_positives = counter[StatisticalClassification.TRUE_POSITIVE]
        return RecallInfo(true_positives, self.number_of_truth_probes)


class PrecisionInfo(CalculatorInfo):
    def __init__(self, true_positives: float, number_of_calls: float):
        super().__init__(true_positives, number_of_calls)
        self.precision = self.ratio

class PrecisionCalculator(Calculator):
    def __init__(self, reports: Iterable[pd.DataFrame]):
        super().__init__(reports)
        self.create_gt_conf_column_from("query_probe_header")

    def calculate_precision(self, conf_threshold: float = 0.0) -> PrecisionInfo:
        confident_classifications = self.get_confident_classifications(conf_threshold)
        true_positives = sum(confident_classifications)
        number_of_calls = len(confident_classifications)
        return PrecisionInfo(true_positives, number_of_calls)
