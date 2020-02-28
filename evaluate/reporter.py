from typing import Iterable, TextIO
import pandas as pd
from evaluate.classifier import Classifier


class Reporter:
    def __init__(self, classifiers: Iterable[Classifier], delim: str = "\t"):
        self.classifiers = classifiers
        self.delim = delim
        self.columns = [
            "sample",
            "query_probe_header",
            "ref_probe_header",
            "classification",
        ]

    def _generate_report(self,
                         fixed_info_to_add_to_query_probe_header: str = None,
                         fixed_info_to_add_to_ref_probe_header: str = None) -> pd.DataFrame:
        report_entries = []
        for classifier in self.classifiers:
            classifications = classifier.classify()
            for classification in classifications:
                assessment = classification.assessment()

                query_probe_header = str(classification.query_probe.header)
                if fixed_info_to_add_to_query_probe_header is not None:
                    query_probe_header += fixed_info_to_add_to_query_probe_header

                ref_probe_header = str(classification.ref_probe.header)
                if fixed_info_to_add_to_ref_probe_header is not None:
                    ref_probe_header += fixed_info_to_add_to_ref_probe_header

                report_entries.append(
                    [classifier.name, query_probe_header, ref_probe_header, assessment]
                )

        return pd.DataFrame(data=report_entries, columns=self.columns)

    def save_report(self, report: pd.DataFrame, file_handle: TextIO) -> None:
        report.to_csv(file_handle, sep=self.delim, header=True, index=False)


class PrecisionReporter(Reporter):
    def generate_report(self) -> pd.DataFrame:
        return self._generate_report()


class RecallReporter(Reporter):
    def generate_report(self, ref_gt_conf: int) -> pd.DataFrame:
        return self._generate_report(fixed_info_to_add_to_ref_probe_header=f"GT_CONF={ref_gt_conf};")

