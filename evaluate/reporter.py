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

    def generate_report(self, fixed_GT_conf = None) -> pd.DataFrame:
        report_entries = []
        for classifier in self.classifiers:
            classifications = classifier.classify()
            for classification in classifications:
                assessment = classification.assessment()

                query_probe_header = str(classification.query_probe.header)
                query_probe_header = self._add_suffix_to_header_if_not_present(query_probe_header, "GT_CONF", fixed_GT_conf)

                ref_probe_header = str(classification.ref_probe.header)
                ref_probe_header = self._add_suffix_to_header_if_not_present(ref_probe_header, "GT_CONF", fixed_GT_conf)

                report_entries.append(
                    [classifier.name, query_probe_header, ref_probe_header, assessment]
                )

        return pd.DataFrame(data=report_entries, columns=self.columns)

    def save_report(self, report: pd.DataFrame, file_handle: TextIO) -> None:
        report.to_csv(file_handle, sep=self.delim, header=True, index=False)

    def _get_header_suffix(self, field, value):
        header_suffix = ""
        if value is not None:
            header_suffix = f"{field}={value};"
        return header_suffix

    def _add_suffix_to_header_if_not_present(self, header, field, value):
        header_suffix = self._get_header_suffix(field, value)
        if field not in header:
            header += header_suffix
        return header


class RecallReporter(Reporter):
    pass


class PrecisionReporter(Reporter):
    pass
