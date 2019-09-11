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

    def generate_report(self) -> pd.DataFrame:
        report_entries = []
        for classifier in self.classifiers:
            classifications = classifier.classify()
            for classification in classifications:
                assessment = classification.assessment()
                query_probe_header = str(classification.query_probe.header)
                ref_probe_header = str(classification.ref_probe.header)
                report_entries.append(
                    [classifier.name, query_probe_header, ref_probe_header, assessment]
                )

        return pd.DataFrame(data=report_entries, columns=self.columns)

    def save(self, file_handle: TextIO) -> pd.DataFrame:
        report = self.generate_report()
        report.to_csv(file_handle, sep=self.delim, header=True, index=False)
        return report


class RecallReporter(Reporter):
    pass


class PrecisionReporter(Reporter):
    pass
