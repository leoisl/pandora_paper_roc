from typing import Iterable, List

import pysam

from evaluate.classification import (
    Classification,
    RecallClassification,
    PrecisionClassification,
)


class Classifier:
    def __init__(self, sam: Iterable[pysam.AlignedSegment] = None, name: str = ""):
        if sam is None:
            sam = []
        self.sam = sam
        self.name = name

    def classify(self) -> List[Classification]:
        classifications = []
        for record in self.sam:
            classification = self.make_classification(record=record)
            classifications.append(classification)

        return classifications

    def make_classification(self, record):
        return Classification(record)


class RecallClassifier(Classifier):
    def make_classification(self, record):
        return RecallClassification(record)


class PrecisionClassifier(Classifier):
    def make_classification(self, record):
        return PrecisionClassification(record)
