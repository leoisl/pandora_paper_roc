from abc import abstractmethod
from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import math
from typing import TextIO, List, Tuple, Deque
from collections import deque


class FixVCF:
    @abstractmethod
    def correct_sample_names(self, line: str, suffix: str) -> str:
        pass

    @abstractmethod
    def get_gt_confs(self, record: str) -> Deque[float]:
        pass

    @abstractmethod
    def set_gt_confs(self, record: str, gt_confs: Deque[float]) -> str:
        pass


    @staticmethod
    def get_header_and_record_lines(vcf_filehandler: TextIO) -> Tuple[List[str], List[str]]:
        header_lines = []
        record_lines = []
        for line in vcf_filehandler:
            line = line.strip()
            is_header = line.startswith("#")

            if is_header:
                header_lines.append(line)
            else:
                record_lines.append(line)
        return header_lines, record_lines


    def correct_headers(self, headers: List[str], suffix: str) -> List[str]:
        corrected_headers = []
        for header in headers:
            is_header_with_sample_names = header.startswith("#CHROM")
            if is_header_with_sample_names:
                header = self.correct_sample_names(header, suffix)
            corrected_headers.append(header)
        return corrected_headers


    def get_all_gt_confs(self, records: List[str]) -> Deque[float]:
        all_gt_confs = deque()
        for record in records:
            gt_confs_for_this_record = self.get_gt_confs(record)
            all_gt_confs.extend(gt_confs_for_this_record)
        return all_gt_confs


    @staticmethod
    def get_log_gt_confs(gt_confs: Deque[float]) -> Deque[float]:
        log_gt_confs = deque()
        for gt_conf in gt_confs:
            assert gt_conf >= 0.0, f"Error: gt_conf is negative: {gt_conf}"
            gt_conf += 1.0  # avoids calculating log of values between 0.0 and 1.0 (which can get exponentially small)
            log_gt_conf = math.log2(gt_conf)
            log_gt_confs.append(log_gt_conf)
        return log_gt_confs


    @staticmethod
    def get_normalized_gt_confs(gt_confs: Deque[float]) -> Deque[float]:
        normalized_gt_confs = deque()
        min_gt_conf = min(gt_confs)
        max_gt_cont = max(gt_confs)
        for gt_conf in gt_confs:
            assert gt_conf >= 0.0, f"Error: log_gt_conf is negative: {gt_conf}"
            normalized_gt_conf = (gt_conf-min_gt_conf) / (max_gt_cont-min_gt_conf)
            assert 0.0 <= normalized_gt_conf <= 1.0, f"Error: normalized_gt_conf is not between 0.0 and 1.0: {normalized_gt_conf}"
            normalized_gt_confs.append(normalized_gt_conf)
        return normalized_gt_confs


    @staticmethod
    def percentile_gt_confs(gt_confs: Deque[float]) -> Deque[float]:
        return deque([round(gt_conf*100, 1) for gt_conf in gt_confs])


    def get_log_normalized_percentiled_gt_confs(self, gt_confs: Deque[float]) -> Deque[float]:
        gt_confs = self.get_log_gt_confs(gt_confs)
        gt_confs = self.get_normalized_gt_confs(gt_confs)
        gt_confs = self.percentile_gt_confs(gt_confs)
        return gt_confs


    def correct_gt_confs(self, records: List[str], corrected_gt_confs: Deque[float]) -> List[str]:
        records_corrected = []
        for record in records:
            corrected_record = self.set_gt_confs(record, corrected_gt_confs)
            records_corrected.append(corrected_record)
        assert len(corrected_gt_confs) == 0, f"Error: in correct_gt_confs(), not all corrected_gt_confs were used: {corrected_gt_confs}"
        return records_corrected


    def correct_records(self, records: List[str]) -> List[str]:
        all_gt_confs = self.get_all_gt_confs(records)
        log_normalized_percentiled_gt_confs = self.get_log_normalized_percentiled_gt_confs(all_gt_confs)
        corrected_records = self.correct_gt_confs(records, log_normalized_percentiled_gt_confs)
        return corrected_records
