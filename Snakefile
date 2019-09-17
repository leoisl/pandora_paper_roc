import logging
import pysam
from pathlib import Path
from snakemake.utils import min_version, validate
import pandas as pd

min_version("5.1.0")

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"

validate(config, "analysis/schemas/config.schema.yaml")
samples = pd.read_csv(config["samples"])
validate(samples, "analysis/schemas/samples.schema.yaml")
samples.rename(columns={"reference": "reference_assembly"}, inplace=True)

variant_calls = pd.read_csv(config["variant_calls"])
validate(variant_calls, "analysis/schemas/variant_calls.schema.yaml")
variant_calls.rename(columns={"reference": "vcf_reference"}, inplace=True)


# ======================================================
# Global variables
# ======================================================
data = pd.merge(variant_calls, samples, on="sample_id")

files=[]
for index, row in data.iterrows():
    sample_id = row["sample_id"]
    files.extend([f"analysis/truth_probesets/{sample_id}.truth_probeset.fa"])



# ======================================================
# Rules
# ======================================================
rule all:
    input: files



rules_dir = Path("rules/")
include: str(rules_dir / "common.smk")
