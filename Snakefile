from pathlib import Path
import itertools
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
data: pd.DataFrame = pd.merge(variant_calls, samples, on="sample_id")
data = data.set_index(["sample_id", "coverage", "tool"], drop=False)
samples = samples.set_index(["sample_id"], drop=False)

files=[]

# Common files
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    files.extend([f"analysis/variant_calls_probesets/{sample_id}/{coverage}/{tool}.variant_calls_probeset.fa"])

#Precision files
for index, row in data.iterrows():
    sample_id, coverage, tool = row["sample_id"], row["coverage"], row["tool"]
    files.extend([f"analysis/precision/variant_calls_probesets_mapped_to_refs/{sample_id}/{coverage}/{tool}/variant_calls_probeset_mapped.sam"])


# ======================================================
# Rules
# ======================================================
rule all:
    input: files

rules_dir = Path("analysis/rules/")
include: str(rules_dir / "common.smk")
#include: str(rules_dir / "recall.smk")
include: str(rules_dir / "precision.smk")
