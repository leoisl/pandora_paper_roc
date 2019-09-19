import pandas as pd
import matplotlib.pyplot as plt

error_rate_and_recall_df = pd.read_csv(
    snakemake.input.error_rate_and_recall_file, sep="\t"
)

fig, ax = plt.subplots()

for label, group in error_rate_and_recall_df.groupby("label"):
    ax = group.plot(
        ax=ax, kind="line", x="error_rate", y="recall", label=label, linewidth=0.3
    )
ax.set_ylabel("recall")
ax.set_xlim(0.0, 0.02)
ax.set_ylim(0.0, 1.0)

plt.legend(loc="best")
# plt.legend(error_rate_and_recall_df["label"].unique(), frameon=False, loc='lower left', bbox_to_anchor=(1.02, 0))
plt.savefig(snakemake.output.plot)
