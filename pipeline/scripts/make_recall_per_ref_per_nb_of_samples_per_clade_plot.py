from pathlib import Path
import sys
sys.path.append(str(Path().absolute()))
import logging
log_level = "INFO"
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=log_level,
    format="[%(asctime)s]:%(levelname)s: %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S %p",
)

from snakemake import shell

csv_data = snakemake.input.csv_data
plots = snakemake.output.plots
gif = snakemake.output.gif
list_with_number_of_samples = snakemake.params.list_with_number_of_samples
tools=snakemake.wildcards.tools_to_keep


# call R to generate the several plots
for input_csv, output_plot, nb_of_samples in zip(csv_data, plots, list_with_number_of_samples):
    args = f"{input_csv} 0.0 recall_{tools}_number_of_samples={nb_of_samples} {output_plot}"
    shell(f"Rscript eda/recall_per_ref_per_nb_of_samples_per_clade/clade_plots.R {args}")

# generate the gif
import imageio
images = [imageio.imread(plot) for plot in plots]
imageio.mimsave(gif, images, duration=2)
