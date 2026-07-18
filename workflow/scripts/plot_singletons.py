import sys
from pathlib import Path

from spectrseqtools.parsers import PreprocessingOptions
from spectrseqtools.plotting.plot_singletons import plot_singletons

if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        singleton_plot = plot_singletons(
            preprocessing_options=PreprocessingOptions(
                input=Path(smk.input["raw_data"]),
                meta=Path(smk.input["meta"]),
                alphabet=Path(smk.input["alphabet"]),
            ),
            scan_dir=Path(smk.params["scan_dir"]),
        )

        singleton_plot.save(Path(smk.output["all"]))


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
