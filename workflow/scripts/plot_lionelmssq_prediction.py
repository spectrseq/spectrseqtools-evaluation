import sys
import polars as pl

from lionelmssq.plotting import plot_prediction
from lionelmssq.prediction import Prediction


if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        prediction = Prediction.from_files(
            sequence_path=smk.input.pred_seq,
            fragments_path=smk.input.pred_fragments,
        )

        true_seq = smk.wildcards.seq
        simulation = pl.read_csv(smk.input.sim, separator="\t")

        chart = plot_prediction(
            prediction,
            true_seq,
            simulation if smk.wildcards.modus == "simulation" else None,
        )

        chart.save(smk.output[0])


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
