import sys
import polars as pl

from spectrseqtools.common import parse_nucleosides
from spectrseqtools.plotting import plot_prediction
from spectrseqtools.prediction import Prediction


if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        prediction = Prediction.from_files(
            sequence_path=smk.input.pred_seq,
            fragments_path=smk.input.pred_fragments,
        )

        true_seq = parse_nucleosides(smk.wildcards.seq)
        # simulation = pl.read_csv(smk.input.sim, separator="\t")

        charts = plot_prediction(
            prediction,
            true_seq,
            # simulation if smk.wildcards.modus == "simulation" else None,
            None,
        )

        charts[0].save(smk.output["start"])
        charts[1].save(smk.output["end"])
        charts[2].save(smk.output["internal"])
        charts[3].save(smk.output["any"])

        (charts[0] | charts[1] | charts[2]).save(smk.output["all"])


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
