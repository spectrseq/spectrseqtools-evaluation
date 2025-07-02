import sys

sys.stderr = open(snakemake.log[0], "w")

import polars as pl
from lionelmssq.plotting import plot_prediction
from lionelmssq.prediction import Prediction

prediction = Prediction.from_files(
    sequence_path=snakemake.input.pred_seq,
    fragments_path=snakemake.input.pred_fragments,
)

true_seq = snakemake.wildcards.seq
simulation = pl.read_csv(snakemake.input.sim, separator="\t")

chart = plot_prediction(prediction, true_seq, simulation)

chart.save(snakemake.output[0])
