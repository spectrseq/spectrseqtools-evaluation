Optimization over {{snakemake.wildcards.parameter}}

Prediction results for stratified simulation over different parameters.
Default settings are a 75th percentile intensity cutoff, 60s for the full
and 5s for the short solving of a linear program.
For the simulation, all parameters are set to their default value except the one
being evaluated

The results are divided into six categories: fully correct prediction ("identical"), cor-
rect except switches between G and 55U, which are within the allowed tolerance τprediction ("identical
(minus 55U/G)"), correct composition in the wrong order ("correct composition"), correct length
but wrong composition ("failed prediction"), wrong selected length ("wrong length"), and an empty
prediction due to insufficient evidence to either build a full sequence from the two ladders or resolve
at least one position within the given time frame of the linear optimization ("no prediction").