Comparison of Simulated Data over {{snakemake.wildcards.parameter}}

Prediction results for stratified simulation over different parameters. Default settings are a
phantom rate of 0%, 100 replicates per sequence, a maximum of 10 false positive singletons, and a
mutation rate of 10%. For the simulation, all parameters are set to their default value except the one
being evaluated

The results are divided into six categories: fully correct prediction ("identical"), cor-
rect except switches between G and 55U, which are within the allowed tolerance τprediction ("identical
(minus 55U/G)"), correct composition in the wrong order ("correct composition"), correct length
but wrong composition ("failed prediction"), wrong selected length ("wrong length"), and an empty
prediction due to insufficient evidence to either build a full sequence from the two ladders or resolve
at least one position within the given time frame of the linear optimization ("no prediction").