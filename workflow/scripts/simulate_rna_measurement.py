import sys

sys.stderr = open(snakemake.log[0], "w")

import polars as pl
import random, re
from scipy.stats import norm

# regex for separating given sequence into nucleosides
nucleoside_re = re.compile(r"\d*[ACGU]")

# input
true_sequence = nucleoside_re.findall(snakemake.wildcards.seq)
n_fragments = int(snakemake.wildcards.n_fragments)
nucleoside_masses = pl.read_csv(snakemake.input[0], separator="\t")

# TODO also simulate different tails:
phosphorus = 30.973761998
oxygen = 31.989829239
d_ribose = 150.05282342
tail_5prime = [d_ribose, oxygen, 2 * oxygen + phosphorus, oxygen]


# helper functions
def random_fragment():
    l = random.randint(0, len(true_sequence) - 1)
    r = random.randint(l + 1, len(true_sequence))
    return l, r


# sample random fragments from true sequence
fragments = pl.from_records(
    [random_fragment() for _ in range(n_fragments)],
    schema=["left", "right"],
    orient="row",
)

# calculate true masses based on known nucleoside masses
true_fragment_masses = [
    sum(
        nucleoside_masses.filter(pl.col("nucleoside") == b)
        .select(pl.col("monoisotopic_mass"))
        .item()
        for b in true_sequence[fragment["left"] : fragment["right"]]
    )
    for fragment in fragments.iter_rows(named=True)
]

# simulate observed masses with noise
# TODO include breaking at different tails
fragment_noise = norm.rvs(scale=0.5, size=n_fragments)
observed_fragment_masses = [
    max(true_mass + noise, 0.0)
    for noise, true_mass in zip(fragment_noise, true_fragment_masses)
]

# compile final dataframe
fragments = fragments.with_columns(
    (pl.col("left") == 0).alias("is_start"),
    (pl.col("right") == len(true_sequence)).alias("is_end"),
    pl.Series(true_fragment_masses).alias("true_mass"),
    pl.Series(observed_fragment_masses).alias("observed_mass"),
)

# write simulation output
fragments.write_csv(snakemake.output[0], separator="\t")
