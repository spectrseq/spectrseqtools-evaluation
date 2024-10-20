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
max_n_parts = int(snakemake.wildcards.parameters.max_n_parts)

# TODO also simulate different tails:
phosphorus = 30.973761998
oxygen = 31.989829239
d_ribose = 150.05282342
tail_5prime = [d_ribose, oxygen, 2 * oxygen + phosphorus, oxygen]


# helper functions
#def random_fragment():
#    l = random.randint(0, len(true_sequence) - 1)
#    r = random.randint(l + 1, len(true_sequence))
#    return l, r

# Divide the sequence into a maximum given number of "max_n_parts". 
# At high energies the sequence can be divided into multiple parts! Sample uniformly how many parts it gets broken into.
# TODO: Implement the proper addition of backbone masses with the nucleoside masses
# TODO: Implement that in some cases there is no base pair generated, but only the backbone with sugar etc.

def random_fragment(num_parts = max_n_parts):
    if num_parts == 1: #If it breaks into one part, the entire sequence remains intact.
        breakagepoints = [int(0)]
    else:
        breakagepoints = sorted(random.sample(range(1, len(true_sequence)), num_parts - 1))
        #Implement the condition that the num_parts cannot be greater than the sequence length!
    return breakagepoints

#Select how many parts the sequence gets broken into.
#It can happen that the entire sequence remains intact.
breakagepoints = [random_fragment(num_parts = random.randint(1,max_n_parts)) for _ in range(n_fragments)]

# Generate "left" and "right" neuclotides based on the breakagepoints, 
# which are the exact indices of the nucleotides in the generated fragments:
# NOTE: The change here from before: now left and right are the exact indices of the python list of nucleotides,
# with both left and right included in the fragment. If left = right, then its a single nucleotide.
def generate_left_right(breakagepoints):
    left,right = [],[]
    for exp in breakagepoints:
        left.append(0)
        for part in exp:
            if part != 0:
                right.append(part-1)
                left.append(part)
        right.append(len(true_sequence)-1)
    return left,right

l_breakage,r_breakage = generate_left_right(breakagepoints) 

# sample random fragments from true sequence
fragments = pl.from_records(
    list(zip(l_breakage,r_breakage)), 
    schema=["left", "right"],
    orient="row",
)

# calculate true masses based on known nucleoside masses
true_fragment_masses = [
    sum(
        nucleoside_masses.filter(pl.col("nucleoside") == b)
        .select(pl.col("monoisotopic_mass"))
        .item()
        for b in true_sequence[fragment["left"] : fragment["right"] + 1]
    )
    for fragment in fragments.iter_rows(named=True)
]

# simulate observed masses with noise
# TODO include breaking at different tails
fragment_noise = norm.rvs(scale=0.5, size=len(true_fragment_masses))
observed_fragment_masses = [
    max(true_mass + noise, 0.0)
    for noise, true_mass in zip(fragment_noise, true_fragment_masses)
]

# compile final dataframe
fragments = fragments.with_columns(
    (pl.col("left") == 0).alias("is_start"),
    (pl.col("right") == (len(true_sequence))-1).alias("is_end"),
    (pl.col("right") == pl.col("left")).alias("single_nucleotide"),
    pl.Series(true_fragment_masses).alias("true_mass"),
    pl.Series(observed_fragment_masses).alias("observed_mass"),
)

# write simulation output
fragments.write_csv(snakemake.output[0], separator="\t")
