import sys

sys.stderr = open(snakemake.log[0], "w")

import polars as pl
import random, re
from scipy.stats import norm, uniform

# regex for separating given sequence into nucleosides
nucleoside_re = re.compile(r"\d*[ACGU]")

# input
true_sequence = nucleoside_re.findall(snakemake.wildcards.seq)
n_fragments = int(snakemake.wildcards.n_fragments)
nucleoside_masses = pl.read_csv(snakemake.input[0], separator="\t")

max_n_parts = int(snakemake.config["fragmentation_parameters"][0]["max_n_parts"])
rel_error_rate = snakemake.config["fragmentation_parameters"][1]["rel_error_rate"]
breakageline = snakemake.config["fragmentation_parameters"][2]["breakage_line"]
fragmentation_process = snakemake.config["fragmentation_parameters"][3][
    "fragmentation_process"
]
noise_distrubution = snakemake.config["fragmentation_parameters"][4][
    "noise_distrubution"
]


backbone_masses = snakemake.config["backbone_weights"]
backbone_masses = {
    k: v for d in backbone_masses for k, v in d.items()
}  # Convert the list of dictionaries to a single dictionary.

# helper functions
# def random_fragment():
#    l = random.randint(0, len(true_sequence) - 1)
#    r = random.randint(l + 1, len(true_sequence))
#    return l, r

# Divide the sequence into a given number of "max_n_parts".


def random_fragment(num_parts=max_n_parts):
    if num_parts == 1:
        # If it breaks into one part, the entire sequence remains intact.
        breakagepoints = [int(0)]

    elif num_parts <= len(true_sequence):
        # If the number of parts is less than or equal to the sequence length, then the sequence gets broken into the number of parts.
        breakagepoints = sorted(
            random.sample(range(1, len(true_sequence)), num_parts - 1)
        )
        # Change this 1 to 2 if there are not enough statistics for the sequence brekage!

    elif num_parts > len(true_sequence):
        # The number of parts cannot be greater than the sequence length.
        raise ValueError(
            "The number of parts cannot be greater than the sequence length!"
        )

    return breakagepoints


# TODO: Implement that in some cases there is no base pair generated, but only the backbone with sugar etc?

# Select how many parts the sequence gets broken into.
if fragmentation_process == "random":
    # At high energies the sequence can be divided into multiple parts! Sample uniformly how many parts it gets broken into.
    # Each fragmentation process (n_fragments) results in (1,max_n_parts) number of fragments (uniformly).
    # It can happen that the entire sequence remains intact.
    breakagepoints = [
        random_fragment(num_parts=random.randint(1, max_n_parts))
        for _ in range(n_fragments)
    ]
elif fragmentation_process == "exact":
    # Each fragmentation process (n_fragments) results in exactly max_n_parts number of fragments.
    breakagepoints = [
        random_fragment(num_parts=max_n_parts) for _ in range(n_fragments)
    ]


# Generate "left" and "right" neuclotides based on the breakagepoints,
# which are the exact indices of the nucleotides in the generated fragments:


def generate_left_right(breakagepoints):
    left, right = [], []

    for exp in breakagepoints:
        left.append(0)

        for part in exp:
            if part != 0:
                right.append(part)
                left.append(part)
        right.append(len(true_sequence))

    return left, right


l_breakage, r_breakage = generate_left_right(breakagepoints)

# sample random fragments from true sequence
fragments = pl.from_records(
    list(zip(l_breakage, r_breakage)),
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

# Calculate observed masses WITHOUT the backbone masses
if noise_distrubution == "normal":
    fragment_noise = norm.rvs(scale=rel_error_rate, size=len(true_fragment_masses))
elif noise_distrubution == "uniform":
    fragment_noise = (
        uniform.rvs(scale=rel_error_rate, size=len(true_fragment_masses))
        - rel_error_rate / 2
    )
observed_fragment_masses_without_backbone = [
    max(mass + mass * noise, 0.0)
    for noise, mass in zip(fragment_noise, true_fragment_masses)
]


# Add backbone masses to the fragments, based on the position of the fragment in the sequence!
def add_backbone_masses(fragment, backbone_masses, fragment_masses, breakageline="c-y"):
    if breakageline == "c-y":
        # TODO: Implement breaking at other points, this is currently c-y breakgae
        # Middle -> + PO2 - H
        # 5'end (start)  -> correect (-H) (We add an additional H at the end)
        # 3'end (end) cyclic phosphate -> + PO2  - 2H
        backbone_middle_addition = (
            backbone_masses["phosphorous"]
            + 2 * backbone_masses["oxygen"]
            - backbone_masses["hydrogen"]
        )
        backbone_start_addition = -backbone_masses["hydrogen"]
        backbone_end_addition = (
            backbone_masses["phosphorous"]
            + 2 * backbone_masses["oxygen"]
            - 2 * backbone_masses["hydrogen"]
        )  # Cyclic phosphate

    elif breakageline == "a-w":
        raise NotImplementedError("Breaking at 'a-w' is not implemented yet.")

    elif breakageline == "b-x":
        raise NotImplementedError("Breaking at 'a-w' is not implemented yet.")

    elif breakageline == "d-z":
        raise NotImplementedError("Breaking at 'a-w' is not implemented yet.")

    for iter, fragment in enumerate(fragments.iter_rows(named=True)):
        if fragment["left"] == 0:
            fragment_masses[iter] += (
                backbone_start_addition + backbone_masses["Tag5prime"]
            )  # Can add a 5'tag mass if applicable here!
        else:
            fragment_masses[iter] += (
                backbone_end_addition  # The end here means that its the end of a fragment which is NOT the 5' end of the sequence itself.
            )

        if fragment["right"] == len(true_sequence):
            fragment_masses[iter] += +backbone_masses[
                "Tag3prime"
            ]  # Can add a 3'tag mass if applicable here!
            # Note that if the other end of one of the fragment corressponds to the 3' end, then there is no special mass addition to be considered.

        fragment_masses[iter] += backbone_middle_addition * max(
            fragment["right"] - fragment["left"] - 1, 0
        )

        # Add a terminal Hydrogen
        fragment_masses[iter] += backbone_masses["hydrogen"]

    return fragment_masses


# Copying the fragment masses to a new list to add the backbone masses to the fragments.
true_fragment_masses_with_backbone = [elem for elem in true_fragment_masses]
true_fragment_masses_with_backbone = add_backbone_masses(
    fragments, backbone_masses, true_fragment_masses_with_backbone
)


# simulate observed masses with noise -- we use the relative error rate to simulate the noise!
if noise_distrubution == "normal":
    fragment_noise = norm.rvs(
        scale=rel_error_rate, size=len(true_fragment_masses_with_backbone)
    )
elif noise_distrubution == "uniform":
    fragment_noise = (
        uniform.rvs(scale=rel_error_rate, size=len(true_fragment_masses_with_backbone))
        - rel_error_rate / 2
    )
observed_fragment_masses = [
    max(mass + mass * noise, 0.0)
    for noise, mass in zip(fragment_noise, true_fragment_masses_with_backbone)
]

# Get the fragment sequences
fragment_sequences = [
    "".join(true_sequence[fragment["left"] : fragment["right"]])
    for fragment in fragments.iter_rows(named=True)
]

# compile final dataframe
fragments = fragments.with_columns(
    # (pl.col("left") == 0).alias("is_start"),
    # (pl.col("right") == (len(true_sequence))).alias("is_end"),
    ((pl.col("left") == 0) & (~(pl.col("right") == (len(true_sequence))))).alias("is_start"),
    ((pl.col("right") == (len(true_sequence))) & (~(pl.col("left") == 0))).alias("is_end"),
    ((pl.col("left") == 0) & (pl.col("right") == (len(true_sequence)))).alias("is_start_end"),
    ((~(pl.col("left") == 0)) & (~(pl.col("right") == (len(true_sequence))))).alias("is_internal"),
    (pl.col("right") == pl.col("left") + 1).alias("single_nucleoside"),
    pl.Series(fragment_sequences).alias("sequence"),
    pl.Series(true_fragment_masses).alias("true_nucleoside_mass"),
    pl.Series(observed_fragment_masses_without_backbone).alias(
        "observed_mass_without_backbone"
    ),
    pl.Series(true_fragment_masses_with_backbone).alias("true_mass_with_backbone"),
    pl.Series(observed_fragment_masses).alias("observed_mass"),
)

# write simulation output
fragments.write_csv(snakemake.output[0], separator="\t")
