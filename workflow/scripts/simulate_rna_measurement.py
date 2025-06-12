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
nucleoside_masses = pl.read_csv(snakemake.input["nucleosides"], separator="\t")
element_masses = pl.read_csv(snakemake.input["elements"], separator="\t")

random.seed(snakemake.config["fragmentation_params"]["random_seed"])
max_n_parts = int(snakemake.config["fragmentation_params"]["max_n_parts"])
rel_error_rate = snakemake.config["fragmentation_params"]["rel_error_rate"]
breakage_line = snakemake.config["fragmentation_params"]["breakage_line"]
frag_process = snakemake.config["fragmentation_params"]["fragmentation_process"]
noise_dist = snakemake.config["fragmentation_params"]["noise_distribution"]
mass_5_prime = snakemake.config["fragmentation_params"]["mass_5_prime"]
mass_3_prime = snakemake.config["fragmentation_params"]["mass_3_prime"]

# Convert dataframe of element weights to dict
element_masses = {row[element_masses.get_column_index("symbol")]:
                  row[element_masses.get_column_index("mass")] for
                  row in element_masses.iter_rows()}

# helper functions
# def random_fragment():
#    l = random.randint(0, len(true_sequence) - 1)
#    r = random.randint(l + 1, len(true_sequence))
#    return l, r

# Divide the sequence into a given number of "max_n_parts".


def select_fragmentation_sites(num_parts=max_n_parts):
    # Ensure there is a positive number of parts
    if num_parts <= 0:
        raise ValueError(
            "The number of parts cannot be less than one!"
        )

    # Ensure the number of parts is not greater than the sequence length
    if num_parts > len(true_sequence):
        raise ValueError(
            "The number of parts cannot be greater than the sequence length!"
        )

    # If the sequence "breaks" into one part, it remains intact
    if num_parts == 1:
        return [int(0)]

    # Return randomly sampled breakage positions in the sequence
    return sorted(random.sample(range(1, len(true_sequence)), num_parts - 1))
    # LCK: I do not understand the comment below
    # Change this 1 to 2 if there are not enough statistics for the sequence breakage!


# TODO: Implement that in some cases there is no base pair generated, but only the backbone with sugar etc?

# Select how many parts the sequence gets broken into.
match frag_process:
    case "random":
        # At high energies the sequence can be divided into multiple parts! Sample uniformly how many parts it gets broken into.
        # Each fragmentation process (n_fragments) results in (1,max_n_parts) number of fragments (uniformly).
        # It can happen that the entire sequence remains intact.
        breakage_points = [
            select_fragmentation_sites(num_parts=random.randint(1, max_n_parts))
            for _ in range(n_fragments)
        ]
    case "exact":
        # Each fragmentation process (n_fragments) results in exactly max_n_parts number of fragments.
        breakage_points = [
            select_fragmentation_sites(num_parts=max_n_parts) for _ in range(n_fragments)
        ]
    case _:
        raise NotImplementedError(
            f"There is no option for the fragmentation process called '{frag_process}'."
        )


# Generate "left" and "right" nucleotides based on the breakage points,
# which are the exact indices of the nucleotides in the generated fragments:


def generate_left_right(breakage_points):
    left, right = [], []

    for exp in breakage_points:
        left.append(0)

        for part in exp:
            if part != 0:
                right.append(part)
                left.append(part)
        right.append(len(true_sequence))

    return left, right


l_breakage, r_breakage = generate_left_right(breakage_points)

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
match noise_dist:
    case "normal":
        fragment_noise = norm.rvs(scale=rel_error_rate, size=len(true_fragment_masses))
    case "uniform":
        fragment_noise = (
            uniform.rvs(scale=rel_error_rate, size=len(true_fragment_masses))
            - rel_error_rate / 2
        )
    case _:
        raise NotImplementedError(
            f"There is no option for the noise distribution called '{noise_dist}'."
        )

observed_fragment_masses_without_backbone = [
    max(mass + mass * noise, 0.0)
    for noise, mass in zip(fragment_noise, true_fragment_masses)
]

# METHOD: Consider each base in the form of a standard unit, which can be
# combined arbitrarily to build any sequence, and only adapt the masses of the
# fragment ends (either based on a tag or fragmentation/breakage).
def add_backbone_masses(
    fragments, element_masses, fragment_masses, breakage_line="c/y",
    mass_tag_5_prime=mass_5_prime, mass_tag_3_prime=mass_3_prime
):
    # Determine breakage-specific masses for the 5'- and 3'-ends of a fragment
    match breakage_line:
        # TODO: Implement breaking at other points, this is currently c/y-breakage
        case "c/y":
            mass_breakage_5_prime = element_masses["H+"]
            mass_breakage_3_prime = -element_masses["H+"]  # assuming cyclization
            # mass_breakage_3_prime = -element_masses["H+"]  # assuming no cyclization (?)
        case "a/w":
            raise NotImplementedError("Breaking at 'a/w' is not implemented yet.")
        case "b/x":
            raise NotImplementedError("Breaking at 'b/x' is not implemented yet.")
        case "d/z":
            raise NotImplementedError("Breaking at 'd/z' is not implemented yet.")
        case _:
            raise NotImplementedError(
                f"There is no breakage option called '{breakage_line}'."
            )

    for idx, fragment in enumerate(fragments.iter_rows(named=True)):
        # Turn nucleoside masses into those of the corresponding standard units
        fragment_masses[idx] += (fragment["right"]-fragment["left"])*(
            element_masses["P"] + 2 * element_masses["O"] -
            element_masses["H+"]
        )

        # Adapt 5'-end of fragment
        if fragment["left"] == 0:
            # Remove OH and add START tag for terminal fragments
            fragment_masses[idx] += (
                mass_tag_5_prime - element_masses["O"] - element_masses["H+"]
            )
        else:
            # Add breakage-specific mass for internal fragments
            fragment_masses[idx] += mass_breakage_5_prime

        # Adapt 3'-end of fragment
        if fragment["right"] == len(true_sequence):
            # Remove PO3H2 and add END tag for terminal fragments
            fragment_masses[idx] += (
                mass_tag_3_prime - element_masses["P"] -
                3 * element_masses["O"] -
                2 * element_masses["H+"]
            )
        else:
            # Add breakage-specific mass for internal fragments
            fragment_masses[idx] += mass_breakage_3_prime

    return fragment_masses


# Copying the fragment masses to a new list to add the backbone masses to the fragments.
true_fragment_masses_with_backbone = [elem for elem in true_fragment_masses]
true_fragment_masses_with_backbone = add_backbone_masses(
    fragments, element_masses, true_fragment_masses_with_backbone
)


# simulate observed masses with noise -- we use the relative error rate to simulate the noise!
match noise_dist:
    case "normal":
        fragment_noise = norm.rvs(
            scale=rel_error_rate, size=len(true_fragment_masses_with_backbone)
        )
    case "uniform":
        fragment_noise = (
            uniform.rvs(scale=rel_error_rate, size=len(true_fragment_masses_with_backbone))
            - rel_error_rate / 2
        )
    case _:
        raise NotImplementedError(
            f"There is no option for the noise distribution called '{noise_dist}'."
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
    ((pl.col("left") == 0) & (~(pl.col("right") == (len(true_sequence))))).alias(
        "is_start"
    ),
    ((pl.col("right") == (len(true_sequence))) & (~(pl.col("left") == 0))).alias(
        "is_end"
    ),
    ((pl.col("left") == 0) & (pl.col("right") == (len(true_sequence)))).alias(
        "is_start_end"
    ),
    ((~(pl.col("left") == 0)) & (~(pl.col("right") == (len(true_sequence))))).alias(
        "is_internal"
    ),
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
