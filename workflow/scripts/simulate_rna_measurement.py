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
def simulate(frag_process, max_n_parts, true_sequence, nucleoside_masses,
             element_masses, n_fragments, rel_error_rate, breakage_line,
             noise_dist, mass_5_prime, mass_3_prime):
    # TODO: Implement that in some cases there is no base pair generated, but only the backbone with sugar etc?

    # Sample random fragments from true sequence
    frag_sites = [
        select_fragmentation_sites(
            select_num_of_breaks(frag_process, max_n_parts))
        for _ in range(n_fragments)
    ]
    fragments = pl.from_records(
        compute_fragment_tuples(frag_sites, len(true_sequence)),
        schema=["left", "right"],
        orient="row",
    )

    # Calculate true masses based on known nucleoside masses
    true_fragment_masses = [
        sum(
            nucleoside_masses.filter(pl.col("nucleoside") == b)
            .select(pl.col("monoisotopic_mass"))
            .item()
            for b in true_sequence[fragment["left"]: fragment["right"]]
        )
        for fragment in fragments.iter_rows(named=True)
    ]

    # Calculate observed masses WITHOUT the backbone masses
    observed_fragment_masses_without_backbone = induce_noise(
        noise_dist, rel_error_rate, true_fragment_masses
    )

    # Copying the fragment masses to a new list to add the backbone masses to the fragments
    true_fragment_masses_with_backbone = add_backbone_masses(
        fragments, element_masses, true_fragment_masses.copy(),
        mass_5_prime, mass_3_prime
    )

    # Simulate observed masses with noise using the relative error rate
    observed_fragment_masses = induce_noise(
        noise_dist, rel_error_rate, true_fragment_masses_with_backbone
    )

    # Get the fragment sequences
    fragment_sequences = [
        "".join(true_sequence[fragment["left"] : fragment["right"]])
        for fragment in fragments.iter_rows(named=True)
    ]

    # Compile final dataframe
    return fragments.with_columns(
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

def select_num_of_breaks(frag_process, max_n_parts):
    # Select how many parts the sequence gets broken into.
    match frag_process:
        case "random":
            # At high energies the sequence can be divided into multiple parts! Sample uniformly how many parts it gets broken into.
            # Each fragmentation process (n_fragments) results in (1,max_n_parts) number of fragments (uniformly).
            # It can happen that the entire sequence remains intact.
            return random.randint(1, max_n_parts)
        case "exact":
            # Each fragmentation process (n_fragments) results in exactly max_n_parts number of fragments.
            return max_n_parts
        case _:
            raise NotImplementedError(
                f"There is no option for the fragmentation process called '{frag_process}'."
            )

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

# Generate tuples of start and end index for each fragments
def compute_fragment_tuples(frag_sites, seq_len):
    tuples = []

    for seq_copy in frag_sites:
        if seq_copy[0] != 0:
            seq_copy.insert(0,0)
        if seq_copy[-1] != seq_len:
            seq_copy.append(seq_len)

        tuples += list(zip(seq_copy[:-1], seq_copy[1:]))

    return tuples

def induce_noise(distribution_method, error_rate, mass_list):
    match distribution_method:
        case "normal":
            noise_list = norm.rvs(scale=error_rate, size=len(mass_list))
        case "uniform":
            noise_list = (uniform.rvs(scale=error_rate, size=len(mass_list)) -
                      error_rate / 2)
        case _:
            raise NotImplementedError(
                f"There is no option for the noise distribution called '{noise_dist}'."
            )

    return [max(mass * (1 + noise), 0.0) for mass, noise in zip(mass_list, noise_list)]

# METHOD: Consider each base in the form of a standard unit, which can be
# combined arbitrarily to build any sequence, and only adapt the masses of the
# fragment ends (either based on a tag or fragmentation/breakage).
def add_backbone_masses(
    fragments, element_masses, fragment_masses, mass_tag_5_prime,
    mass_tag_3_prime, breakage_line="c/y"
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

# Write simulation output
simulated_fragments = simulate(
    frag_process, max_n_parts, true_sequence, nucleoside_masses,
    element_masses, n_fragments, rel_error_rate, breakage_line, noise_dist,
    mass_5_prime, mass_3_prime)
simulated_fragments.write_csv(snakemake.output[0], separator="\t")
