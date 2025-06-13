import sys
import polars as pl
import random
import re
import numpy as np
from scipy.stats import norm, uniform

# regex for separating given sequence into nucleosides
nucleoside_re = re.compile(r"\d*[ACGU]")

if 'snakemake' in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        """Simulate mass-spectrometry data via Snakemake."""
        random.seed(smk.config["fragmentation_params"]["random_seed"])
        np.random.seed(0)

        # Input
        true_sequence = nucleoside_re.findall(smk.wildcards.seq)
        n_fragments = int(smk.wildcards.n_fragments)
        nucleoside_masses = pl.read_csv(smk.input["nucleosides"], separator="\t")
        element_masses = pl.read_csv(smk.input["elements"], separator="\t")

        max_n_parts = int(smk.config["fragmentation_params"]["max_n_parts"])
        rel_error_rate = smk.config["fragmentation_params"]["rel_error_rate"]
        breakage_line = smk.config["fragmentation_params"]["breakage_line"]
        frag_process = smk.config["fragmentation_params"]["fragmentation_process"]
        noise_dist = smk.config["fragmentation_params"]["noise_distribution"]
        mass_5_prime = smk.config["fragmentation_params"]["mass_5_prime"]
        mass_3_prime = smk.config["fragmentation_params"]["mass_3_prime"]

        # Convert dataframe of element weights to dict
        element_masses = {row[element_masses.get_column_index("symbol")]:
                          row[element_masses.get_column_index("mass")] for
                          row in element_masses.iter_rows()}

        # Simulate fragments
        simulated_fragments = simulate(
            frag_process, max_n_parts, true_sequence, nucleoside_masses,
            element_masses, n_fragments, rel_error_rate, breakage_line,
            noise_dist, mass_5_prime, mass_3_prime
        )

        # Write simulation data to file
        simulated_fragments.write_csv(smk.output[0], separator="\t")

# helper functions
# def random_fragment():
#    l = random.randint(0, len(true_sequence) - 1)
#    r = random.randint(l + 1, len(true_sequence))
#    return l, r

# TODO: Implement that in some cases there is no base pair generated, but only the backbone with sugar etc?

# Divide the sequence into a given number of "max_n_parts".
def simulate(frag_process, max_n_parts, true_sequence, nucleoside_masses,
             element_masses, n_fragments, rel_error_rate, breakage_line,
             noise_dist, mass_5_prime, mass_3_prime):
    # Sample random fragments from true sequence
    frag_sites = [
        select_fragmentation_sites(
            select_num_of_breaks(frag_process, max_n_parts), len(true_sequence))
        for _ in range(n_fragments)
    ]

    # Build fragment dataframe
    fragments = pl.from_records(
        compute_fragment_tuples(frag_sites, len(true_sequence)),
        schema=["left", "right"],
        orient="row",
    )

    # Add columns with boolean values to fragment dataframe
    fragments = fragments.with_columns(
        ((pl.col("left") == 0) & (~(pl.col("right") == len(true_sequence)))).alias("is_start"),
        ((~(pl.col("left") == 0)) & (pl.col("right") == len(true_sequence))).alias("is_end"),
        ((pl.col("left") == 0) & (pl.col("right") == len(true_sequence))).alias("is_start_end"),
        ((~(pl.col("left") == 0)) & (~(pl.col("right") == len(true_sequence)))).alias("is_internal"),
        (pl.col("right") == (pl.col("left") + 1)).alias("single_nucleoside"),
    )

    # Add column with true sequences to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("left", "right").map_elements(
            lambda x: "".join(true_sequence[x["left"]:x["right"]]),
            return_dtype=str
        ).alias("sequence")
    )

    # Add column with exact nucleoside masses to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("left", "right").map_elements(
            lambda x: sum(
                nucleoside_masses.filter(pl.col("nucleoside") == base)
                .select(pl.col("monoisotopic_mass"))
                .item()
                for base in
                true_sequence[x["left"]:x["right"]]
            ),
            return_dtype=float
        ).alias("true_nucleoside_mass")
    )

    # Add column with noisy nucleoside masses to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("*").map_elements(
            lambda x: induce_noise(noise_dist, rel_error_rate, x["true_nucleoside_mass"]),
            return_dtype=float
        ).alias("observed_mass_without_backbone")
    )

    # Add column with exact nucleotide masses to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("*").map_elements(
            lambda x: add_backbone_mass(
                x, element_masses, x["true_nucleoside_mass"],
                len(true_sequence), mass_5_prime, mass_3_prime, breakage_line
            ),
            return_dtype=float
        ).alias("true_mass_with_backbone")
    )

    # Add column with noisy nucleotide masses to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("*").map_elements(
            lambda x: induce_noise(noise_dist, rel_error_rate,
                                   x["true_mass_with_backbone"]),
            return_dtype=float
        ).alias("observed_mass")
    )

    # Return final fragment dataframe
    return fragments

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

def select_fragmentation_sites(num_parts, seq_len):
    # Ensure there is a positive number of parts
    if num_parts <= 0:
        raise ValueError(
            "The number of parts cannot be less than one!"
        )

    # Ensure the number of parts is not greater than the sequence length
    if num_parts > seq_len:
        raise ValueError(
            "The number of parts cannot be greater than the sequence length!"
        )

    # If the sequence "breaks" into one part, it remains intact
    if num_parts == 1:
        return [int(0)]

    # Return randomly sampled breakage positions in the sequence
    return sorted(random.sample(range(1, seq_len), num_parts - 1))
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

def induce_noise(distribution_method, error_rate, mass):
    match distribution_method:
        case "normal":
            noise = norm.rvs(scale=error_rate)
        case "uniform":
            noise = uniform.rvs(scale=error_rate) - error_rate / 2
        case _:
            raise NotImplementedError(
                f"There is no option for the noise distribution called '{distribution_method}'."
            )

    return max(mass * (1 + noise), 0.0)

# METHOD: Consider each base in the form of a standard unit, which can be
# combined arbitrarily to build any sequence, and only adapt the masses of the
# fragment ends (either based on a tag or fragmentation/breakage).
def add_backbone_mass(
        fragment, element_masses, mass, seq_len, mass_tag_5_prime,
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

    # Turn nucleoside mass into the one of the corresponding standard units
    mass += (fragment["right"]-fragment["left"]) * (
            element_masses["P"]+2 * element_masses["O"]-element_masses["H+"]
    )

    # Adapt 5'-end of fragment
    if fragment["left"] == 0:
        # Remove OH and add START tag for terminal fragments
        mass += mass_tag_5_prime-element_masses["O"]-element_masses["H+"]
    else:
        # Add breakage-specific mass for internal fragments
        mass += mass_breakage_5_prime

    # Adapt 3'-end of fragment
    if fragment["right"] == seq_len:
        # Remove PO3H2 and add END tag for terminal fragments
        mass += (
                mass_tag_3_prime-element_masses["P"]-
                3 * element_masses["O"]-
                2 * element_masses["H+"]
        )
    else:
        # Add breakage-specific mass for internal fragments
        mass += mass_breakage_3_prime

    return mass

if __name__ == '__main__':
    if "snakemake" in locals():
        main()
    else:
        print('Not Defined.')