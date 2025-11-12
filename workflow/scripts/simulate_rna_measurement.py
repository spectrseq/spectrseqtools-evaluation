import sys
import os
import polars as pl
import random
import re
import yaml
import numpy as np
from pathlib import Path
from scipy.stats import norm, uniform

# Regex expression to separate given sequence into nucleosides
_NUCLEOSIDE_RE = re.compile(r"\d*[ACGU]")
GHOST_FRAGMENT_MAGNITUDE = 1000
NO_FRAGMENTATION_PROBABILITY = 0.05

if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        """Simulate mass-spectrometry data via Snakemake."""
        # Set random seeds
        seed_path = "/".join(smk.output["meta"].split("/")[:-1] + ["seed.txt"])
        if not os.path.isfile(seed_path):
            seed = 0
        else:
            with open(seed_path, "r") as seed_file:
                seed = int(seed_file.readline().rstrip("\n"))
        random.seed(seed)
        np.random.seed(seed)

        # Build meta dict
        true_sequence = _NUCLEOSIDE_RE.findall(smk.wildcards.seq)
        meta = {
            "identity": "simulated data",
            "label_mass_3T": smk.config["fragmentation_params"]["mass_3_prime"],
            "label_mass_5T": smk.config["fragmentation_params"]["mass_5_prime"],
            "true_sequence": "".join(true_sequence),
        }

        # Build dict with extra masses
        extra_mass_dict = build_extra_mass_dict(
            element_mass_path=smk.input["elements"],
            breakage_line=smk.config["fragmentation_params"]["breakage_line"],
            mass_5_prime=meta["label_mass_5T"],
            mass_3_prime=meta["label_mass_3T"],
        )

        # Add sequence mass to meta dict
        meta["sequence_mass"] = (
            get_seq_weight(
                seq=true_sequence,
                masses=pl.read_csv(smk.input["nucleosides"], separator="\t"),
            )
            + len(true_sequence) * extra_mass_dict["to_standard_unit"]
            + extra_mass_dict["3_prime_terminal"]
            + extra_mass_dict["5_prime_terminal"]
        )

        # Simulate fragments
        simulated_fragments = simulate(
            true_sequence=true_sequence,
            nucleoside_masses=pl.read_csv(smk.input["nucleosides"], separator="\t"),
            n_fragments=int(smk.wildcards.n_fragments),
            ghost_rate=float(smk.config["fragmentation_params"]["ghost_rate"]),
            rel_error_rate=smk.config["fragmentation_params"]["rel_error_rate"],
            noise_dist=smk.config["fragmentation_params"]["noise_distribution"],
            extra_mass_dict=extra_mass_dict,
        )

        # Write simulation data to file
        simulated_fragments.write_csv(smk.output["fragments"], separator="\t")

        with open(smk.output["meta"], "w") as f:
            yaml.safe_dump(meta, f)


# METHOD: Consider each base in the form of a standard unit, which can be
# combined arbitrarily to build any sequence, and only adapt the masses of the
# fragment ends (either based on a tag or fragmentation/breakage).
def build_extra_mass_dict(
    breakage_line: str,
    element_mass_path: Path,
    mass_5_prime: float,
    mass_3_prime: float,
) -> dict:
    # Build dict of elemental masses
    element_masses = pl.read_csv(element_mass_path, separator="\t")
    element_masses = {
        row[element_masses.get_column_index("symbol")]: row[
            element_masses.get_column_index("mass")
        ]
        for row in element_masses.iter_rows()
    }

    # Initialize dict with universal masses
    extra_mass_dict = {
        # Mass needed to turn a nucleoside to a standard unit (SU)
        "to_standard_unit": (
            element_masses["P"] + 2 * element_masses["O"] - element_masses["H+"]
        ),
        # Remove O from SU and add START tag (-H) for 5'-end of terminal fragments
        "5_prime_terminal": mass_5_prime - element_masses["O"] - element_masses["H+"],
        # Remove PO3H from SU and add END tag (-H) for 3'-end of terminal fragments
        "3_prime_terminal": (
            mass_3_prime
            - element_masses["P"]
            - 3 * element_masses["O"]
            - 2 * element_masses["H+"]
        ),
    }

    # Add breakage-specific masses for 5'- and 3'-ends of a fragment to dict
    match breakage_line:
        case "a/w":  # assuming double bond for 3'-end
            extra_mass_dict["5_prime_internal"] = (
                element_masses["P"] + 3 * element_masses["O"] + 2 * element_masses["H+"]
            )
            extra_mass_dict["3_prime_internal"] = -(
                element_masses["P"] + 3 * element_masses["O"] + 2 * element_masses["H+"]
            )
        case "b/x":
            extra_mass_dict["5_prime_internal"] = (
                element_masses["P"] + 2 * element_masses["O"]
            )
            extra_mass_dict["3_prime_internal"] = -(
                element_masses["P"] + 2 * element_masses["O"]
            )
        case "c/y":  # assuming cyclization for 3'-end
            extra_mass_dict["5_prime_internal"] = element_masses["H+"]
            extra_mass_dict["3_prime_internal"] = -element_masses["H+"]
        case "d/z":  # assuming double bond for 5'-end
            extra_mass_dict["5_prime_internal"] = -(
                element_masses["O"] + element_masses["H+"]
            )
            extra_mass_dict["3_prime_internal"] = (
                element_masses["O"] + element_masses["H+"]
            )
        case _:
            raise NotImplementedError(
                f"There is no breakage option called '{breakage_line}'."
            )

    return extra_mass_dict


def get_seq_weight(seq: list, masses: dict) -> float:
    seq_df = pl.DataFrame(data=seq, schema=["name"])
    seq_df = seq_df.with_columns(
        pl.col("name")
        .map_elements(
            lambda x: masses.filter(pl.col("nucleoside") == x)
            .get_column("monoisotopic_mass")
            .to_list()[0],
            return_dtype=pl.Float64,
        )
        .alias("mass")
    )
    return seq_df.select("mass").sum().item()


def simulate(
    true_sequence,
    nucleoside_masses,
    n_fragments,
    ghost_rate,
    rel_error_rate,
    noise_dist,
    extra_mass_dict,
):
    # Sample random fragments from true sequence
    seq_len = len(true_sequence)
    frag_sites = [
        select_fragmentation_sites(select_num_breaks(seq_len=seq_len), seq_len)
        for _ in range(round(n_fragments * (1 + ghost_rate)))
    ]

    # Build fragment dataframe
    fragments = pl.from_records(
        compute_fragment_tuples(frag_sites, len(true_sequence)),
        schema=["left", "right"],
        orient="row",
    )

    # Add columns with boolean values to fragment dataframe
    fragments = fragments.with_columns(
        ((pl.col("left") == 0) & (~(pl.col("right") == len(true_sequence)))).alias(
            "is_start"
        ),
        ((~(pl.col("left") == 0)) & (pl.col("right") == len(true_sequence))).alias(
            "is_end"
        ),
        ((pl.col("left") == 0) & (pl.col("right") == len(true_sequence))).alias(
            "is_start_end"
        ),
        ((~(pl.col("left") == 0)) & (~(pl.col("right") == len(true_sequence)))).alias(
            "is_internal"
        ),
        (pl.col("right") == (pl.col("left") + 1)).alias("single_nucleoside"),
    )

    # Add column with true sequences to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("left", "right")
        .map_elements(
            lambda x: "".join(true_sequence[x["left"] : x["right"]]), return_dtype=str
        )
        .alias("sequence")
    )

    # Add column with exact nucleoside masses to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("left", "right")
        .map_elements(
            lambda x: sum(
                nucleoside_masses.filter(pl.col("nucleoside") == base)
                .select(pl.col("monoisotopic_mass"))
                .item()
                for base in true_sequence[x["left"] : x["right"]]
            ),
            return_dtype=float,
        )
        .alias("true_nucleoside_mass")
    )

    # Add column with noisy nucleoside masses to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("*")
        .map_elements(
            lambda x: induce_noise(
                noise_dist, rel_error_rate, x["true_nucleoside_mass"]
            ),
            return_dtype=float,
        )
        .alias("observed_mass_without_backbone")
    )

    # Add column with exact nucleotide masses to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("*")
        .map_elements(
            lambda x: add_backbone_mass(
                fragment=x,
                mass=x["true_nucleoside_mass"],
                seq_len=len(true_sequence),
                mass_dict=extra_mass_dict,
            ),
            return_dtype=float,
        )
        .alias("true_mass_with_backbone")
    )

    # Add column with noisy nucleotide masses to fragment dataframe
    fragments = fragments.with_columns(
        pl.struct("*")
        .map_elements(
            lambda x: induce_noise(
                noise_dist, rel_error_rate, x["true_mass_with_backbone"]
            ),
            return_dtype=float,
        )
        .alias("observed_mass")
    )

    # Select ghost (i.e. invalid) fragments
    fragments = fragments.with_columns(
        pl.struct("*")
        .map_elements(
            lambda x: True if uniform.rvs() < ghost_rate else False,
            return_dtype=bool,
        )
        .alias("is_ghost_fragment")
    )

    # Update classification for ghost fragments by setting all to internal
    fragments = fragments.with_columns(
        (pl.col("is_start") & ~pl.col("is_ghost_fragment")).alias("is_start"),
        (pl.col("is_end") & ~pl.col("is_ghost_fragment")).alias("is_end"),
        (pl.col("is_start_end") & ~pl.col("is_ghost_fragment")).alias("is_start_end"),
        (pl.col("is_internal") | pl.col("is_ghost_fragment")).alias("is_internal"),
    )

    # Update observed mass for ghost fragments by adjusting it randomly
    fragments = fragments.with_columns(
        pl.struct("*")
        .map_elements(
            lambda x: x["observed_mass"]
            + int(x["is_ghost_fragment"])
            * uniform.rvs(
                loc=-GHOST_FRAGMENT_MAGNITUDE, scale=2 * GHOST_FRAGMENT_MAGNITUDE
            ),
            return_dtype=float,
        )
        .alias("observed_mass")
    )

    # Return final fragment dataframe
    return fragments


# METHOD: Consider fragments without any breakage, i.e. complete fragments,
# separately (randomly select based on given probability); if the sequence
# does break, use a geometric distribution to determine the number of breaks
# while approximating the true distribution of fragment lengths observed in
# experimental data (exponential distribution with many small and few larger
# fragments, which gets sharper with increasing sequence length)
def select_num_breaks(seq_len: int) -> int:
    if np.random.rand() < NO_FRAGMENTATION_PROBABILITY:
        return 0
    # Note that p = factor/seq_len with factor = seq_len/alpha
    # thus using p = seq_len/alpha/seq_len = 1/alpha
    return min(np.random.geometric(p=0.3), seq_len - 1)


# TODO: Implement that in some cases there is no base pair generated, but only the backbone with sugar etc?
def select_fragmentation_sites(num_breaks, seq_len):
    # Ensure there is a positive number of parts (i.e. number of breaks + 1)
    if num_breaks < 0:
        raise ValueError("The number of parts cannot be less than one!")

    # Ensure the number of parts is not greater than the sequence length
    if num_breaks + 1 > seq_len:
        raise ValueError(
            "The number of parts cannot be greater than the sequence length!"
        )

    # If the sequence has zero breaks, it remains intact
    if num_breaks == 0:
        return [int(0)]

    # Return randomly sampled breakage positions in the sequence
    return sorted(set(np.random.choice(range(1, seq_len), num_breaks)))


def compute_fragment_tuples(frag_sites, seq_len):
    tuples = []

    # Generate tuples of start and end index for each fragments
    for seq_copy in frag_sites:
        if seq_copy[0] != 0:
            seq_copy.insert(0, 0)
        if seq_copy[-1] != seq_len:
            seq_copy.append(seq_len)

        tuples += list(zip(seq_copy[:-1], seq_copy[1:]))

    return tuples


def induce_noise(distribution_method, error_rate, mass):
    match distribution_method:
        case "normal":
            noise = norm.rvs(scale=error_rate)
        case "uniform":
            noise = uniform.rvs(loc=-error_rate, scale=2 * error_rate)
        case _:
            raise NotImplementedError(
                f"There is no option for the noise distribution called '{distribution_method}'."
            )

    return max(mass * (1 + noise), 0.0)


def add_backbone_mass(fragment, mass, seq_len, mass_dict):
    # Turn nucleoside mass into the one of the corresponding standard units
    mass += (fragment["right"] - fragment["left"]) * mass_dict["to_standard_unit"]

    # Adapt 5'-end of fragment
    mass += (
        mass_dict["5_prime_terminal"]
        if fragment["left"] == 0
        else mass_dict["5_prime_internal"]
    )

    # Adapt 3'-end of fragment
    mass += (
        mass_dict["3_prime_terminal"]
        if fragment["right"] == seq_len
        else mass_dict["3_prime_internal"]
    )

    return mass


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
