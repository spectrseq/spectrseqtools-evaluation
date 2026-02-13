import polars as pl
import sys
from typing import List

from spectrseqtools.masses import EXPLANATION_MASSES


NUC_REPS = {
    **{
        nuc: row[EXPLANATION_MASSES.get_column_index("nucleoside")]
        for row in EXPLANATION_MASSES.rows()
        for nuc in row[EXPLANATION_MASSES.get_column_index("nucleoside_list")]
    }
}

STATUS_ORDER = [
    "identical",
    "identical (minus 55U/G)",
    "correct composition",
    "failed prediction",
    "wrong length",
    "no prediction",
]


if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        results = collect_results(smk.input)

        if smk.wildcards:
            params = select_params(smk)
            for key in params:
                if key == smk.wildcards.parameter:
                    continue
                results = results.with_columns(pl.lit(params[key][0]).alias(key))
            results = results.rename({"comp_val": smk.wildcards.parameter})

        results.write_csv(smk.output[0], separator="\t")


def select_params(smk) -> dict:
    # print(type(smk))
    match smk.params["mode"]:
        case "optimization":
            return smk.config["optimization"][smk.wildcards.parameter]
        case "simulation":
            return smk.config["comparison"]["studies"][smk.wildcards.parameter]
        case _:
            return {}


def collect_results(input_dict: dict) -> List[str]:
    results = []
    for benchmark, fragments in zip(input_dict["benchmarks"], input_dict["fragments"]):
        new_data = pl.read_csv(benchmark, separator="\t")
        # frag = pl.read_csv(fragments, separator="\t")
        # if "ppm_group" in frag.columns:
        #      num_entries = len(frag.select("ppm_group").unique())
        # else:
        with open(fragments, "r") as f:
            num_entries = len(f.readlines())
        new_data.insert_column(0, pl.Series("num_frag", [num_entries]))
        results.append(new_data)

    return pl.concat(results)


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
