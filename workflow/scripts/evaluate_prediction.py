import polars as pl
import sys
from typing import List

from lionelmssq.common import parse_nucleosides
from lionelmssq.masses import EXPLANATION_MASSES


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
        # results = results.with_columns(pl.lit(0.1).alias("mutation_rate"))
        # num_copies = smk.input[0].split("/")[-1].split(".")[0]
        # results = results.with_columns(pl.lit(num_copies).alias("num_copies"))
        # num_mods = smk.config["fragmentation_params"]["max_singletons"]
        # results = results.with_columns(pl.lit(num_mods).alias("num_singletons"))
        # print(results)
        # results.write_csv(f"results.{num_copies}copies.{num_mods}mods.tsv",
        #                   separator="\t")
        results.write_csv(smk.output[0], separator="\t")


def collect_results(files: List[str]) -> List[str]:
    results = []
    true_sequences = []
    pred_sequences = []
    true_lengths = []
    pred_lengths = []
    for file_path in files:
        true_seq = file_path.split("/")[-2]
        with open(file_path, "r") as f:
            f.readline()
            pred_seq = f.readline().rstrip("\n")
        print(
            len(parse_nucleosides(true_seq)),
            len(parse_nucleosides(pred_seq)),
        )
        print("true:", true_seq)
        print("pred:", pred_seq)
        result = compare_sequences(
            parse_nucleosides(true_seq),
            parse_nucleosides(pred_seq),
        )
        print("result:", result)
        print()
        results.append(result)
        true_sequences.append(true_seq)
        pred_sequences.append(pred_seq)
        true_lengths.append(len(parse_nucleosides(true_seq)))
        pred_lengths.append(len(parse_nucleosides(pred_seq)))

    return pl.DataFrame(
        {
            "true_sequence": true_sequences,
            "pred_sequence": pred_sequences,
            "true_len": true_lengths,
            "pred_len": pred_lengths,
            "result": results,
            "order": [STATUS_ORDER.index(res) for res in results],
        }
    )


def compare_sequences(true_seq: List[str], pred_seq: List[str]) -> str:
    true_seq = [NUC_REPS[nuc] if nuc in NUC_REPS else nuc for nuc in true_seq]
    pred_seq = [NUC_REPS[nuc] if nuc in NUC_REPS else nuc for nuc in pred_seq]
    if len(pred_seq) < 1:
        return "no prediction"
    if len(true_seq) != len(pred_seq):
        return "wrong length"
    if true_seq == pred_seq:
        return "identical"
    if true_seq == ["G" if nuc == "55U" else nuc for nuc in pred_seq]:
        return "identical (minus 55U/G)"
    if sorted(true_seq) == sorted(pred_seq):
        return "correct composition"
    return "failed prediction"


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
