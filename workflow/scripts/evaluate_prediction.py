import polars as pl
import sys
from typing import List

from spectrseqtools.dataclasses import Sequence
from spectrseqtools.error_calculator import ErrorUnderL1Norm
from spectrseqtools.nucleotide_alphabet import NucleotideAlphabet

NUCLEOTIDE_DF = NucleotideAlphabet.from_file(error=ErrorUnderL1Norm()).to_dataframe()
NUC_REPS = {
    **{
        nuc: row[NUCLEOTIDE_DF.get_column_index("names")][0]
        for row in NUCLEOTIDE_DF.rows()
        for nuc in row[NUCLEOTIDE_DF.get_column_index("names")]
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
    match smk.params["mode"]:
        case "optimization":
            return smk.config["optimization"][smk.wildcards.parameter]
        case "simulation":
            return smk.config["comparison"]["studies"][smk.wildcards.parameter]
        case _:
            return {}


def collect_results(files: List[str]) -> List[str]:
    comp_values = []
    results = []
    true_sequences = []
    pred_sequences = []
    true_lengths = []
    pred_lengths = []
    for file_path in files:
        comp_value = file_path.split("/")[-3]
        true_seq = file_path.split("/")[-2]
        with open(file_path, "r") as f:
            f.readline()
            pred_seq = f.readline().rstrip("\n")
        print(
            len(Sequence.from_str(true_seq).sequence),
            len(Sequence.from_str(pred_seq).sequence),
        )
        print("true:", true_seq)
        print("pred:", pred_seq)
        result = compare_sequences(
            Sequence.from_str(true_seq).sequence,
            Sequence.from_str(pred_seq).sequence,
        )
        print("result:", result)
        print()
        results.append(result)
        true_sequences.append(true_seq)
        pred_sequences.append(pred_seq)
        true_lengths.append(len(Sequence.from_str(true_seq).sequence))
        pred_lengths.append(len(Sequence.from_str(pred_seq).sequence))
        comp_values.append(comp_value)

    return pl.DataFrame(
        {
            "true_sequence": true_sequences,
            "pred_sequence": pred_sequences,
            "true_len": true_lengths,
            "pred_len": pred_lengths,
            "comp_val": comp_values,
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
