import re
import pandas as pd
import altair as alt
import sys
from typing import List

from lionelmssq.masses import EXPLANATION_MASSES


_NUCLEOSIDE_RE = re.compile(r"(\d*[ACGU])")
_REP_IDX = EXPLANATION_MASSES.get_column_index("nucleoside")
_LIST_IDX = EXPLANATION_MASSES.get_column_index("nucleoside_list")
NUC_REPS = {
    **{
        nuc: row[_REP_IDX]
        for row in EXPLANATION_MASSES.rows()
        for nuc in row[_LIST_IDX]
    }
}

STATUS_COLORS = {
    "identical": "#8cb110",
    "identical (minus 55U/G)": "yellow",
    "correct composition": "purple",
    # "failed prediction": "#be0a26",
    "failed prediction": "#990000",
    "wrong length": "gray",
}


if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        results = collect_results(smk.input)
        donut_chart = plot_results(results)
        donut_chart.save(smk.output[0])


def plot_results(results):
    results = pd.DataFrame(results)
    return (
        alt.Chart(results)
        .mark_arc(innerRadius=32, outerRadius=50)
        .encode(
            theta=alt.Theta("count(result):Q"),
            color=alt.Color(
                "result:N",
                scale=alt.Scale(
                    domain=STATUS_COLORS.keys(),
                    range=[
                        STATUS_COLORS[stat] if STATUS_COLORS.get(stat) else stat
                        for stat in STATUS_COLORS.keys()
                    ],
                ),
            ),
            tooltip=["result", "count(result)"],
        )
    )


def collect_results(files: List[str]) -> List[str]:
    results = []
    sequences = []
    for file_path in files:
        true_seq = file_path.split("/")[-2]
        with open(file_path, "r") as f:
            f.readline()
            pred_seq = f.readline().rstrip("\n")
        print(
            len(re.findall(_NUCLEOSIDE_RE, true_seq)),
            len(re.findall(_NUCLEOSIDE_RE, pred_seq)),
        )
        print("true:", true_seq)
        print("pred:", pred_seq)
        result = compare_sequences(
            re.findall(_NUCLEOSIDE_RE, true_seq),
            re.findall(_NUCLEOSIDE_RE, pred_seq),
        )
        print("result:", result)
        print()
        results.append(result)
        sequences.append(true_seq)
    return pd.DataFrame({"sequence": sequences, "result": results})


def compare_sequences(true_seq: List[str], pred_seq: List[str]) -> str:
    true_seq = [NUC_REPS[nuc] if nuc in NUC_REPS else nuc for nuc in true_seq]
    pred_seq = [NUC_REPS[nuc] if nuc in NUC_REPS else nuc for nuc in pred_seq]
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
