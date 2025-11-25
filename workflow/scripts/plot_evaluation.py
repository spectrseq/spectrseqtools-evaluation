import altair as alt
import polars as pl
import sys
from typing import List

from lionelmssq.common import parse_nucleosides
from lionelmssq.masses import EXPLANATION_MASSES


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
    "no prediction": "black",
}

STATUS_ORDER = [
    "identical",
    "identical (minus 55U/G)",
    "correct composition",
    "failed prediction",
    "wrong length",
    "no prediction",
]

LEGEND_PARAMS = {
    "padding": 10,
    "strokeColor": "black",
    "cornerRadius": 5,
    "fillColor": "white",
}


if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        results = collect_results(smk.input)
        donut_chart = plot_results(results)
        donut_chart.save(smk.output["donut"])

        bar_chart = create_stacked_barplot(results)
        bar_chart.save(smk.output["bar"])


def plot_results(results):
    results = pl.DataFrame(results)
    return (
        alt.Chart(results)
        .mark_arc(innerRadius=32, outerRadius=50)
        .encode(
            theta=alt.Theta("count(result):Q", sort=STATUS_ORDER),
            color=alt.Color(
                "result:N",
                scale=alt.Scale(
                    domain=STATUS_ORDER,
                    range=[
                        STATUS_COLORS[stat] if STATUS_COLORS.get(stat) else stat
                        for stat in STATUS_ORDER
                    ],
                ),
            ),
            order=alt.Order("order:O", sort="descending"),
            tooltip=["result", "count(result)"],
        )
    )


def create_stacked_barplot(data: pl.DataFrame) -> alt.Chart:
    """Create stacked barplot over prediction status.

    Parameters
    ----------
    data : polars.Dataframe
        Dataframe containing prediction status data.

    Returns
    -------
    altair.Chart
        Stacked barplot over prediction status.

    """
    chart = (
        alt.Chart(data, title="Status Assessment")
        .mark_bar()
        .encode(
            x=alt.X("num_singletons:N", title="Maximum Number of Singletons"),
            y=alt.Y("count(result):Q", sort=STATUS_ORDER, title="Number of Sequences"),
            color=alt.Color(
                "result:N",
                scale=alt.Scale(
                    domain=STATUS_ORDER,
                    range=[
                        STATUS_COLORS[stat] if STATUS_COLORS.get(stat) else stat
                        for stat in STATUS_ORDER
                    ],
                ),
                legend=alt.Legend(
                    **LEGEND_PARAMS, orient="right", title="Prediction Status"
                ),
                sort=STATUS_ORDER,
            ),
            order=alt.Order("order:O", sort="descending"),
            tooltip=["result", "count(result)"],
        )
    )
    return chart


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
