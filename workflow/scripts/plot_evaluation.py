import altair as alt
import polars as pl
import sys


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
        results = pl.read_csv(smk.input[0], separator="\t")

        donut_chart = create_donut_plot(results)
        donut_chart.save(smk.output["donut"])

        bar_chart = create_stacked_barplot(results)
        bar_chart.save(smk.output["bar"])


def create_donut_plot(results):
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
            x=alt.X("true_len:N", title="Sequence Length"),
            # x=alt.X("num_singletons:N", title="Maximum Number of Singletons"),
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


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
