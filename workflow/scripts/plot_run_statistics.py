import altair as alt
import polars as pl
import sys


STATUS_COLORS = {
    "identical": "#8cb110",
    "identical (minus 55U/G)": "yellow",
    "correct composition": "orange",
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

WIDTH = 850

if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        mode = smk.params["mode"]

        sim_results = pl.read_csv(smk.input["sim"], separator="\t")
        sim_chart = create_scatterplot(data=sim_results, color="#990000", mode=mode)

        exp_results = pl.read_csv(smk.input["exp"], separator="\t")
        exp_chart = create_scatterplot(
            data=exp_results, color="purple", size=25, mode=mode
        )

        chart = sim_chart + exp_chart

        chart.save(smk.output[0])


def create_scatterplot(data: pl.DataFrame, color: str, mode: str, size: int = 10):
    base_chart = (
        alt.Chart(data)
        .encode(
            x=alt.X(
                "num_frag:Q", title="Number of fragments", scale=alt.Scale(padding=0)
            ),
            y=select_y_axis(mode=mode),
        )
        .properties(width=WIDTH)
    )

    scatterplot = base_chart.mark_circle(size=size).encode(
        color=alt.value(color),
        tooltip=["s", "num_frag"],
    )

    return scatterplot


def select_y_axis(mode: str):
    match mode:
        case "runtime":
            return alt.Y(
                "s:Q",
                title="Runtime (in sec)",
                scale=alt.Scale(type="linear"),
            )
        case "memory":
            return alt.Y(
                "max_rss:Q",
                title="Memory",
                scale=alt.Scale(type="linear"),
            )
        case _:
            return alt.Y(
                "s:Q",
                title="Runtime (in sec)",
                scale=alt.Scale(type="linear"),
            )


def select_tooltip(mode: str):
    match mode:
        case "runtime":
            return ["s", "num_frag"]
        case "memory":
            return ["max_rss", "num_frag"]
        case _:
            return ["s", "num_frag"]


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
