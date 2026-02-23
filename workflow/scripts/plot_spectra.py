import altair as alt
import polars as pl
import sys


STATUS_COLORS = {
    "true": "#990000",
    # "true": "#2b73b5",
    "false": "#808285",
}

LEGEND_PARAMS = {
    "padding": 10,
    "strokeColor": "black",
    "cornerRadius": 5,
    "fillColor": "white",
}

MAX_INTENSITY = 0.2
HISTOGRAM_WIDTH = 850

if "snakemake" in locals():
    smk = snakemake
    sys.stderr = open(smk.log[0], "w")

    def main() -> None:
        pred_frags = pl.read_csv(smk.input["pred_fragments"], separator="\t")
        raw_frags = pl.read_csv(
            smk.input["raw_fragments"], separator="\t"
        ).with_row_index(name="orig_index")

        create_spectra_plot(data_raw=raw_frags, data_pred=pred_frags).save(
            smk.output[0]
        )


def create_spectra_plot(data_raw: pl.DataFrame, data_pred: pl.DataFrame) -> alt.Chart:
    """Create histogram over fragments (raw vs. used in prediction)..

    Parameters
    ----------
    data_raw : polars.Dataframe
        Dataframe containing all raw fragments.
    data_pred : polars.Dataframe
        Dataframe containing all fragments used for prediction.

    Returns
    -------
    altair.Chart
        Histogram over fragments (mass vs. intensity).

    """
    # Select highest measured intensity peak
    max_intensity = data_raw["intensity"].max()

    # Adapt intensity for both dataframes
    data_raw = adjust_intensity_data(data=data_raw, max_intensity=max_intensity)
    data_pred = adjust_intensity_data(data=data_pred, max_intensity=max_intensity)

    # Create histogram for raw fragments
    base_chart = (
        alt.Chart(data_raw)
        .mark_bar(size=1)
        .encode(
            x=alt.X(
                "observed_mass:Q",
                title="Observed mass",
                scale=alt.Scale(padding=0),
                axis=alt.Axis(grid=False),
            ),
            y=alt.Y(
                "plot_intensity:Q",
                title="Relative intensity",
                axis=alt.Axis(grid=False),
                stack=False,
            ),
            tooltip=["orig_index", "observed_mass", "intensity", "rel_intensity"],
            color=alt.value("gray"),
        )
        .properties(width=HISTOGRAM_WIDTH)
    )

    # Create histogram for prediction fragment
    pred_chart = (
        alt.Chart(data_pred)
        .mark_bar(size=1.5)
        .encode(
            x=alt.X(
                "observed_mass:Q",
                title="Observed mass",
                scale=alt.Scale(padding=0),
                axis=alt.Axis(grid=False),
            ),
            y=alt.Y(
                "plot_intensity:Q",
                title="Relative intensity",
                axis=alt.Axis(grid=False),
                stack=False,
            ),
            tooltip=["orig_index", "observed_mass", "intensity", "rel_intensity"],
            color=alt.value("#990000"),
        )
        .properties(width=HISTOGRAM_WIDTH)
    )

    # Add index values for predicted fragments
    text = pred_chart.mark_text(align="center", dy=-5, fontWeight="bold").encode(
        text="orig_index:N", color=alt.value("black")
    )

    # Return combined charts
    return (
        (base_chart + text + pred_chart)
        .configure_view(strokeWidth=0)
        .resolve_scale("shared")
    )


def adjust_intensity_data(data: pl.DataFrame, max_intensity: float):
    """Adapt intensity data to fit histogram.

    Set intensity data to be relative to maximum peak and reformat it to fit
    the maximum displayed relative intensity in the histogram.

    Parameters
    ----------
    data : pl.DataFrame
        Dataframe containing fragments.
    max_intensity : float
        Maximum measured intensity peak.

    Returns
    -------
    pl.DataFrame
        Fragment dataframe with updated intensity data (relative and fit plot)

    """
    data = data.with_columns(
        (pl.col("intensity").cast(float) / max_intensity).alias("rel_intensity")
    )
    return data.with_columns(
        pl.min_horizontal(pl.col("rel_intensity"), MAX_INTENSITY).alias(
            "plot_intensity"
        )
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
