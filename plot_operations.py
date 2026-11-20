import plotly.graph_objects as go
import pandas as pd


def create_manhattan_plot(df: pd.DataFrame, x_axis_name:str="Plant",  y_axis_name: str="-logp", threshold: float = 7.3):
    fig = go.Figure()

    colors = [
    "red",
    "violet",
    "blue",
    "green",
    "yellow",
    "orange",]

    for i, chrom in enumerate(sorted(df["chr_num"].unique())):
        subset = df[df["chr_num"] == chrom]

        fig.add_trace(go.Scattergl(
            x=subset[x_axis_name],
            y=subset[y_axis_name],
            mode="markers",
            name=f"Chr {chrom}",
            marker=dict(
                color=colors[i % len(colors)],
                size=4
            )
        ))

    # Add threshold line
    fig.add_hline(
        y=threshold,
        line_dash="dash",
        line_color="orange",
        annotation_text=f"significantie niveau: {threshold}",
    )

    fig.update_layout(
        title="Manhattan Plot",
        xaxis_title=x_axis_name,
        yaxis_title=y_axis_name,
    )

    fig.show()