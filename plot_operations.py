import plotly.graph_objects as go
import pandas as pd
import numpy as np

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
    fig.write_html("plots/manhattan_plot.html")

def create_genomic_heatmap(df: pd.DataFrame,
                           metadata: list = ["Markers", "Chromosome", "location", "chr_num"]):
    df = df.sort_values(["chr_num", "location"]).reset_index(drop=True)
    genotype = df.drop(columns=metadata).apply(pd.to_numeric, errors="coerce")
    n_rows, n_cols = genotype.shape
    z = genotype.values

    # Hover text: a plain 2D list of strings (robust to serialize)
    label_map = {0: "Parent 1", 1: "Parent 2", 2: "Heterozygous"}
    labels = np.full(z.shape, "Missing", dtype=object)
    for code, name in label_map.items():
        labels[z == code] = name
    markers = df["Markers"].astype(str).tolist()
    text = [[f"Marker: {m}<br>Genotype: {g}" for g in row]
            for m, row in zip(markers, labels)]

    fig = go.Figure(data=go.Heatmap(
        z=z,
        x=genotype.columns.astype(str).tolist(),
        y=list(range(n_rows)),
        text=text,
        colorscale=[
            [0.0, "#2166ac"], [1/3, "#2166ac"],
            [1/3, "#bdbdbd"], [2/3, "#bdbdbd"],
            [2/3, "#b2182b"], [1.0, "#b2182b"],
        ],
        zmin=-0.5,
        zmax=2.5,
        colorbar=dict(
            title="Genotype",
            tickvals=[0, 1, 2],
            ticktext=list(label_map.values()),
        ),
        hovertemplate="Plant: %{x}<br>%{text}<extra></extra>",
    ))

    # Chromosome blocks
    blocks = (
        df.reset_index()
          .groupby("chr_num", sort=True)
          .agg(start=("index", "min"), end=("index", "max"), name=("Chromosome", "first"))
          .reset_index()
    )

    strip_colors = ["#4d4d4d", "#a6a6a6"]
    tickvals, ticktext = [], []

    for i, row in enumerate(blocks.itertuples(index=False)):
        top, bottom = row.start - 0.5, row.end + 0.5

        fig.add_shape(
            type="rect", xref="paper", yref="y",
            x0=-0.012, x1=-0.002, y0=top, y1=bottom,
            fillcolor=strip_colors[i % 2], line_width=0,
        )
        if i > 0:
            fig.add_shape(
                type="line", xref="paper", yref="y",
                x0=0, x1=1, y0=top, y1=top,
                line=dict(color="black", width=1.5),
            )

        tickvals.append((row.start + row.end) / 2)
        ticktext.append(str(row.name).replace("SL4.0ch", ""))

    fig.update_layout(
        title="SNP Genotype Heatmap",
        xaxis_title="Plants",
        yaxis_title="Chromosome",
        template="plotly_white",
        height=750,
        margin=dict(l=110, b=110),
    )

    fig.update_yaxes(
        range=[n_rows - 0.5, -0.5],
        tickmode="array",
        tickvals=tickvals,
        ticktext=ticktext,
        showgrid=False,
        ticks="outside",
        ticklen=22,                      # pushes labels left, past the strip
        tickcolor="rgba(0,0,0,0)",       # makes the tick marks themselves invisible
    )

    # Show roughly 60 plant labels at most, so they don't overlap
    step = max(1, n_cols // 60)
    fig.update_xaxes(
        type="category",
        showgrid=False,
        showticklabels=True,
        tickangle=-90,
        tickfont=dict(size=9),
        dtick=step,
    )

    return fig