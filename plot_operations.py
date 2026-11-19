import plotly.graph_objects as go
import pandas as pd

def create_manhattan_plot(df: pd.DataFrame, y_axis_name: str, x_axis_name:str="Plant", threshold: float = 7.3):
    fig = go.Figure()

    # Add datapoints
    fig.add_trace(go.Scatter(
        x=df[x_axis_name],
        y=df[y_axis_name],
        mode="markers",
        name="data"
        ))

    # Add threshold line
    fig.add_hline(
        y=threshold,
        line_dash="dash",
        line_color="orange",
        annotation_text=f"significantie niveau",
    )

    fig.update_layout(
        title="Manhattan Plot",
        xaxis_title=x_axis_name,
        yaxis_title=y_axis_name,
    )

    fig.show()