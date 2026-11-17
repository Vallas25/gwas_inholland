import plotly.express as px
import pandas as pd

def create_manhattan_like_plot(df: pd.DataFrame, y_axis_name: str, threshold: float = 7.3):
    fig = px.scatter(df, x="Plant", y=y_axis_name, color=y_axis_name)
    fig.show()