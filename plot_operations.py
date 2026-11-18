import plotly.express as px
import pandas as pd

def create_manhattan_like_plot(df: pd.DataFrame, y_axis_name: str, x_axis_name:str="Plant", threshold: float = 7.3):
    fig = px.scatter(df, x=x_axis_name, y=y_axis_name, color=y_axis_name)
    fig.show()