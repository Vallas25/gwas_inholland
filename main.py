import pandas as pd
from dataframe_operations import *
from plot_operations import *
import numpy as np

def main():
    df = pd.read_csv("data/merged_df.csv")
    result = group_by_location(df)
    create_manhattan_like_plot(df=result[5]["df"], y_axis_name="Pollen Viability", x_axis_name=result[5]["position"])

if __name__ == "__main__":
    main()
