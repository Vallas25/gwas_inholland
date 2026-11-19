import pandas as pd
from dataframe_operations import *
from plot_operations import *
import numpy as np

def main():
    df = pd.read_csv("data/merged_df.csv")
    result = group_by_location(df)
    p_values = get_p_values(result)
    p_values["-logp"] = -np.log10(p_values["p_value"])
    create_manhattan_like_plot(df=p_values, y_axis_name="-logp", x_axis_name="position")

if __name__ == "__main__":
    main()
