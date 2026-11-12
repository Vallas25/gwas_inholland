import pandas as pd

def change_series_type(df: pd.DataFrame, series_name: str) -> pd.DataFrame:
    # changes the seriestype to numeric
    df.copy()
    df[series_name] = pd.to_numeric(df[series_name], errors="coerce").astype("float64")
    return df

def get_average_per_plant(df: pd.DataFrame, series1_name: str, series2_name: str) -> pd.DataFrame:
    return (
        df.groupby(series1_name)[series2_name].mean().reset_index()
    )

def is_above_mean(df: pd.DataFrame, series_name: str) -> pd.DataFrame:
    df = df.copy()
    mean = df[series_name].mean()

    df[series_name + "_is_high"] = mean < df[series_name]

    return df


if __name__ == "__main__":
    df = pd.read_excel("data/Phenotype_data_PennxMB_population.xlsx")
    df = df.drop(["Unnamed: 10","Unnamed: 11","Unnamed: 12","cellcount/seconds/flowers*120(=seconds to measure complete sample)"],
                  axis=1)
    result = get_average_per_plant(df, "Plant", "#Pollen per flower")
    print(result)