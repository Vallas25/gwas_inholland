import pandas as pd
import scipy.stats as stats

def is_above_mean(df: pd.DataFrame, series_name: str) -> pd.DataFrame:
    # retruns df where is above mean is added
    df = df.copy()
    mean = df[series_name].mean()

    df[series_name + "_is_high"] = mean < df[series_name]

    return df

def t_test_M(df: pd.DataFrame) -> pd.DataFrame:
    # returns dataframe where all non parental plants are compared to the M parental plant and a indipended t-test is performed

    df = df.copy()
    
    # creating dataframe containing only the mesurements of the M plants
    df["Plant"] = df["Plant"].astype(str)
    df_m = df[(df["Plant"].str.contains("M"))].reset_index()

    # creating dataframe not containing plants from either the M or N parrent
    df_other = df[~(df["Plant"].str.contains("M")) & ~(df["Plant"].str.contains("N"))].reset_index()

    df_other["Plant"] = df_other["Plant"].astype(int)

    # grouping by plant, performing a t-test on the pollen viability column
    grouped = df_other.groupby("Plant").agg(
        p_value=("Pollen Viability",
                 lambda plant: stats.ttest_ind(plant, df_m["Pollen Viability"]).pvalue), 
    )

    return grouped.reset_index().sort_values("Plant").reset_index(drop=True)

def check_merge(df_merged: pd.DataFrame):
    #checks merges
    print(df_merged.shape)
    print(df_merged.head())
    print(df_merged["Plant"].duplicated().sum())
    print(df_merged.isna().sum())

def get_sub_loci(df: pd.DataFrame, chr: int = 4, threshold: float | int =4) -> pd.DataFrame:
    # returns df of sublocus position with certain paramerters

    return df[(df["chr_num"] == chr) & (df["-logp"] >= threshold)]

def get_top_performers(df: pd.DataFrame, top_percent: int | float = 10, asc=False) -> pd.DataFrame:
    # gets top performers of plants

    df["Plant"] = df["Plant"].astype(str)
    df = df[~(df["Plant"].str.contains("M")) & ~df["Plant"].str.contains("N")]

    grouped = df.groupby("Plant")["Pollen Viability"].mean()
    grouped = grouped.sort_values(ascending=asc)

    top_n = round(len(grouped) * top_percent / 100)

    return grouped.head(top_n)
