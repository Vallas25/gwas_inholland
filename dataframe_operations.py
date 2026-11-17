import pandas as pd
import scipy.stats as stats

def get_average_per_plant(df: pd.DataFrame, series1_name: str, series2_name: str) -> pd.DataFrame:
    # groups series, then gets average by grouped.
    return (
        df.groupby(series1_name)[series2_name].mean().reset_index()
    )

def clean_phenotype_data(df: pd.DataFrame) -> pd.DataFrame:
    # removes unnecacarry rows from phenotype_data dataset
    return df.drop(["Unnamed: 10","Unnamed: 11","Unnamed: 12","cellcount/seconds/flowers*120(=seconds to measure complete sample)"],
                     axis=1)

def is_above_mean(df: pd.DataFrame, series_name: str) -> pd.DataFrame:
    # retruns df where is above mean is added
    df = df.copy()
    mean = df[series_name].mean()

    df[series_name + "_is_high"] = mean < df[series_name]

    return df

def filter_failed_hetrozygote(df: pd.DataFrame) -> pd.DataFrame:
    # returns df where all failed mesurements and hetrozygote loci are removed
    df = df.copy()

    allele_columns = ["Parent1.1", "Parent1.2", "Parent2.1", "Parent2.2"]

    mask_heterozygote = (
        df[allele_columns].eq(df["Parent1.1"], axis=0).all(axis=1)
    )

    mask_failed = (
        df.eq("failed").any(axis=1)
    )

    df = df.loc[~mask_heterozygote]
    return df.loc[~mask_failed]

def t_test_M(df: pd.DataFrame) -> pd.DataFrame:
    # returns dataframe where all non parental plants are compared to the M parental plant and a indipended t-test is performed

    df = df.copy()
    
    df["Plant"] = df["Plant"].astype(str)
    df_m = df[(df["Plant"].str.contains("M"))]

    df_other = df[~(df["Plant"].str.contains("M")) & ~(df["Plant"].str.contains("N"))]

    grouped = df_other.groupby("Plant").agg(
        p_value=("Pollen Viability", lambda plant: stats.ttest_ind(plant, df_m["Pollen Viability"]).pvalue), 
    )
    return grouped.reset_index()


if __name__ == "__main__":
    df = pd.read_excel("data/Phenotype_data_PennxMB_population.xlsx")
    df = clean_phenotype_data(df)
    df["Plant"] = df["Plant"].astype(str)
    result = df[~(df["Plant"].str.contains("M")) & ~(df["Plant"].str.contains("N"))]
    print(len(result))