import pandas as pd
import scipy.stats as stats

def get_average_per_plant(df: pd.DataFrame, series1_name: str, series2_name: str) -> pd.DataFrame:
    # groups series, then gets average by grouped.
    return (
        df.groupby(series1_name)[series2_name].mean().reset_index()
    )

def clean_phenotype_data(df: pd.DataFrame) -> pd.DataFrame:
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
    print(df is None)
    df = df.copy()
    # use mask with .str.contains("M") instead of this convoluted way
    df = df.sort_values(
        ["Plant"],
        key=lambda col: col.astype(str).map(
            lambda x: (0, int(x)) if x.isdigit() else (1,x)
        )
        ).reset_index()
    df_m = df.iloc[646:685]
    df_other = df.iloc[:645] # use mask ~.str.contains("M")|~.str.contains("N")

    grouped = df_other.groupby("Plant").agg(
        p_value="Pollen Viability", 
        lambda plant: stats.ttest_ind(plant, df_other).pvalue 
    )
    return grouped.reset_index()


if __name__ == "__main__":
    df = pd.read_excel("data/Phenotype_data_PennxMB_population.xlsx")
    df = clean_phenotype_data(df)
    result = t_test_M(df)
    print(result)