import pandas as pd
import scipy.stats as stats
import sys

def get_average_per_plant(df: pd.DataFrame, series_name: str) -> pd.DataFrame:
    # groups series, then gets average by grouped.
    return (
        df.groupby("Plant")[series_name].mean().reset_index()
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
    return df.loc[~mask_failed].reset_index(drop=True)

def t_test_M(df: pd.DataFrame) -> pd.DataFrame:
    # returns dataframe where all non parental plants are compared to the M parental plant and a indipended t-test is performed

    df = df.copy()
    
    df["Plant"] = df["Plant"].astype(str)
    df_m = df[(df["Plant"].str.contains("M"))].reset_index()

    df_other = df[~(df["Plant"].str.contains("M")) & ~(df["Plant"].str.contains("N"))].reset_index()

    df_other["Plant"] = df_other["Plant"].astype(int)

    grouped = df_other.groupby("Plant").agg(
        p_value=("Pollen Viability",
                 lambda plant: stats.ttest_ind(plant, df_m["Pollen Viability"]).pvalue), 
    )

    return grouped.reset_index().sort_values("Plant").reset_index(drop=True)

def merge_geno_and_fenotype(df_geno: pd.DataFrame, df_feno: pd.DataFrame) -> pd.DataFrame:
    df_geno = df_geno.copy()

    df_feno_minimal = df_feno[["Plant", "Pollen Viability"]]

    df_feno_minimal = get_average_per_plant(df_feno_minimal, "Pollen Viability")
    
    marker_info = df_geno[["Markers", "Chromosome", "location"]]

    df_geno_transposed = (
        df_geno
        .drop(columns=["Plant", "Chromosome", "location"])
        .set_index("Markers")
        .T
        .rename_axis("Plant")
    .   reset_index()
    )

    df_geno_transposed["Plant"] = df_geno_transposed["Plant"].astype(str).str.strip()
    df_feno_minimal["Plant"] = df_feno_minimal["Plant"].astype(str).str.strip()

    return pd.merge(df_geno_transposed, df_feno_minimal, on="Plant", how="inner"), marker_info

def check_merge(df_merged: pd.DataFrame):
    #checks merges
    print(df_merged.shape)
    print(df_merged.head())
    print(df_merged["Plant"].duplicated().sum())
    print(df_merged.isna().sum())

def group_by_location(df: pd.DataFrame, phenotype:str="Pollen Viability") -> list:
    df = df.copy()

    snp_cols = []
    for col in df.columns:
        if col != "Plant" and col != phenotype and col != "Unnamed: 0":
            snp_cols.append(col)

    result = []

    for col in snp_cols:
        result.append({
            "position": col,
            "df": df[[col,phenotype]]
        })

    return result

def get_p_values(data_frames: list, phenotype:str="Pollen Viability") -> pd.DataFrame:
    positions = []
    p_values = []
    
    for item in data_frames:
        position = item["position"]
        df = item["df"]
        unique = df[position].unique()

        grouped = df.groupby(df[position])
        groups = []

        if len(unique) < 2:
            continue
            
        for combination in unique:
            groups.append(grouped.get_group(combination)[phenotype])

        p_value = stats.f_oneway(*groups).pvalue

        positions.append(position)
        p_values.append(p_value)

    return pd.DataFrame({
        "position": positions,
        "p_value": p_values
        })

if __name__ == "__main__":
    df = pd.read_csv("data/merged_df.csv")
    result = group_by_location(df)
    print(len(result))
    p_values = get_p_values(result)
    print(p_values)