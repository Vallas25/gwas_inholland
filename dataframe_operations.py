import pandas as pd
import scipy.stats as stats
import numpy as np

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

def filter_failed_homozygote(df: pd.DataFrame) -> pd.DataFrame:
    # returns df where all failed mesurements and hetrozygote loci are removed
    df = df.copy()

    # names of allele columns
    allele_columns = ["Parent1.1", "Parent1.2", "Parent2.1", "Parent2.2"]

    # masks for removing rows containing a failed mesurement or both parents are heterozygotes
    mask_homozygote = (
        df[allele_columns].eq(df["Parent1.1"], axis=0).all(axis=1)
    )

    mask_failed = (
        df.eq("failed").any(axis=1)
    )

    # returning dataframe with only homologous parental markers and non failed mesurements
    df = df.loc[~mask_homozygote]
    return df.loc[~mask_failed].reset_index(drop=True)

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

def merge_geno_and_fenotype(df_geno: pd.DataFrame, df_pheno: pd.DataFrame) -> pd.DataFrame:
    # returns two dataframes where the geno and pheno type are merged and another with extra information from the genotype dataframe
    df_geno = df_geno.copy()

    # creating a minimal version of the phenotype df and grouping by plant
    df_pheno_minimal = df_pheno[["Plant", "Pollen Viability"]]

    df_pheno_minimal = get_average_per_plant(df_pheno_minimal, "Pollen Viability")
    
    # creating a seperate df for the marker info
    marker_info = df_geno[["Markers", "Chromosome", "location"]]

    # transposing the geno df (switching columns and rows)
    df_geno_transposed = (
        df_geno
        .drop(columns=[ "Chromosome", "location"])
        .set_index("Markers")
        .T
        .rename_axis(index="Plant", columns=None)
        .reset_index()
    )

    print("Transposed")

    # formatting everything correctly, so merge will go correctly
    df_geno_transposed["Plant"] = df_geno_transposed["Plant"].astype(str).str.strip()
    df_pheno_minimal["Plant"] = df_pheno_minimal["Plant"].astype(str).str.strip()

    # merging
    return pd.merge(df_geno_transposed, df_pheno_minimal, on="Plant", how="inner"), marker_info

def check_merge(df_merged: pd.DataFrame):
    #checks merges
    print(df_merged.shape)
    print(df_merged.head())
    print(df_merged["Plant"].duplicated().sum())
    print(df_merged.isna().sum())

def group_by_location(df: pd.DataFrame, phenotype:str="Pollen Viability") -> list:
    #gets SNP locations and returns list containing dictionary, wich contains position and df with geno and phenotype

    # getting vallid SNP columns
    snp_cols = []
    for col in df.columns:
        if col != "Plant" and col != phenotype and col != "Unnamed: 0":
            snp_cols.append(col)

    # getting list containing the location and data for said location
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
        unique = list(df[position].unique())

        #grouping by basepairs
        grouped = df.groupby(df[position])
        groups = []

        if "failed" in unique:
            unique.remove("failed")

        # if there are less than 2 unique genotypes, annova cannot be performed
        if len(unique) < 2:
            continue
            
        # makes list for every base pair combination in data
        for combination in unique:
            groups.append(grouped.get_group(combination)[phenotype])

        # performs ANOVA
        p_value = stats.f_oneway(*groups).pvalue

        positions.append(position)
        p_values.append(p_value)

    # creates dataframe containing info
    return pd.DataFrame({
        "Markers": positions,
        "p_value": p_values
        })

def merge_pvalues_and_genotype_info(p_values:pd.DataFrame, genotype_info:pd.DataFrame) -> pd.DataFrame:
    # merges p-values and genotype info dataframes
    merged = pd.merge(
        p_values,
        genotype_info,
        how="inner",
        on="Markers",
        indicator=True
    )
    return merged

def prep_df_for_plotting(df: pd.DataFrame) -> pd.DataFrame:
    df["-logp"] = -np.log10(df["p_value"])
    df["chr_num"] = df["Chromosome"].str.extract(r"ch(\d+)$")[0].astype(int)

    grouped = df.groupby("chr_num")["location"].max()

    chr_offset = grouped.cumsum().shift(1, fill_value=0)

    df["location_offset"] = (df["location"] + df["chr_num"].map(chr_offset))

    return df

if __name__ == "__main__":
    df = pd.read_csv("data/merged_df.csv")
    result = group_by_location(df)
    print(len(result))
    p_values = get_p_values(result)
    print(p_values)