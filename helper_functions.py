import pandas as pd
import scipy.stats as stats
import numpy as np
from plot_operations import create_genomic_heatmap
from dataframe_operations import *

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

def serialise_genotype(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    cols = df.columns.to_list()
    df_serialised = df[["Chromosome", "location", "Markers"]]
    to_remove = ["Parent1.1", "Parent1.2", "Parent2.1", "Parent2.2", "Chromosome", "location", "Markers"]
    for item in to_remove:
        cols.remove(item)

    encoded = {}

    for col in cols:
        homo_p1 = (
            (df[col] == df["Parent1.1"]) &
            (df[col] == df["Parent1.2"])
        )

        homo_p2 = (
            (df[col] == df["Parent2.1"]) &
            (df[col] == df["Parent2.2"])
        )

        hetero = ~(homo_p1 | homo_p2)

        encoded[col] = np.select(
            [homo_p1, homo_p2, hetero],
            [0, 1, 2]
        )

    df_serialised = pd.concat(
        [df_serialised.drop(columns=cols, errors="ignore"),
        pd.DataFrame(encoded, index=df.index)],
        axis=1
    )

    return df_serialised

def sort_cols (df_geno: pd.DataFrame, df_pheno: pd.DataFrame) -> pd.DataFrame:
    metadata = ["Chromosome", "location", "Markers"]

    df_pheno["Plant"] = df_pheno["Plant"].astype(str)
    df_pheno = df_pheno[~(df_pheno["Plant"].str.contains("M|N"))]
    df_pheno["Plant"] = df_pheno["Plant"].astype(int)

    grouped = df_pheno.groupby("Plant")["Pollen Viability"].mean()
    
    order = grouped.sort_values().index.astype(str).to_list()
    
    return df_geno.loc[:, metadata + order]

def sort_for_plotting(df: pd.DataFrame) -> pd.DataFrame:

    df["chr_num"] = df["Chromosome"].str.extract(r"ch(\d+)$")[0].astype(int)

    return df.sort_values(
        by=["chr_num", "location"],
        ascending=[True,True]
    ).reset_index(drop=True)

if __name__ == "__main__":
    df = pd.read_excel("data/Sequencing_data_Inholland.xlsx")
    df_pheno = pd.read_excel('data/Phenotype_data_PennxMB_population.xlsx')
    df_pheno = clean_phenotype_data(df_pheno)
    df_clean = filter_failed_homozygote(df)
    serialised = serialise_genotype(df)
    sorted = sort_cols(df_geno=serialised, df_pheno=df_pheno)
    plotable= sort_for_plotting(sorted)
    print(plotable)
    plot = create_genomic_heatmap(plotable)
    plot.write_html("plots/heatmap.html")
    