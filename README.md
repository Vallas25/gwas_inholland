# GWAS

## tools
- vscode
- wsl (ubuntu)
- plotly
- pandas
- openpyxl
- numpy
- scipy

## data
De bestanden [Phenotype_data_PennxMB_population.xlsx](/data/Phenotype_data_PennxMB_population.xlsx) en [Sequencing_data_Inholland.xlsx](/data/Sequencing_data_Inholland.xlsx) horen bij elkaar.

## code
De functie ```get_average_per_plan()``` accepteerd een dataframe en de naam van twee series(colom namen). Deze functie maakt vervolgens een dataframe met series 1 gegroupeerd (planten) en het gemiddelde van series 2

De functie ```clean_phenotype_data()``` accepteerd een dataframe (phenotype data) en haalt de colommen weg die geen nuttige informatie hebben.

De functie ```filter_failed_hetrozygote()``` accepteerd een dataframe (sequencing_data) en haalt de rijen weg waarbij alle loci van de ouders (Parent) heterozygoot zijn, maar ook als er maar 1 van de readings gefaald (failed) is.
De heterozygoot test is [hier](/data/heterozygote_test.csv) te vinden en de heterozygoot + failed test is [hier](/data/failed_test.csv) te vinden.