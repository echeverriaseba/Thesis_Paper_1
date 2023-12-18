# Thesis_Paper_1: "Ecological trade-offs as a product of water saving irrigation strategies in rice fields"
1. R Project: Thesis_Paper_1_git 
2. /src: R scripts.
3. /data: data inputs (e.g. .csv, .xlsx)
4. /outputs: analysis outputs (e.g., plots, .csv, etc.)

## /src: R scripts:
Here the main scripts used for data preparation, plotting and statistical analyses are described. Other scripts can be found in this folder containing coding tests. 

### Macroinv_2022: Working with biodiversity data. 
- Field sampling data in "/data/BIO/Macroinv_2022.csv". 
- Taxonomic corrections afer collaboration with Dani Boix (UdG).
- Taxonomic assignations: Assigns Family_SubFamily, Genus and Species of individuals identified until previous taxonomical level in regard to those identified up to this assignation level (assignation is weighted according to the proportion of identified individuals and only considers possible Family_SubFamily-Genus-Species combinations). 
- Diversity Analysis - Hill Numbers: Calculates diversity metrics according to Hill Numbers using iNEXT package (Chao et al., 2014).
- Plotting diversity metrics.
- Creating dataframe merging macroinvertebrate and fish-amphibean data (latter from: "/data/BIO/Macrofauna_2022.csv").
- Plotting accumulated abundance.

### Macroinv_2022_stats: Statistical analyses for biodiversity data.
- Preparing biodiversity-physicochemical_data dataframe ("Hills_Physchem").
- Data validation: According to "Analysing the impact of multiple stressors in aquatic biomonitoring data: A ‘cookbook’ with applications in R" - Feld et al., 2016.
   - Outlier analysis.
   - Testing correlations: Spearman rank correlation, Variance inflation factors (VIF) and Data Dendrogram
 - GLMMs:
   - Checking Hill metrics vs Samplings: Testing for linear or quadratic relation. Cuadratic relation identified for q0 and q1.
   - Testing models: Testing for best fitting model for q0 and q1 considering different variables, with/without outliers, different combinations of Treat*Sampling interactions, different disfribution families, etc.
   - q0: Post-hoc tests do not proceed due to signifficant effect of Treat-Sampling2 in selected model. Post-hoc tests not applicable for Factor-Contiuous_variable interactions. Plotting q0 vs Treat for each Sampling using empirical (not predicted) data.
   - q1:  Post-hoc tests do proceed because of a signifficant effect of Treat variable. Post-hoc tests using emmeans package to test differences among Treat.
  
### GHG_rates_2022: Calculating original CH<sub>4</sub>emission rates from sampled GHG concentration through gas chromatography.
- Creating "Chrom_results_2022" dataframe, which calculates GHG concentrations in mg m<sup>-2</sup> from C-ppm chromatography results.
- Calculating emission rates through lm(), also R<sup>2</sup> for each lm() to apply posterior model corrections (see "Emission_rates_2022" dataframe).

### GHG_rates_2022_w_corrections: Applying corrections to GHG (CH<sub>4</sub>, N<sub>2</sub>O and CO<sub>2</sub>) emission rates.   
- Same procedure as "GHG_rates_2022" script but fitting 4 alternative "3-values" models (each one removing one time-step) and then selecting that which achieves higher R<sup>2</sup> and positive rate (lm slope).
- Output results saved in "/outputs/GHG/2022/Rates_corrected/Emission_rates_w_corrections_2022.RData" and in "/outputs/GHG/2022/Rates_corrected/Results_CH4_w_corrections.xlsx". 

### Emissions_2022: Working with GHG emission rates data.
- Creating "Master_GHG_2022" dataframe. Emission rates were previously calculated and corrected in script "GHG_rates_2022_w_corrections" and loaded from "/outputs/GHG/2022/Rates_corrected/Emission_rates_w_corrections_2022.RData". This dataframe merges GHG emissions, physicochemical and water level data. Water level is "corrected", meaning that it takes the closest level data to the GHG sampling (see script for notes on script for correction conditions).
- Plot Analyses: plotting Emissions vs water level. See outputs in "/outputs/Plots/GHG")

### Prod_2022: Working with yield data.
- Harvest data from "/data/PROD/Prod_2022.csv".
- Plotting yield vs Treat. 
