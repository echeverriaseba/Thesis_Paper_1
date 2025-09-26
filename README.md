# Thesis_Paper_1: "Water-Saving Irrigation Can Mitigate Climate Change But Entails Negative Side Effects on Biodiversity in Rice Paddy Fields"
Open access article: [Agriculture, Ecosystems and Environment](https://doi.org/10.1016/j.agee.2025.109719)

Contact: secheverriap@gmail.com (do not hesitate to contact in case you have any doubts using this material)

## Index:
1. R Project: Thesis_Paper_1_git 
2. /src: R scripts.
3. /data: data inputs (e.g. .csv, .xlsx). All data colected across the sampling campaign 2022.
4. /outputs: analysis outputs (e.g., plots, .csv, etc.)

## /src: R scripts
Description of the main scripts used for data preparation, plotting and statistical analyses: 

### Macroinv_2022: Working with biodiversity data. 
- Field sampling data in "/data/BIO/Macroinv_2022.csv". 
- Taxonomic corrections afer collaboration with Prof. Dani Boix (UdG).
- Taxonomic assignations: Assigns Family_SubFamily, Genus and Species of individuals identified until previous taxonomical level in regard to those identified up to this assignation level (assignation is weighted according to the proportion of identified individuals and only considers possible Family_SubFamily-Genus-Species combinations). 
- Diversity Analysis - Hill Numbers: Calculates diversity metrics according to Hill Numbers using iNEXT package (Chao et al., 2014), implementing iNEXT::iNEXT() within a new "extrar_divmetrics" function that calculates and filters for diversity metrics: "Species richness", "Shannon diversity" and "Simpson diversity". The function calculates metrics for each siteID (Plot_Sampling_Treat) data subset.
- Plotting diversity metrics. Plot files: "outputs/Plots". 
- Creating dataframe merging macroinvertebrate and fish-amphibian data (latter from: "/data/BIO/Macrofauna_2022.csv"). Merged dataframe: "Abundance_2022"
- Plotting accumulated abundance. Plot files: "outputs/Plots". 

### Macroinv_2022_stats: Statistical analyses for biodiversity data.
- Preparing biodiversity-physicochemical_data dataframe ("Hills_Physchem").
- Data validation: According to "Analysing the impact of multiple stressors in aquatic biomonitoring data: A ‘cookbook’ with applications in R" - Feld et al., 2016.
   - Outlier analysis.
   - Testing correlations: Spearman rank correlation, Variance inflation factors (VIF) and Data Dendrogram
 - GLMMs:
   - Checking Hill metrics vs Samplings: Testing for linear or quadratic relation. Cuadratic relation identified for q0 and q1.
   - Testing models: Testing for best fitting model for q0 and q1 considering different variables, with/without outliers, different combinations of Treat*Sampling interactions, different disfribution families, etc. This script contains only selected models, for all alternative tested models see. script "Macroinv_2022_stats_alt_models".
   - q0: Post-hoc tests do not proceed due to signifficant effect of Treat-Sampling2 in selected model. Post-hoc tests not applicable for Factor-Contiuous_variable interactions. Plotting q0 vs Treat for each Sampling using empirical (not predicted) data.
   - q1:  Post-hoc tests do proceed because of a signifficant effect of Treat variable. Post-hoc tests using emmeans package to test differences among Treat.
   - Abundance models.

### Macroinv_2022_stats_alt_models: All alternative tested models.
  
### GHG_rates_2022: Calculating original CH<sub>4</sub>emission rates from sampled GHG concentration through gas chromatography.
- Creating "Chrom_results_2022" dataframe, which calculates GHG concentrations in mg m<sup>-2</sup> from C-ppm chromatography results.
- Calculating emission rates through lm(), also R<sup>2</sup> for each lm() to apply posterior model corrections (see "Emission_rates_2022" dataframe).

### GHG_rates_2022_w_corrections: Applying corrections to GHG (CH<sub>4</sub>, N<sub>2</sub>O and CO<sub>2</sub>) emission rates.   
- Same procedure as "GHG_rates_2022" script but fitting 4 alternative "3-values" models (each one removing one time-step) and then selecting that which achieves higher R<sup>2</sup> and positive rate (lm slope).
- Output results saved in "/outputs/GHG/2022/Rates_corrected/Emission_rates_w_corrections_2022.RData" and in "/outputs/GHG/2022/Rates_corrected/Results_CH4_w_corrections.xlsx". 

### Emissions_2022: Working with GHG emission rates data.
- Creating "Master_GHG_2022" dataframe. Emission rates were previously calculated and corrected in script "GHG_rates_2022_w_corrections" and loaded from "/outputs/GHG/2022/Rates_corrected/Emission_rates_w_corrections_2022.RData". This dataframe merges GHG emissions, physicochemical and water level data. Water level is "corrected", meaning that it takes the closest level data to the GHG sampling (see script for notes on script for correction conditions).
- Plot Analyses: plotting Emissions vs water level. See outputs in "/outputs/Plots/GHG")

### Emissions_2022_stats: Statistical analyses for GHG emission rates data.
- Creates a merged dataframe for calculated fluxes and field physicochemical data.
- Data validation: According to "Analyzing the impact of multiple stressors in aquatic biomonitoring data: A ‘cookbook’ with applications in R" - Feld et al., 2016:
  - Outlier analysis.
  - Testing correlations:
    - Spearman rank correlation.
    - Variance inflation factors (VIF).
    - Data Dendrogram.
- Fitting GLMMs.   
- Other analyses:
  - Comparing overall emission averages per Treat.
  - Calculating global cumulative C-CH4 emissions (kg C-CH4 ha-1).
  - Model only for CON - Phys. drivers analysis.

### Prod_2022: Working with yield data.
- Harvest data from "/data/PROD/Prod_2022.csv".
- Plotting yield vs Treat.

### Prod_2022_stats: Statistical analyses for yield data.
- Linear model to test treatment effect.
