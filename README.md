# Thesis_Paper_1: "Ecological trade-offs as a product of water saving irrigation strategies in rice fields"
1. R Project: Thesis_Paper_1_git 
2. /src: R scripts.
3. /data: data inputs (e.g. .csv, .xlsx)
4. /outputs: analysis outputs (e.g., plots, .csv, etc.)

## /src: R scripts:

### Macroinv_2022: Working with biodiversity data. 
- Field sampling data in "data/BIO/Macroinv_2022.csv". 
- Taxonomic corrections afer collaboration with Dani Boix (UdG).
- Taxonomic assignations: Assigns Family_SubFamily, Genus and Species of individuals identified until previous taxonomical level in regard to those identified up to this assignation level (assignation is weighted according to the proportion of identified individuals and only considers possible Family_SubFamily-Genus-Species combinations). 
- Diversity Analysis - Hill Numbers: Calculates diversity metrics according to Hill Numbers using iNEXT package (Chao et al., 2014).
- Plotting diversity metrics.
- Creating dataframe merging macroinvertebrate and fish-amphibean data (latter from: "data/BIO/Macrofauna_2022.csv").
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



