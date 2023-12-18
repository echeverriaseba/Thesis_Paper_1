# Thesis_Paper_1: "Ecological trade-offs as a product of water saving irrigation strategies in rice fields"
1. R Project: Thesis_Paper_1_git 
2. /src: R scripts.
3. /data: data inputs (e.g. .csv, .xlsx)
4. /outputs: analysis outputs (e.g., plots, .csv, etc.)

## /src: R scripts:

### Macroinv_2022: Working with macroinvertebrate data. 
- Field sampling data in "data/BIO/Macroinv_2022.csv". 
- Taxonomic corrections afer collaboration with Dani Boix (UdG).
- Taxonomic assignations: Assigns Family_SubFamily, Genus and Species of individuals identified until previous taxonomical level in regard to those identified up to this assignation level (assignation is weighted according to the proportion of identified individuals and only considers possible Family_SubFamily-Genus-Species combinations). 
- Diversity Analysis - Hill Numbers: Calculates diversity metrics according to Hill Numbers using iNEXT package (Chao et al., 2014).
- Plotting diversity metrics.
- Plotting accumulated abundance.

