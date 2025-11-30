# GDAR

# Code associated with: Variability, drivers, and utility of genetic diversity-area relationships in terrestrial vertebrates 
## Schmidt et al. 2025 EcoEvoRxiv (doi.org/10.32942/X2Q34P)

File descriptions: <br>
- <b>1_GDAR_script_git.R:</b> generate artificial area samples to construct and plot GDARs and estimate scaling exponents (z) <br>
- <b>2_ interaction_models_git.R:</b> fit and plot models relating genetic diversity to area, FST, and their interaction <br>
- <b>3_trait_data_git.R:</b> compile species trait data <br>
- <b>4_training_test_git.R:</b> test the predictability of FST and scaling exponents z using species traits <br>
- <b>functions_git.R:</b> functions needed to run analyses in other scripts <br>
- <b> GDAR_env_figs.Rdata:</b> R environment object including data needed to run scripts <br>
- <b> zfst_long_se.rds:</b> a saved R object produced in script 1 <br>
- <b> mammal_traits.csv:</b> trait data for species included in this analysis <br>
   Column descriptions:<br>
  - filename: name of file containing genotypes <br>
  - species: Genus_species<br>
  - individuals: number of individuals sampled<br>
  - num_sites: number of sites (locations) sampled<br>
  - global_fst: estimated global FST (population differentiation)<br>
  - zAC, zAR, zGD: estimated scaling values (z) for allele count (AC), allelic richness (AR), and gene diversity (GD)<br>
  - zACse, zARse, zGDse: standard errors for scaling values (z)	<br>
  - class, uid, family, genus, order: other taxonomic information	<br>
  - mass_g, log_mass: species mass in grams, and log10 of this value<br>
  - range_area_km2, log_range: species range area in square kilometers, and log10 of this value<br>
  - mean_HR_km2, log_HR: species mean home range area in square kilometers, and log10 of this value<br>
  - area_km2, log_area: area of sample extent (per study) in square kilometers, and log10 of this value<br>
- <b>mpg_test_data.csv:</b> trait data for mammal species from the MacroPopGen database <br>
   Column descriptions: <br>
  - species_ref: species binomial and MacroPopGen reference ID <br>
  - species: Genus_species <br>
  - order, family, genus: other taxonomic information <br>
  - RefID: MacroPopGen reference ID<br>
  - global_fst: global FST value per study/species combination<br>
  - mass_g, log_mass: species mass in grams, and log10 of this value<br>
  - range_area_km2, log_range: species range area in square kilometers, and log10 of this value<br>
  - mean_HR_km2, log_HR: species mean home range area in square kilometers, and log10 of this value<br>
  - area_km2, log_area: area of sample extent (per study) in square kilometers, and log10 of this value<br>
