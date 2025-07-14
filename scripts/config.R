# 
# configuration 
# 
config <- list(
  debug = FALSE,
  # debug = TRUE,
  verbose = FALSE,
  wanted_bacteria = c(
    "Staphylococcus aureus",
    "Staphylococcus epidermidis",
    "Enterococcus faecalis",
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa"
  ),
  major_bacteria = c(
    "Staphylococcus aureus",
    "Escherichia coli",
    "Pseudomonas aeruginosa"
  ),
  minor_bacteria = c(
    "Enterococcus faecalis",
    "Klebsiella pneumoniae",
    "Staphylococcus epidermidis"
  ),
  # STRINGdb data
  species_mapping = c(
    "Staphylococcus aureus" = 93061,
    "Staphylococcus epidermidis" = 176279,
    "Enterococcus faecalis" = 1260356,
    "Escherichia coli" = 199310,
    "Klebsiella pneumoniae" = 272620,
    "Pseudomonas aeruginosa" = 208964
  ),
  # Parameters
  nr_best_pockets = 100,
  nr_permutations = 1000,
  threshold_comp = 0.05,
  threshold_zzscore = 0.00,
  threshold_zscore_mol = 0.00,
  threshold_zscore_poc = 0.00,
  threshold_zzscore_barplot = 2.00,
  threshold_zscore_mol_barplot = 1.00,
  threshold_zscore_poc_barplot = 1.00,
  threshold_nrP = 1000,
  threshold_pvalue = 0.05
)

# NOTE: in these data we mostly have only two classes,
# oxygenated_monoterpenoids and in minor extent, phenylpropanoids
# So, we assign for all the same colour, but the script is ready
# to assign different colors to different phytochemical classes
colors_phytoclasses <- list(
  oxygenated_monoterpenoids = "#AEE7F8",
  monoterpenoids = "#AEE7F8",
  sesquiterpenes = "#AEE7F8",
  phenylpropanoids = "#AEE7F8",
  other_organic_compounds = "#FAB972"
)

