#
# Script to prepare the read data about bacteria. 
# Thew XLSX file was manually curated
# 

library(dplyr)
library(magrittr)


repo_dir <- getwd()

# Input
data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")
bacteriadata_filename <- file.path(data_dir, "BACTERIA_data", "BACTERIADATA.xlsx")

# Output
bacteria_data <- file.path(mapp_dir, "bacteria_data.RDS")

wanted_bacteria <- c(
  "Staphylococcus aureus",
  "Staphylococcus epidermidis",
  "Enterococcus faecalis",
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Pseudomonas aeruginosa"
)



# Input info on bacteria include codes on different databases (ATCC, KEGG)
df_bacteria <- readxl::read_excel(path=bacteriadata_filename) %>%
  dplyr::filter(bacteria %in% wanted_bacteria)


saveRDS(object=df_bacteria, file=bacteria_data)

