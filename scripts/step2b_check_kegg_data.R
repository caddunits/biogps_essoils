
library(dplyr)
library(magrittr)
library(writexl)
library(tidyr)

repo_dir <- getwd()

# Settings are set once, in a configuration file
source(file.path(repo_dir, "scripts", "config.R"))


# Input
data_dir <- file.path(repo_dir, "data")
kegg_dir <- file.path(data_dir, "KEGG_data")

# Output
output_dir <- file.path(repo_dir, "output", "kegg_details")

# Create the ouput directory if it does not exist
dir.create(output_dir, showWarnings = FALSE)

# 
# Read from each subdirectory all the .RDS files 
# 

kegg_files <- list.files(kegg_dir, pattern=".RDS$", recursive=T, full.names=T)


kegg_info_df <- do.call('rbind', lapply(kegg_files, function(f) {
  kegg_data_singlefile <- readRDS(file = f)
  single_df <- data.frame(
    bact_dir = basename(dirname(f)),
    bact_code = gsub(pattern='.RDS$', replacement='', x=basename(f)),
    nr_uniprot = length(unique(kegg_data_singlefile$uniprot[!is.na(kegg_data_singlefile$uniprot)])),
    nr_NCBI = length(unique(kegg_data_singlefile$NCBI[!is.na(kegg_data_singlefile$NCBI)])),
    nr_pathways = length(unique(kegg_data_singlefile$pw_value[!is.na(kegg_data_singlefile$pw_value)])),
    stringsAsFactors = FALSE
  )
}))

writexl::write_xlsx(x = kegg_info_df, path = file.path(output_dir, "kegg_details.xlsx"))
