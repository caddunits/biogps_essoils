#
# Script to prepare the ATCC mapping
# Data is read from original CSV files and saved as .RDS object
#
# Below we describe the procedure to extract data from the ATCC website
# 1. Login into website ATCC
# 2. At the page https://genomes.atcc.org insert the ATCC code in the field
# "Search for a genome"
# 3. Click on the "View" button in correspondence to the column "Genomic Data"
# 4. Click on the tab "Genome Browser"
# 5. Click on the button "Download table CSV"
# Note that the uniprot codes are in the J column of the downloaded files
# 

library(dplyr)
library(magrittr)

repo_dir <- getwd()

# Input
data_dir <- file.path(repo_dir, "data")
atcc_dir <- file.path(data_dir, "ATCC_data")

# Output
mapp_dir <- file.path(data_dir, "Mapping_data")
atcc_mapping <- file.path(mapp_dir, "atcc_bacteria_uniprot.RDS")

atcc_correspondences <- data.frame(
  bact_atcc = c(
    '27853',
    '700603',
    '29212',
    '35984',
    '29213',
    '25922'
  ),
  bacteria = c(
    'Pseudomonas aeruginosa',
    'Klebsiella pneumoniae',
    'Enterococcus faecalis',
    'Staphylococcus epidermidis',
    'Staphylococcus aureus',
    'Escherichia coli'
  ),
  stringsAsFactors = FALSE
)


atcc_files <- list.files(atcc_dir, pattern=".csv$", full.names=T)
cat("READ DATA FROM THE FOLLOWING CSV FILES:\n")
print(atcc_files)


if (length(atcc_files) > 0) {
  atcc_bacteria_uniprot_df <- do.call('rbind', lapply(atcc_files, function(f) {
    atcc_code <- gsub(pattern='ATCC_|.csv', replacement='', x=basename(f))
    read.delim2(file=f, sep = ",", header = TRUE) %>% 
      dplyr::rename('uniprot'='Uniprot.ID', 'genename'='Name',
                    'description'='Product', 'aa_sequence'='Amino.Acid.Sequence') %>%
      dplyr::select(uniprot, genename, description, aa_sequence) %>%
      dplyr::distinct() %>%
      dplyr::mutate(bact_atcc=atcc_code)
    
  })) %>% 
    dplyr::left_join(atcc_correspondences, by="bact_atcc") %>% 
    dplyr::filter(description != "hypothetical protein",
                  uniprot != "",
                  genename != "") %>% 
    dplyr::select(-aa_sequence)
  
  # Save as .RDS object
  saveRDS(object=atcc_bacteria_uniprot_df, file=atcc_mapping)
  
} else {
  cat("\nSome problems in the input files\n\n")
}
