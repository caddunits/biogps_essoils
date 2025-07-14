# biogps_essoils
BioGPS application for phytocomplexes from essential oils

### Prepare Data

#### Data mapping Uniprot <-> ATCC
1. Run the script step1_prepare_data_atcc_mapping.R

This script is used to create a file named "atcc_bacteria_uniprot.RDS"
which is stored in the directory data/Mapping_data

NOTE: according to ATCC website, the download of CSV files is available only to 
"ATCC Genome Portal supporting members or those who have purchased a 
 corresponding physical product".
Therefore, in the official repository of the current project we cannot put
CSV files, but only a modified version, as .RDS file with only the necessary 
columns for running the other scripts. 

However, We keep in the repository the script to create the file RDS from the
set of CSV files.



