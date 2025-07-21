#
# Script to prepare a table with the summary of bacterial data, with information
# from KEGG, ATCC and BioGPS.
# Data is read from original .RDS files and saved as unique .CSV file
# which corresponds to the table for the paper
# It can only work if the file info.sql is available (that was provided by
# Molecular Discovery Ltd)
# 

library(dplyr)
library(magrittr)
library(DBI)
library(RSQLite)

# Change this line according to the path of your repository
repo_dir <- getwd()

# Input
data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")
kegg_bacteria_unified_df <- readRDS(file = file.path(mapp_dir, "kegg_bacteria_unified.RDS"))
atcc_bacteria_uniprot_df <- readRDS(file = file.path(mapp_dir, "atcc_bacteria_uniprot.RDS"))

# NOTE: This database was gently provided by Molecular Discovery Ltd
db_dir <- file.path(repo_dir, "db")
dbfilename_biogps_info <- file.path(db_dir, "info.sql")

wanted_bacteria <- c(
  "Staphylococcus aureus",
  "Staphylococcus epidermidis",
  "Enterococcus faecalis",
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Pseudomonas aeruginosa"
)


# Output
finaltable_filename <- file.path(repo_dir, "output", "finaltable2.csv")

#
# Functions to extract data from BioGPS database 
#
# Note: in some cases it might be convenient to filter the results 
# --AND dat.exp_method = 'X-RAY DIFFRACTION' (Filter not used in this project)
# --AND dat.resolution <= 2.5 (Filter not used in this project)
# 

extract_biogps_pockets <- function(con, organism_name) {
  
  query <- paste0("SELECT DISTINCT poc.inchi
, cha.protein_code
, dat.resolution
, dat.organism
, dat.exp_method
 FROM pockets poc
 LEFT JOIN pocket_chains cha 
   ON cha.pocket_inchi = poc.inchi
 LEFT JOIN chain_data dat 
   ON dat.protein_code = cha.protein_code 
   AND dat.chain_code = cha.chain_code
WHERE dat.organism LIKE '%", organism_name, "%'")
  
  df <- dbGetQuery(con, query)
  return(df)
}


extract_biogps_pdb <- function(con, organism_name) {
  
  query <- paste0("SELECT DISTINCT cha.protein_code
, dat.resolution
, dat.organism
, dat.exp_method
 FROM pockets poc
 LEFT JOIN pocket_chains cha 
   ON cha.pocket_inchi = poc.inchi
 LEFT JOIN chain_data dat 
  ON dat.protein_code = cha.protein_code 
  AND dat.chain_code = cha.chain_code
WHERE dat.organism LIKE '%", organism_name, "%'")
  
  df <- dbGetQuery(con, query)
  return(df)
}


extract_biogps_uniprot <- function(con, organism_name) {
  
  query <- paste0("SELECT DISTINCT uni.uniprot
, dat.organism
 FROM pockets poc
 LEFT JOIN pocket_chains cha 
   ON cha.pocket_inchi = poc.inchi
 LEFT JOIN chain_data dat 
   ON dat.protein_code = cha.protein_code 
   AND dat.chain_code = cha.chain_code
 LEFT JOIN chain_uniprot uni 
   ON uni.protein_code = cha.protein_code 
   AND uni.chain_code = cha.chain_code
 WHERE dat.organism LIKE '%", organism_name, "%'")
  
  df <- dbGetQuery(con, query)
  return(df)
}


extract_bacteria_data_pockets <- function(con, bact_codes) {
  
  bact_codes_string <- paste0(bact_codes, collapse = "', '")
  query <- sprintf("SELECT DISTINCT  
  poc.name 
  , pdb.pdb
  , uni.uniprot
  , org.organism
  FROM pockets poc
  LEFT JOIN pdb ON pdb.id = poc.id_pdb
  LEFT JOIN mapping_pdb_uniprot mpu ON mpu.id_pdb = pdb.id
  LEFT JOIN uniprot uni ON uni.id = mpu.id_uniprot
  LEFT JOIN organisms org ON org.id = uni.id_organism
  WHERE org.organism IN ('%s')", bact_codes_string)
  
  df <- dbGetQuery(con, query)
  
  return(df)
}


extract_bacteria_data_pdb <- function(con, bact_codes) {
  
  bact_codes_string <- paste0(bact_codes, collapse = "', '")
  query <- sprintf("SELECT DISTINCT  
  pdb.pdb
  , uni.uniprot
  , org.organism
  FROM pockets poc
  LEFT JOIN pdb ON pdb.id = poc.id_pdb
  LEFT JOIN mapping_pdb_uniprot mpu ON mpu.id_pdb = pdb.id
  LEFT JOIN uniprot uni ON uni.id = mpu.id_uniprot
  LEFT JOIN organisms org ON org.id = uni.id_organism
  WHERE org.organism IN ('%s')", bact_codes_string)
  
  df <- dbGetQuery(con, query)
  
  return(df)
}


extract_bacteria_data_uniprot <- function(con, bact_codes) {
  
  bact_codes_string <- paste0(bact_codes, collapse = "', '")
  query <- sprintf("SELECT DISTINCT  
  uni.uniprot
  , org.organism
  FROM pockets poc
  LEFT JOIN pdb ON pdb.id = poc.id_pdb
  LEFT JOIN mapping_pdb_uniprot mpu ON mpu.id_pdb = pdb.id
  LEFT JOIN uniprot uni ON uni.id = mpu.id_uniprot
  LEFT JOIN organisms org ON org.id = uni.id_organism
  WHERE org.organism IN ('%s')", bact_codes_string)
  
  df <- dbGetQuery(con, query)
  
  return(df)
}

assign_bacteria_codes <- function(bact_name, single=TRUE) {
  if (single) {
    if (bact_name == "Staphylococcus aureus") {
      codes <- c("STAAU")
    } else if (bact_name == "Staphylococcus epidermidis") {
      codes <- c("STAEP")
    } else if (bact_name == "Enterococcus faecalis") {
      codes <- c("ENTFL")
    } else if (bact_name == "Escherichia coli") {
      codes <- c("ECOLX")
    } else if (bact_name == "Klebsiella pneumoniae") {
      codes <- c("KLEPN")
    } else if (bact_name == "Pseudomonas aeruginosa") {
      codes <- c("PSEAI")
    }
  } else {
    if (bact_name == "Staphylococcus aureus") {
      codes <- c("STAAU", "STAAC", "STAAE", "STAAN", "STAAM", "STAAR", "STAAS", 
                 "STAAW", "STAA1", "STAAB", "STAA3", "STAA8")
    } else if (bact_name == "Staphylococcus epidermidis") {
      codes <- c("STAEP", "STAES", "STAEQ")
    } else if (bact_name == "Enterococcus faecalis") {
      codes <- c("ENTFL", "ENTFA")
    } else if (bact_name == "Escherichia coli") {
      codes <- c("ECOLX", "ECOLI", "ECOL6", "ECOLC", "ECOBD", "ECOCB", "ECODH", 
                 "ECOH1")
    } else if (bact_name == "Klebsiella pneumoniae") {
      codes <- c("KLEPN")
    } else if (bact_name == "Pseudomonas aeruginosa") {
      codes <- c("PSEAI", "PSEA8", "PSEA7", "PSEAB")
    }
  }
  return(codes)
}



#
# Rows for the summary table concerning data from BioGPS database,
# with numerical values for different bacteria about: 
# pockets, PDB codes and Uniprot codes
# 

db_con <- DBI::dbConnect(RSQLite::SQLite(), dbfilename_biogps_info)

# 
# Pockets
# 
tablerow_nr_biogps_pockets <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    df_bact_pockets <- extract_biogps_pockets(con=db_con, organism_name=sel_bacteria)
    data.frame(bacteria=sel_bacteria, nr_biogps_pockets=nrow(df_bact_pockets))
  })) %>% 
  dplyr::arrange(bacteria)

# 
# PDB codes
# 
tablerow_nr_biogps_pdb <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    df_bact_pdb <- extract_biogps_pdb(con=db_con, organism_name=sel_bacteria)
    data.frame(bacteria=sel_bacteria, nr_biogps_pdb=nrow(df_bact_pdb))
  })) %>% 
  dplyr::arrange(bacteria)

# 
# Uniprot codes
# 
tablerow_nr_biogps_uniprot <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    df_bact_uniprot <- extract_biogps_uniprot(con=db_con, organism_name=sel_bacteria)
    data.frame(bacteria=sel_bacteria, nr_biogps_uniprot=nrow(df_bact_uniprot))
  })) %>% 
  dplyr::arrange(bacteria)

DBI::dbDisconnect(conn = db_con)


# 
# Row for the summary table, with numerical values for different bacteria about: 
# number of org_code (with at least one uniprot)
tablerow_nr_org_codes <- kegg_bacteria_unified_df %>% 
  dplyr::select(bacteria, bact_kegg) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(bacteria) %>% 
  dplyr::add_count() %>% 
  dplyr::rename(nr_kegg_org_code = n) %>% 
  dplyr::select(bacteria, nr_kegg_org_code) %>% 
  dplyr::distinct()


# 
# Row for the summary table, with numerical values for different bacteria about: 
# number of uniprot, considering all the org_code values
tablerow_nr_uniprot_allcodes <- kegg_bacteria_unified_df %>% 
  dplyr::select(bacteria, uniprot) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(bacteria) %>% 
  dplyr::add_count() %>% 
  dplyr::rename(nr_kegg_uniprot = n) %>% 
  dplyr::select(bacteria, nr_kegg_uniprot) %>% 
  dplyr::distinct()


# 
# Row for the summary table, with numerical values for different bacteria about: 
# number of uniprot, considering all the ATCC codes
tablerow_nr_uniprot_atcc <- atcc_bacteria_uniprot_df %>% 
  dplyr::select(bacteria, uniprot) %>%
  dplyr::distinct() %>% 
  dplyr::group_by(bacteria) %>% 
  dplyr::add_count() %>% 
  dplyr::rename(nr_atcc_uniprot = n) %>% 
  dplyr::select(bacteria, nr_atcc_uniprot) %>% 
  dplyr::distinct()



#
# Read data from the other database, where we stored data on pockets and scores
# 
db_con <- dbConnect(RSQLite::SQLite(), file.path(db_dir, "biogps.db"))


# 
# Data for uniprot
# 
tablerow_nr_bactuniprot_singlecode <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    sel_bact_codes <- assign_bacteria_codes(bact_name=sel_bacteria, single=TRUE)
    df_bact_data <- extract_bacteria_data_uniprot(con=db_con, bact_codes=sel_bact_codes)
    data.frame(bacteria=sel_bacteria, nr_uniprot_singlecode=nrow(df_bact_data))
  })) %>% 
  dplyr::arrange(bacteria)


tablerow_nr_bactuniprot_multiplecode <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    sel_bact_codes <- assign_bacteria_codes(bact_name=sel_bacteria, single=FALSE)
    df_bact_data <- extract_bacteria_data_uniprot(con=db_con, bact_codes=sel_bact_codes)
    data.frame(bacteria=sel_bacteria, nr_uniprot_multiplecode=nrow(df_bact_data))
  })) %>% 
  dplyr::arrange(bacteria)

dbDisconnect(db_con)


tablerow_coverage_singlecode <- tablerow_nr_uniprot_atcc %>% 
  dplyr::left_join(tablerow_nr_bactuniprot_singlecode, by="bacteria") %>% 
  dplyr::mutate(cov_singlecode = nr_uniprot_singlecode/nr_atcc_uniprot) %>% 
  dplyr::select(bacteria, cov_singlecode) %>% 
  dplyr::arrange(bacteria)

tablerow_coverage_multiplecode <- tablerow_nr_uniprot_atcc %>% 
  dplyr::left_join(tablerow_nr_bactuniprot_multiplecode, by="bacteria") %>% 
  dplyr::mutate(cov_multiplecode = nr_uniprot_multiplecode/nr_atcc_uniprot) %>% 
  dplyr::select(bacteria, cov_multiplecode) %>% 
  dplyr::arrange(bacteria)


# 
# Create the final table and save as .CSV file
# 

tablerow_nr_kegg_additional_identifiers <- as.numeric(t(tablerow_nr_org_codes)[2,]) - 1

finaltable_df <- rbind(
  t(tablerow_nr_biogps_pockets)[2,],
  t(tablerow_nr_biogps_pdb)[2,],
  t(tablerow_nr_biogps_uniprot)[2,],
  tablerow_nr_kegg_additional_identifiers,
  t(tablerow_nr_uniprot_allcodes)[2,],
  t(tablerow_nr_uniprot_atcc)[2,],
  t(tablerow_nr_bactuniprot_singlecode)[2,],
  t(tablerow_nr_bactuniprot_multiplecode)[2,],
  t(tablerow_coverage_singlecode)[2,],
  t(tablerow_coverage_multiplecode)[2,]
)
colnames(finaltable_df) <- t(tablerow_nr_org_codes)[1,]

rownames(finaltable_df) <- c(
  "nr_biogps_pockets",
  "nr_biogps_pdb",
  "nr_biogps_uniprot",
  "nr_kegg_additional_orgcode",
  "nr_uniprot_allcodes",
  "nr_uniprot_atcc",
  "nr_bactuniprot_singlecode",
  "nr_bactuniprot_multiplecode",
  "coverage_singlecode",
  "coverage_multiplecode"
)


write.csv2(x = finaltable_df, file = finaltable_filename)

