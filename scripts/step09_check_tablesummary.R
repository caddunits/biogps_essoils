#
# Script to prepare a table with the summary of bacterial data, with information
# from KEGG, ATCC and BioGPS.
# Data is read from original .RDS files and saved as unique .CSV file
# which corresponds to the table for the paper
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
dbfilename_biogps_info <- "C:/BIOGPS_INFO/info.sql"

wanted_bacteria <- c(
  "Staphylococcus aureus",
  "Staphylococcus epidermidis",
  "Enterococcus faecalis",
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Pseudomonas aeruginosa"
)


# Output
finaltable_filename <- file.path(repo_dir, "output", "finaltable.csv")

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

print(tablerow_nr_org_codes)

# 
# Manual curation for some rows with the column "bacteria" non specified
# (but clearly assignable)
# 
# kegg_bacteria_unified_df %>% 
#   dplyr::filter(is.na(bacteria)) %>% 
#   dplyr::select(organism) %>% 
#   dplyr::distinct()
# 
# # They are all STAAU 
# # => We can assign in the column bacteria the value
# kegg_bacteria_unified_df$bacteria[is.na(kegg_bacteria_unified_df$bacteria)] <- "Staphylococcus aureus"

# tablerow_nr_org_codes <- kegg_bacteria_unified_df %>% 
#   dplyr::select(bacteria, bact_kegg) %>% 
#   dplyr::distinct() %>% 
#   dplyr::group_by(bacteria) %>% 
#   dplyr::add_count() %>% 
#   dplyr::rename(nr_kegg_org_code = n) %>% 
#   dplyr::select(bacteria, nr_kegg_org_code) %>% 
#   dplyr::distinct()


# 
# Row for the summary table, with numerical values for different bacteria about: 
# number of uniprot, considering only the main org_code 
tablerow_nr_uniprot_singlecodes <- kegg_bacteria_unified_df %>% 
  dplyr::filter(bact_kegg %in% c('sac', 'ser', 'efq', 'eco', 'kpn', 'pae')) %>% 
  dplyr::select(bacteria, uniprot) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(bacteria) %>% 
  dplyr::add_count() %>% 
  dplyr::rename(nr_kegg_uniprot = n) %>% 
  dplyr::select(bacteria, nr_kegg_uniprot) %>% 
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
# number of genes, considering all the org_code values
# tablerow_nr_genes_allcodes <- kegg_bacteria_unified_df %>% 
#   dplyr::select(bacteria, gene) %>% 
#   dplyr::distinct() %>% 
#   dplyr::group_by(bacteria) %>% 
#   dplyr::add_count() %>% 
#   dplyr::rename(nr_kegg_genes = n) %>% 
#   dplyr::select(bacteria, nr_kegg_genes) %>% 
#   dplyr::distinct()


# 
# Row for the summary table, with numerical values for different bacteria about: 
# number of genes, considering all the ATCC codes
tablerow_nr_genes_atcc <- atcc_bacteria_uniprot_df %>% 
  dplyr::select(bacteria, genename) %>%
  dplyr::distinct() %>% 
  dplyr::group_by(bacteria) %>% 
  dplyr::add_count() %>% 
  dplyr::rename(nr_atcc_genes = n) %>% 
  dplyr::select(bacteria, nr_atcc_genes) %>% 
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
    }))

# 
# PDB codes
# 
tablerow_nr_biogps_pdb <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    df_bact_pdb <- extract_biogps_pdb(con=db_con, organism_name=sel_bacteria)
    data.frame(bacteria=sel_bacteria, nr_biogps_pdb=nrow(df_bact_pdb))
  }))

# 
# Uniprot codes
# 
tablerow_nr_biogps_uniprot <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    df_bact_uniprot <- extract_biogps_uniprot(con=db_con, organism_name=sel_bacteria)
    data.frame(bacteria=sel_bacteria, nr_biogps_uniprot=nrow(df_bact_uniprot))
  }))

DBI::dbDisconnect(conn = db_con)


#
# Read data from the other database, where we stored data on pockets and scores
# 
db_dir <- file.path(repo_dir, "db")
db_con <- dbConnect(RSQLite::SQLite(), file.path(db_dir, "biogps.db"))

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
# Data for pockets
# 
# tablerow_nr_bactpockets_singlecode <- do.call(
#   'rbind', 
#   lapply(wanted_bacteria, function(sel_bacteria) { 
#     sel_bact_codes <- assign_bacteria_codes(bact_name=sel_bacteria, single=TRUE)
#     df_bact_data <- extract_bacteria_data_pockets(con=db_con, bact_codes=sel_bact_codes)
#     data.frame(bacteria=sel_bacteria, nr_pockets_singlecode=nrow(df_bact_data))
#   }))
# 
# tablerow_nr_bactpockets_multiplecode <- do.call(
#   'rbind', 
#   lapply(wanted_bacteria, function(sel_bacteria) { 
#     sel_bact_codes <- assign_bacteria_codes(bact_name=sel_bacteria, single=FALSE)
#     df_bact_data <- extract_bacteria_data_pockets(con=db_con, bact_codes=sel_bact_codes)
#     data.frame(bacteria=sel_bacteria, nr_pockets_multiplecode=nrow(df_bact_data))
#   }))


# 
# Data for pdb
# 
# tablerow_nr_bactpdb_singlecode <- do.call(
#   'rbind', 
#   lapply(wanted_bacteria, function(sel_bacteria) { 
#     sel_bact_codes <- assign_bacteria_codes(bact_name=sel_bacteria, single=TRUE)
#     df_bact_data <- extract_bacteria_data_pdb(con=db_con, bact_codes=sel_bact_codes)
#     data.frame(bacteria=sel_bacteria, nr_pdb_singlecode=nrow(df_bact_data))
#   }))
# 
# tablerow_nr_bactpdb_multiplecode <- do.call(
#   'rbind', 
#   lapply(wanted_bacteria, function(sel_bacteria) { 
#     sel_bact_codes <- assign_bacteria_codes(bact_name=sel_bacteria, single=FALSE)
#     df_bact_data <- extract_bacteria_data_pdb(con=db_con, bact_codes=sel_bact_codes)
#     data.frame(bacteria=sel_bacteria, nr_pdb_multiplecode=nrow(df_bact_data))
#   }))

# 
# Data for uniprot
# 
tablerow_nr_bactuniprot_singlecode <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    sel_bact_codes <- assign_bacteria_codes(bact_name=sel_bacteria, single=TRUE)
    df_bact_data <- extract_bacteria_data_uniprot(con=db_con, bact_codes=sel_bact_codes)
    data.frame(bacteria=sel_bacteria, nr_uniprot_singlecode=nrow(df_bact_data))
  }))

tablerow_nr_bactuniprot_multiplecode <- do.call(
  'rbind', 
  lapply(wanted_bacteria, function(sel_bacteria) { 
    sel_bact_codes <- assign_bacteria_codes(bact_name=sel_bacteria, single=FALSE)
    df_bact_data <- extract_bacteria_data_uniprot(con=db_con, bact_codes=sel_bact_codes)
    data.frame(bacteria=sel_bacteria, nr_uniprot_multiplecode=nrow(df_bact_data))
  }))



dbDisconnect(db_con)


tablerow_coverage_singlecode <- tablerow_nr_uniprot_atcc %>% 
  dplyr::left_join(tablerow_nr_bactuniprot_singlecode, by="bacteria") %>% 
  dplyr::mutate(cov_singlecode = nr_uniprot_singlecode/nr_atcc_uniprot) %>% 
  dplyr::select(bacteria, cov_singlecode)

tablerow_coverage_multiplecode <- tablerow_nr_uniprot_atcc %>% 
  dplyr::left_join(tablerow_nr_bactuniprot_multiplecode, by="bacteria") %>% 
  dplyr::mutate(cov_multiplecode = nr_uniprot_multiplecode/nr_atcc_uniprot) %>% 
  dplyr::select(bacteria, cov_multiplecode)


# 
# Create the final table and save as .CSV file
finaltable_df <- rbind(
  t(tablerow_nr_uniprot_singlecodes)[2,],
  t(tablerow_nr_uniprot_allcodes)[2,],
  # t(tablerow_nr_genes_allcodes)[2,],
  t(tablerow_nr_uniprot_atcc)[2,],
  t(tablerow_nr_genes_atcc)[2,],
  t(tablerow_nr_biogps_pockets)[2,],
  t(tablerow_nr_biogps_pdb)[2,],
  t(tablerow_nr_biogps_uniprot)[2,],
  # t(tablerow_nr_bactpockets_singlecode)[2,],
  # t(tablerow_nr_bactpockets_multiplecode)[2,],
  # t(tablerow_nr_bactpdb_singlecode)[2,],
  # t(tablerow_nr_bactpdb_multiplecode)[2,],
  t(tablerow_nr_bactuniprot_singlecode)[2,],
  t(tablerow_nr_bactuniprot_multiplecode)[2,],
  t(tablerow_coverage_singlecode)[2,],
  t(tablerow_coverage_multiplecode)[2,]
) %>% as.data.frame()
colnames(finaltable_df) <- t(tablerow_nr_uniprot_singlecodes)[1,]
rownames(finaltable_df) <- c(
  "nr_kegg_uniprot_orgcode",
  "nr_kegg_uniprot_allorgcodes",
  "nr_kegg_genes_allorgcodes",
  "nr_atcc_uniprot",
  "nr_atcc_genes",
  "nr_biogps_pockets",
  "nr_biogps_pdb",
  "nr_biogps_uniprot",
  # "nr_bactpockets_singlecode",
  # "nr_bactpockets_multiplecode",
  # "nr_bactpdb_singlecode",
  # "nr_bactpdb_multiplecode",
  "nr_bactuniprot_singlecode",
  "nr_bactuniprot_multiplecode",
  "coverage_singlecode",
  "coverage_multiplecode"
)


write.csv2(x = finaltable_df, file = finaltable_filename)

