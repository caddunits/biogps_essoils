#
# Script to prepare the data from the BioGPS calculations
# 
# Procedure to read data from the BioGPS database (sqlite) where we previously
# stored data by means of other scripts
# 

library(dplyr)
library(magrittr)
library(DBI)


# 
# Functions
# 

extract_uniprot_biogps_data <- function(con) {
  query <- "SELECT DISTINCT zz.id_poc
, poc.name
, pdb.pdb
, mpu.id_uniprot
, uni.uniprot
, org.organism
FROM zzscores_molpoc zz
LEFT JOIN pockets poc ON poc.id = zz.id_poc
LEFT JOIN pdb ON pdb.id = poc.id_pdb
LEFT JOIN mapping_pdb_uniprot mpu ON mpu.id_pdb = poc.id_pdb
LEFT JOIN uniprot uni ON uni.id = mpu.id_uniprot
LEFT JOIN organisms org ON org.id = uni.id_organism
WHERE id_uniprot IS NOT NULL"
  
  df <- dbGetQuery(con, query)
  
  return(df)
  
}


extract_pockets_stat_data <- function(con, norm_set) {
  
  if (!norm_set %in% c('drugcentral')) {
    return(NULL)
  }
  
  query <- sprintf("SELECT * 
    FROM stats_poc WHERE poc_refcode == '%s'", norm_set)
  
  df <- dbGetQuery(con, query)
  
  return(df)
}


extract_zscorepoc_data <- function(con, norm_set) {
  
  if (!norm_set %in% c('drugcentral')) {
    return(NULL)
  }
  
  query_zscores_poc <- sprintf("SELECT DISTINCT *
FROM zscores_poc zp
WHERE zp.poc_refcode == '%s'", norm_set)
  
  df <- dbGetQuery(con, query_zscores_poc)
  return(df)
}


extract_zscoremol_data <- function(con) {
  
  query_zscores_mol <- "SELECT DISTINCT * FROM zscores_mol zp"
  
  df <- dbGetQuery(con, query_zscores_mol)
  return(df)
}


extract_zzscore_data <- function(con, norm_set) {
  
  if (!norm_set %in% c('drugcentral')) {
    return(NULL)
  }
  
  query_zzscores <- sprintf("SELECT zz.id_poc
, poc.name AS pocketname
, poc.id_pdb
, pdb.pdb
, uni.uniprot
, uni.gene
, gsm.genesymbol
, org.organism
, mol.id AS id_mol
, mol.name AS molname
, mol.info
, zz.zzscore_molpoc
, zz.poc_refcode
FROM zzscores_molpoc zz
LEFT JOIN molecules mol ON mol.id = zz.id_mol
LEFT JOIN pockets poc ON poc.id = zz.id_poc
LEFT JOIN pdb ON pdb.id = poc.id_pdb
LEFT JOIN mapping_pdb_uniprot mpu ON mpu.id_pdb = pdb.id
LEFT JOIN uniprot uni ON uni.id = mpu.id_uniprot
LEFT JOIN genesymbols gsm ON gsm.id = uni.id_genesymbol
LEFT JOIN organisms org ON org.id = uni.id_organism
WHERE zz.poc_refcode == '%s'
AND uni.uniprot IS NOT NULL
", norm_set)
  
  df_init <- dbGetQuery(con, query_zzscores)
  return(df_init)
}


refine_zzscore_data <- function(df) {
  
  # Refine the data by taking only the first pocket for each gene
  df_refined <- df %>%
    dplyr::group_by(molname, gene) %>%
    dplyr::arrange(molname, desc(zzscore_molpoc)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Here the product by -1 is used to sort values, then a rank column is added
  df_refined2 <- df_refined %>%
    dplyr::mutate(tmpvar=-1*zzscore_molpoc) %>%
    dplyr::group_by(molname) %>%
    mutate(zzrank=dense_rank(tmpvar)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-tmpvar) %>%
    dplyr::arrange(molname, zzrank) %>%
    dplyr::filter(zzscore_molpoc > 0) %>%
    dplyr::mutate(zzscore = round(x=zzscore_molpoc, digits = 4)) %>%
    dplyr::select(molname, pocketname,
                  pdb, uniprot, genesymbol, organism,
                  zzscore, zzrank, info, poc_refcode,
                  id_mol, id_poc, nr_pockets)
  
  return(df_refined2)
}

# 
# load all the project details
# 

# 
# Targets, Classes, Bacteria
# 
norm_set <- "drugcentral"


# Change this line according to the path of your repository
repo_dir <- getwd()

# Input
data_dir <- file.path(repo_dir, "data")
db_dir <- file.path(repo_dir, "db")
# mapp_dir <- file.path(data_dir, "Mapping_data")


# Output
biogps_data <- file.path(data_dir, "Biogps_zscore_data", "biogps_data.RDS")


# bacteriadata_filename <- file.path(data_dir, "BACTERIADATA.xlsx")
# 
# # Input info on bacteria:
# # Codes on different databases (ATCC, KEGG), from xlsx
# df_bacteria <- readxl::read_excel(path=bacteriadata_filename) %>% 
#   dplyr::filter(bacteria %in% wanted_bacteria)



# 
# Reading ZZ-score data from database
# Read from database all ZZ-score values
# 

# First, extract from the database the ZZ-score data
db_con <- dbConnect(RSQLite::SQLite(), file.path(db_dir, "biogps.db"))
df_zzscore <- extract_zzscore_data(con=db_con, norm_set="drugcentral")


# Then, count how many pockets per gene 
df_gene_pocket_stats <- df_zzscore %>% 
  dplyr::select(pocketname, genesymbol) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(genesymbol) %>% 
  dplyr::add_count() %>% 
  dplyr::rename(nr_pockets = n) %>% 
  dplyr::distinct()


# And finally add the information on the number of pockets for each gene. In
# addition, create a new variable, named zzscore                                (ATTENZIONE: E' NECESSARIO QUESTO O SI POTEVA USARE zzscore_molpoc?)
df_full <- df_zzscore %>% 
  dplyr::mutate(zzscore = round(x=zzscore_molpoc, digits = 4)) %>% 
  dplyr::left_join(df_gene_pocket_stats, by=c('pocketname', 'genesymbol'))

zscore_mol_df <- extract_zscoremol_data(con=db_con)
zscore_poc_df <- extract_zscorepoc_data(con=db_con, norm_set=norm_set)

zscores <- list(
  full=df_full,
  mol=zscore_mol_df,
  poc=zscore_poc_df
)

saveRDS(object=zscores, file=biogps_data)


# 
# ATTENZIONE: SEMBRA CHE QUESTA SI USI SOLO PER ALLUVIAL PLOT, 
# QUINDI NON LO METTIAMO QUA
# 
# Refined data: with the function refine_zzscore_data we select the best value 
# for each gene (we keep only the best pocket-molecule pair) 
# ATTENTION: It uses only zzscore (does not check zscore_mol and zscore_poc)
# df_refined <- df_full %>% 
#   refine_zzscore_data() %>% 
#   dplyr::arrange(dplyr::desc(zzscore)) %>% 
#   dplyr::left_join(df_bacteria, by="organism") %>% 
#   dplyr::left_join(df_composition, by="molname")

dbDisconnect(db_con)
