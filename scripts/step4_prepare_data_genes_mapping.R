#
# Script to prepare the PDB-Uniprot mapping
# 
# Procedure to read data from the file pdbtosp.txt from Uniprot web site
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pdbtosp.txt
# Only some methods are considered: XRAY, NMR, EM, NEUTRON, IR, FIBER, OTHER
# 

library(dplyr)
library(magrittr)

# 
# Functions to handle information from PDB and Uniprot
# 

# This function reads the information of the pdb - uniprot conversion
pdbtosp_read_content <- function(pdbtosp_filename) {
  
  skip_lines_content <- paste0(
    "^ +UniProt - Swiss-Prot"
    ,"|^ +SIB Swiss"
    ,"|^ +European Bioinformatics"
    ,"|^ +Protein Information"
    ,"|.*Number of PDB"
    ,"|.*Number of Swiss"
    ,"|^Description:"
    ,"|^Release:"
    ,"|^The PDB database"
    ,"|^Name"
    ,"|^https"
    ,"|^PDB"
    ,"|^code"
    ,"|^---"
    ,"|^ ---"
    ,"|^___"
    ,"|^ ___"
    ,"|^Distributed"
    ,"|^Copyrighted"
  )
  
  # The input file is first read by using data.table fread function.
  # Some lines at the end of the file are skipped. 
  # Also the first 24 lines are skipped
  # (but we can read the whole file, extract some metadata and skip some rows)
  tmp_df <- data.table::fread(
    file=pdbtosp_filename, header=FALSE, sep="\t", #skip=24,
    col.names=c("origdata"), strip.white=FALSE, blank.lines.skip=TRUE
  ) %>% 
    dplyr::filter(!grepl(skip_lines_content, origdata)) %>% 
    dplyr::mutate(id = row_number())
  
  new_df <- tmp_df %>% 
    dplyr::mutate(origdata = gsub(" +", " ", origdata)) %>% 
    dplyr::mutate(newid = ifelse(substr(origdata, 1, 1) == " ", id-1, id)) %>% 
    dplyr::group_by(newid) %>% 
    dplyr::summarise(origdata = stringr::str_c(origdata, collapse = ""))
  
  id <- 0
  while(nrow(new_df) < nrow(tmp_df)) {
    id <- id+1
    # This re-assignment gets into tmp_df the content previously extracted
    tmp_df <- new_df
    
    # The key step is the check for the first character of "origdata". 
    # It is assumed that each line starts with a no-space character (PDB entry)
    # All the rows starting with an empty character will be concatenated
    # to previous ones. When more lines are present, more iterations are needed
    new_df <- tmp_df %>% 
      dplyr::rename(id=newid) %>% 
      dplyr::mutate(origdata = gsub(" +", " ", origdata)) %>% 
      dplyr::mutate(newid = ifelse(substr(origdata, 1, 1) == " ", id-1, id)) %>% 
      dplyr::group_by(newid) %>% 
      dplyr::summarise(origdata = stringr::str_c(origdata, collapse = ""))
    cat("ITERATION: ", id, "\tRATIO RETAINED LINES", nrow(new_df)/nrow(tmp_df), "\n")
  }
  
  # Finally eliminates the spaces around comma(s)
  final_df <- new_df %>% 
    dplyr::mutate(origdata = gsub(" , |, | ,", ",", origdata)) %>% 
    dplyr::select(-newid)
  
  return(final_df)
  
}


pdbtosp_extract_content <- function(df, methods=c("XRAY", "NMR", "EM", 
                                                  "NEUTRON", "IR", 
                                                  "FIBER", "OTHER")) {
  
  em_nores_df <- NULL
  em_withres_df <- NULL
  nmr_nores_df <- NULL
  xray_nores_df <- NULL
  xray_withres_df <- NULL
  neutron_withres_df <- NULL
  ir_nores_df <- NULL
  fiber_nores_df <- NULL
  fiber_withres_df <- NULL
  other_nores_df <- NULL
  other_withres_df <- NULL
  
  if ("XRAY" %in% methods) {
    xray_withres_df <- df %>% 
      dplyr::filter(grepl("X-ray .*\\..* A", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "XRAY",
        exp_res = trimws(substr(origdata, 11, 17)),
        info_all = trimws(substr(origdata, 19, 10000))
      )
    
    xray_nores_df <- df %>% 
      dplyr::filter(grepl("X-ray -", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "XRAY",
        exp_res = NA,
        info_all = trimws(substr(origdata, 14, 10000))
      )
  }
  
  if ("NMR" %in% methods) {
    nmr_nores_df <- df %>% 
      dplyr::filter(grepl("NMR -", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "NMR",
        exp_res = NA,
        info_all = trimws(substr(origdata, 12, 10000))
      )
  }
  
  if ("EM" %in% methods) {
    em_withres_df <- df %>%
      dplyr::mutate(EM = substr(origdata, 6, 7)) %>% 
      dplyr::filter(EM == "EM") %>%
      dplyr::filter(!grepl("EM -", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "EM",
        exp_res = trimws(substr(origdata, 9, 13)),
        info_all = trimws(substr(origdata, 16, 10000))
      ) %>% 
      dplyr::select(-EM)
    
    em_nores_df <- df %>%
      dplyr::filter(grepl("EM -", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "EM",
        exp_res = NA,
        info_all = trimws(substr(origdata, 11, 10000))
      )
  }
  
  if ("NEUTRON" %in% methods) {
    neutron_withres_df <- df %>%
      dplyr::mutate(Neutron = substr(origdata, 6, 12)) %>% 
      dplyr::filter(Neutron == "Neutron") %>%
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "NEUTRON",
        exp_res = trimws(substr(origdata, 14, 19)),
        info_all = trimws(substr(origdata, 21, 10000))
      ) %>% 
      dplyr::select(-Neutron)
  }
  
  if ("IR" %in% methods) {
    ir_nores_df <- df %>%
      dplyr::filter(grepl("IR -", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "IR",
        exp_res = NA,
        info_all = trimws(substr(origdata, 11, 10000))
      )
  }
  
  if ("FIBER" %in% methods) {
    fiber_nores_df <- df %>% 
      dplyr::filter(grepl("Fiber -", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "FIBER",
        exp_res = NA,
        info_all = trimws(substr(origdata, 14, 10000))
      )
    
    fiber_withres_df <- df %>% 
      dplyr::filter(grepl("Fiber .*\\..* A", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "FIBER",
        exp_res = trimws(substr(origdata, 11, 17)),
        info_all = trimws(substr(origdata, 19, 10000))
      )
  }
  
  if ("OTHER" %in% methods) {
    other_nores_df <- df %>% 
      dplyr::filter(grepl("Other -", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "OTHER",
        exp_res = NA,
        info_all = trimws(substr(origdata, 14, 10000))
      )
    
    other_withres_df <- df %>% 
      dplyr::filter(grepl("Other .*\\..* A", origdata)) %>% 
      dplyr::mutate(
        pdbid = substr(origdata, 0, 4),
        exp_method = "OTHER",
        exp_res = trimws(substr(origdata, 11, 17)),
        info_all = trimws(substr(origdata, 19, 10000))
      )
  }
  
  results_df <- rbind(
    em_nores_df,
    em_withres_df,
    nmr_nores_df,
    xray_nores_df,
    xray_withres_df,
    neutron_withres_df,
    ir_nores_df,
    fiber_nores_df,
    fiber_withres_df,
    other_nores_df,
    other_withres_df
  ) %>% 
    dplyr::mutate(gene_uniprot = gsub(" ^| $", "", info_all)) %>% 
    dplyr::mutate(exp_res = gsub(" +A", "", exp_res)) %>% 
    dplyr::select(-info_all)
  
  return(results_df)
  
}


# repo_dir <- "C:/GITHUB_REPO/eoanalysis"
repo_dir <- getwd()

# Input
data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")

# output
gene_mapping <- file.path(mapp_dir, "gene_uniprot.RDS")
gene_mapping_downloaded <- file.path(mapp_dir, "pdbtosp.txt")

# URL to download updated data
gene_mapping_url <- "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pdbtosp.txt"

# Download in case file does not exist
if (!file.exists(gene_mapping_downloaded)) {
  download.file(gene_mapping_url, gene_mapping_downloaded, method = "curl")
}


cat("READ DATA FROM TXT FILE (the procedure takes some minutes):\n")
tmp_gene_uniprot_df <- pdbtosp_read_content(gene_mapping_downloaded)

selected_methods <- c("XRAY", "NMR", "EM", "NEUTRON", "IR", "FIBER", "OTHER")

pdb2uniprot_df <- pdbtosp_extract_content(df=tmp_gene_uniprot_df,
                                          methods=selected_methods)

gene_uniprot_pdb_df <- pdb2uniprot_df %>%
  tidyr::separate_longer_delim(gene_uniprot, delim = ",") %>%
  dplyr::mutate(info_split = gsub("\\(|\\)", "", gene_uniprot)) %>%
  tidyr::separate_wider_delim(cols=info_split, delim=" ",
                              names=c("gene", "uniprot")) %>%
  tidyr::separate_wider_delim(cols=gene, delim="_",
                              names=c("genesymbol", "organism"),cols_remove = FALSE) %>%
  dplyr::select(-gene_uniprot) %>%
  dplyr::rename(pdb=pdbid)

gene_uniprot_df <- gene_uniprot_pdb_df %>%
  dplyr::select(genesymbol, uniprot, organism) %>%
  dplyr::distinct()


saveRDS(object=gene_uniprot_df, file=gene_mapping)
