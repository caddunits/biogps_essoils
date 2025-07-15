#
# Script to prepare the KEGG mapping
# Data is read from original .RDS files and saved as unique .RDS file
#
# Note: single RDS files are created with a procedure that downloads data 
# from the KEGG website (pathway-uniprot connections) by using the script 
# step2a_read_kegg_data.R
#
# Note: needs internet connection

# 
# In the following KEGG website
#   https://www.kegg.jp/brite/br08611
# we identified the codes of interest for each bacteria.
# For example, by searching "Staphylococcus epidermidis" we have:
# ser: RP62A (MRSE) 
# sep: ATCC 12228
# sepp: PM221 [only NCBI – CDM…, no Uniprot data]
# seps: SEI
# After manual check of all the options, we decided which codes to use
# 

library(dplyr)
library(magrittr)

# BiocManager::install("UniProt.ws")
library(UniProt.ws) 


repo_dir <- getwd()

# Input
data_dir <- file.path(repo_dir, "data")
kegg_dir <- file.path(data_dir, "KEGG_data")

# Output
mapp_dir <- file.path(data_dir, "Mapping_data")
kegg_mapping <- file.path(mapp_dir, "kegg_bacteria_uniprot.RDS")
gene_mapping_kegg <- file.path(mapp_dir, "gene_uniprot_kegg.RDS")


# 
# NOVITA: VEDERE SE PUO SOSTITUIRE GLI ALTRI DUE 
# 
kegg_bacteria_unified <- file.path(mapp_dir, "kegg_bacteria_unified.RDS")


# 'Pseudomonas aeruginosa': 21 files
kegg_paeru <- c(
  'pae', 'paev', 'paei', 'pau', 'pap', 'pag', 'paf', 'pnc', 'paeb', 'pdk',
  'psg', 'prp', 'paep', 'paer', 'paem', 'pael', 'paes', 'paeu', 'paeg', 'paec',
  'sech'
)

# 'Klebsiella pneumoniae': 24 files
kegg_kpneu <- c(
  'kpn', 'kpu', 'kpm', 'kpp', 'kph', 'kpz', 'kpv', 'kpw', 'kpy', 'kpg',
  'kpc', 'kpq', 'kpt', 'kpo', 'kpr', 'kpj', 'kpi', 'kpa', 'kps', 'kpx',
  'kpb', 'kpne', 'kpnu', 'kpnk'
)

# 'Enterococcus faecalis': 8 files
kegg_efaec <- c(
  'efa', 'efl', 'efi', 'efd', 'efs', 'efn', 'efq', 'ene'
)

# 'Staphylococcus epidermidis': 4 files
kegg_aepid <- c(
  'ser', 'sep', 'sepp', 'seps'
)

# 'Staphylococcus aureus': 53 files
kegg_saure <- c(
  'sac', 'saa', 'sab', 'sad', 'sae', 'sah', 'saj', 'sam', 'sams', 'sao',
  'sar', 'sas', 'sau', 'saua', 'saub', 'sauc', 'saud', 'saue', 'sauf', 'saug', 
  'saui', 'sauj', 'sauk', 'saum', 'saun', 'sauq', 'saur', 'saus', 'saut', 'sauu', 
  'sauv', 'sauw', 'saux', 'sauy', 'sauz', 'sav', 'saw', 'sax', 'suc', 'sud', 
  'sue', 'suf', 'sug', 'suj', 'suk', 'suq', 'sut', 'suu', 'suv', 'suw', 
  'sux', 'suy', 'suz'
)

# 'Escherichia coli': 67 files
kegg_ecoli <- c(
  'eco', 'eab', 'ebd', 'ebe', 'ebl', 'ebr', 'ebw', 'ecc', 'ecd', 'ece', 
  'ecf', 'ecg', 'eci', 'ecj', 'eck', 'ecl', 'ecm', 'ecoa', 'ecob', 'ecoc', 
  'ecoh', 'ecoi', 'ecoj', 'ecok', 'ecol', 'ecoo', 'ecos', 'ecp', 'ecq', 'ecr', 
  'ecs', 'ect', 'ecv', 'ecw', 'ecx', 'ecy', 'ecz', 'edh', 'edj', 'eih', 
  'ekf', 'eko', 'elc', 'eld', 'elf', 'elh', 'ell', 'eln', 'elo', 'elp',
  'elr', 'elu', 'elw', 'elx', 'ena', 'eoc', 'eoh', 'eoi', 'eoj', 'eok',
  'ese', 'esl', 'esm', 'eso', 'etw', 'eum', 'eun'
)



kegg_correspondences <- data.frame(
  bact_kegg = c(
    kegg_paeru,
    kegg_kpneu,
    kegg_efaec,
    kegg_aepid,
    kegg_saure,
    kegg_ecoli
  ),
  bacteria = c(
    rep('Pseudomonas aeruginosa', length(kegg_paeru)),
    rep('Klebsiella pneumoniae', length(kegg_kpneu)),
    rep('Enterococcus faecalis', length(kegg_efaec)),
    rep('Staphylococcus epidermidis', length(kegg_aepid)),
    rep('Staphylococcus aureus', length(kegg_saure)),
    rep('Escherichia coli', length(kegg_ecoli))
  ),
  stringsAsFactors = FALSE
)

kegg_files <- list.files(kegg_dir, pattern=".RDS$", full.names=T, recursive=T)
cat("READ DATA FROM THE FOLLOWING RDS FILES:\n")
print(kegg_files)

kegg_bacteria_uniprot_df <- do.call('rbind', lapply(kegg_files, function(f) {
  kegg_code <- gsub(pattern='.RDS', replacement='', x=basename(f))
  tmp_kegg_df <- readRDS(file=f) %>% 
    dplyr::mutate(bact_kegg=kegg_code)
  if ('NCBI' %in% colnames(tmp_kegg_df)) tmp_kegg_df %<>% dplyr::select(-NCBI)
  tmp_kegg_df
})) %>% 
  dplyr::filter(pw_name!="", uniprot!="") %>%
  dplyr::rename(pathway_code=pw_name, pathway_name=pw_value) %>%
  dplyr::select(bact_kegg, uniprot, pathway_code, pathway_name) %>% 
  dplyr::mutate(uniprot = strsplit(as.character(uniprot), " ")) %>% 
  tidyr::unnest(uniprot) %>% 
  dplyr::left_join(kegg_correspondences, by="bact_kegg")

saveRDS(object=kegg_bacteria_uniprot_df, file=kegg_mapping)


# First, extract the list of Uniprot codes
kegg_uniprot_lst <- kegg_bacteria_uniprot_df %>% 
  dplyr::select(uniprot) %>% 
  dplyr::distinct() %>% 
  dplyr::arrange() %>% 
  unlist() %>% 
  as.character()


# Then, use a specific function from UniProt.ws 
# UniProtKB_AC-ID (default) is short for "UniProt accession identifiers"
kegg_conversion_df <- UniProt.ws::mapUniProt(
  from="UniProtKB_AC-ID",
  to='UniProtKB',
  query=kegg_uniprot_lst
)


# check the success rate
kegg_success_rate <- 100 * sum(kegg_uniprot_lst %in% kegg_conversion_df$From) / length(kegg_uniprot_lst)

cat("CONVERTED", kegg_success_rate, "% OF UNIPROT CODES\n")

saveRDS(object=kegg_conversion_df, file=gene_mapping_kegg)



# 
# ATTENZIONE, VALUTARE SE AGGIUNGERE QUESTA QUA INVECE CHE NELLA SCRIPT SUCCESSIVA,
# E SCRIVERE QUINDI UN SOLO FILE .RDS
# (SERVONO I DUE DF SEPARATI? kegg_bacteria_uniprot_df e kegg_conversion_df)
# 

# Join with data from conversion uniprot-gene
kegg_bacteria_unified_df <- kegg_bacteria_uniprot_df %>% 
  dplyr::left_join(kegg_conversion_df, by=c('uniprot'='From')) %>% 
  tidyr::separate_wider_delim(cols='Entry.Name', delim="_", 
                              names=c("gene", "organism")) %>% 
  dplyr::rename(org_details = Organism, 
                protein_names = Protein.names,
                gene_names = Gene.Names) %>% 
  dplyr::select(-Length, -Entry)

saveRDS(object=kegg_bacteria_unified_df, file=kegg_bacteria_unified)



# NOTE SPIEGAZIONE:
# Il df kegg_bacteria_uniprot_df contiene le associazioni uniprot-pathway per 
# tutti i batteri
# Il df kegg_conversion_df contiene i dettagli di tutti i geni per tutti i batteri
# 
# Il df kegg_bacteria_unified_df è un join delle due tabelle
# 
# 
# DA CHIARIRE: 
# * SERVONO TUTTE LE COLONNE? 
# * I DUE DF SEPARATI VENGONO MAI USATI SINGOLARMENTE?
# 



# NOTA: NCBI codes, extracted from KEGG, are not mapped
# https://www.ncbi.nlm.nih.gov/ 
# xxxx <- UniProt.ws::mapUniProt(
#   from="UniProtKB_AC-ID",
#   to='UniProtKB',
#   query=c('CAQ30519')
# )

