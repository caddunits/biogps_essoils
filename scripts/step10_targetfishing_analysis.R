# 
# This script run the first part of the analysis: target fishing
# 
# Input are data from .RDS files
# Output are tabular data (as .XLSX)
# 

library(dplyr)
library(magrittr)
library(circlize)
library(tibble)
library(purrr)
library(writexl)
library(ggpubr)
library(ggplot2)
# BiocManager::install("STRINGdb")
library(STRINGdb)
library(igraph)


# 
# Functions
# 

# Normalize data about node centrality, as in cytoscape
normalize_as_cytoscape <- function(x) { (x - min(x)) / (max(x) - min(x)) }

# Extract percent value from a corresponding array of values, 
# compared to a given threshold
perc_rank <- function(x, x_perc) { sum(x_perc < x) }



# Change this line according to the path of your repository
repo_dir <- getwd()

# Settings are set once, in a configuration file
source(file.path(repo_dir, "scripts", "config.R"))


# Input
data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")

# Output
output_dir <- file.path(repo_dir, "output", "target_fishing")

# Create the ouput directory if it does not exist
dir.create(output_dir, showWarnings = FALSE)




# 
# First, read from .RDS files the data we need in the various steps
# 

biogps_data <- readRDS(file = file.path(data_dir, "Biogps_zscore_data", "biogps_data.RDS"))
df_full <- biogps_data$full
zscore_mol_df <- biogps_data$mol
zscore_poc_df <- biogps_data$poc

df_bacteria <- readRDS(file = file.path(mapp_dir, "bacteria_data.RDS"))

df_composition_new <- readRDS(file = file.path(mapp_dir, "composition_new.RDS"))



string_db_list <- list()
for (bact in config$wanted_bacteria) {
  species_code <- config$species_mapping[[bact]]
  cat("BACTERIUM:", bact, "CODE:", species_code, "\n")
  # cat("BACTERIUM:", bact, "CODE:", species_mapping[bact], "\n")
  string_db_list[bact] <- STRINGdb$new(
    version = "12.0",
    species = species_code,
    # species = species_mapping[bact],
    score_threshold = 400,
    network_type = "full"
  )
}


# Targets for each oil (considered as phytocomplex)
# Note: we set relationship = "many-to-many" to avoid warning
# In case of multiple data for the pair id_mol,id_poc we get the best value
df_targets_full <- df_full %>% 
  dplyr::filter(zzscore > config$threshold_zzscore) %>%
  dplyr::left_join(zscore_mol_df, by=c('id_mol', 'id_poc'),
                   relationship = "many-to-many") %>%
  dplyr::filter(zscore_mol > config$threshold_zscore_mol) %>%
  dplyr::left_join(zscore_poc_df, by=c('id_mol', 'id_poc', 'poc_refcode'),
                   relationship = "many-to-many") %>%
  dplyr::filter(zscore_poc > config$threshold_zscore_poc) %>% 
  dplyr::group_by(molname, gene) %>%
  dplyr::arrange(molname, desc(zzscore_molpoc)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()


# We add information about bacteria
df_targets_bact <- df_targets_full %>% 
  dplyr::left_join(df_bacteria, by="organism")


# We add the data for composition, and filter out data for molecules in traces
# We allow many-to-many relationships because the same molecule can be present 
# in more than one oil, and can have more than one target
df_targets_final <- df_targets_bact %>% 
  dplyr::left_join(df_composition_new, by="molname", 
                   relationship = "many-to-many") %>% 
  dplyr::filter(!is.na(bacteria)) %>% 
  dplyr::filter(perc_comp > config$threshold_comp)


# 
# Write data into an excel file
# 
wanted_columns <- c('oil_eng', 'fullname', 'phytoclass', 'perc_comp',
                    'zzscore_molpoc', 'zscore_mol', 'zscore_poc',
                    'organism', 'bacteria', 'bact_details',
                    'genesymbol', 'uniprot', 'pdb', 'pocketname', 'nr_pockets')
df_targets_filtered <- df_targets_final %>% dplyr::select(all_of(wanted_columns))



#
# If in the analysis we set thresholds to 2, 1, 1 we have very few targets
# Instead, if thresholds are set to 0, 0, 0
# Debug: 16806 rows (16805 records)
# FullData: 173957 rows
# Oregano: 5392 rows
# Thyme: 8481 rows
# Cinnamon: 2935 rows
# 


# 
# Lists of targets for phytocomplexes
# 

df_targets_contributions <- df_targets_final %>%
  dplyr::select(bacteria, oil_eng, fullname, perc_comp, genesymbol, zzscore) %>%
  dplyr::mutate(perc_weight = ifelse(perc_comp>1, 1+log10(perc_comp), perc_comp)) %>%
  dplyr::mutate(contribution = zzscore * perc_weight)

df_targets_phytocomplex <- df_targets_contributions %>%
  dplyr::group_by(genesymbol, oil_eng, bacteria) %>%
  dplyr::summarise(sum_contr = sum(contribution)) %>% 
  dplyr::filter(!is.na(bacteria))

oils <- unique(df_targets_phytocomplex$oil_eng)
bacteria <- unique(df_targets_phytocomplex$bacteria)

# Create a dataframe where to add data along the loop
df_combinations <- tidyr::crossing(
  bacteria = bacteria,
  oil_eng = oils
) %>% 
  dplyr::mutate(
    nr_targets=NA,
    targets_found=NA,
    perc_targets_found=NA
  )

network_targets_list <- list()
# network_details_list <- list()  # ATTENZIONE PROBABILMENTE VA TOLTA
names_to_export <- c()
idx <- 0 # This index is used to fill 

# Loop over the bacteria to extract data for the centrality of targets
for (bact in bacteria) {
  
  # For the given bacterium we extract the data from stringdb
  string_db <- string_db_list[[bact]]
  
  graph <- string_db$get_graph()
  
  # Extract data for the whole set of targets of the givwen bacterium
  # (all the nodes of the network)
  bact_degree <- igraph::degree(graph, V(graph)$name)
  bact_betweenness <- igraph::betweenness(graph, V(graph)$name)
  bact_closeness <- igraph::closeness(graph, V(graph)$name)
  
  # We calculate the quantile
  bact_degree_percent <- quantile(bact_degree, probs=seq(0, 1, 0.01))
  bact_betweenness_percent <- quantile(bact_betweenness, probs=seq(0, 1, 0.01))
  bact_closeness_percent <- quantile(bact_closeness, probs=seq(0, 1, 0.01))
  
  
  # loop over the oils
  for (oil in oils) {
    idx <- idx+1
    
    cat("\nWRITING TARGETS FOR BACTERIA", bact, "AND OIL", oil, "\n")
    sel_targets_df <- df_targets_phytocomplex %>% 
      dplyr::filter(bacteria==bact, oil_eng==oil) %>% 
      dplyr::arrange(desc(sum_contr))
    
    cat("Nr TARGETS:", nrow(sel_targets_df), "\n")
    df_combinations$nr_targets[df_combinations$oil_eng==oil & 
                               df_combinations$bacteria==bact] <- nrow(sel_targets_df)
    
    # 
    # Bioinformatics analysis for centrality
    # 
    df_mapping_orig <- string_db$map(
      as.data.frame(sel_targets_df),
      "genesymbol",
      removeUnmappedRows = TRUE
    ) %>% 
      dplyr::distinct()
    
    cat("\nCheck data:")
    check_names <- sum(!df_mapping_orig$STRING_id %in% V(graph)$name)
    # print(check_names)
    
    # We might have some nodes that are not present in the graph
    if (check_names > 0) {
      error_idx <- which(!df_mapping_orig$STRING_id %in% V(graph)$name)
      print(df_mapping_orig$STRING_id[error_idx])
      nodes_to_map <- df_mapping_orig$STRING_id[-error_idx]
      
    } else {
      nodes_to_map <- df_mapping_orig$STRING_id
      print("Everything is fine")
    }
    

    # bact_degree[nodes_to_map]
    df_degree <- data.frame(
      STRING_id = nodes_to_map,
      degree = as.numeric(bact_degree[nodes_to_map])
    ) %>% 
      dplyr::distinct()
    
    df_betweenness <- data.frame(
      STRING_id = nodes_to_map,
      betweenness = as.numeric(bact_betweenness[nodes_to_map])
    ) %>% 
      dplyr::distinct()
    
    df_closeness <- data.frame(
      STRING_id = nodes_to_map,
      closeness = as.numeric(bact_closeness[nodes_to_map])
    ) %>% 
      dplyr::distinct()
    
    df_perc_betweenness <- bact_betweenness[nodes_to_map] %>% 
      purrr::map_dbl(perc_rank, bact_betweenness_percent) %>% 
      tibble::enframe() %>% 
      dplyr::rename(STRING_id=name, betweenness_percent=value) %>% 
      dplyr::distinct()
    
    df_perc_degree <- bact_degree[nodes_to_map] %>% 
      purrr::map_dbl(perc_rank, bact_degree_percent) %>% 
      tibble::enframe() %>% 
      dplyr::rename(STRING_id=name, degree_percent=value) %>% 
      dplyr::distinct()
    
    df_perc_closeness <- bact_closeness[nodes_to_map] %>% 
      purrr::map_dbl(perc_rank, bact_closeness_percent) %>% 
      tibble::enframe() %>% 
      dplyr::rename(STRING_id=name, closeness_percent=value) %>% 
      dplyr::distinct()

    # Here we join all the dataframes together by using STRING_id
    # It would be the same of doing:
    # df_mapping <- df_mapping_orig %>% 
    #   dplyr::left_join(df_degree, by="STRING_id") %>% 
    #   dplyr::left_join(df_betweenness, by="STRING_id") %>%
    #   dplyr::left_join(df_closeness, by="STRING_id") %>% 
    #   dplyr::left_join(perc_betweenness_df, by="STRING_id") %>% 
    #   ...
    
    df_mapping <- list(df_mapping_orig, 
                    df_degree, 
                    df_betweenness,
                    df_closeness,
                    df_perc_degree,
                    df_perc_betweenness,
                    df_perc_closeness) %>% 
      purrr::reduce(left_join, by="STRING_id") %>% 
      dplyr::group_by(genesymbol) %>%
      dplyr::slice_max(order_by = degree, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()

    # Write data to the dataframe with all the cambinations oil/bact
    df_combinations$targets_found[df_combinations$oil_eng==oil &
                                    df_combinations$bacteria==bact] <- nrow(df_mapping)
    
    # We extract a unique score for defining the centrality
    # which accounts for the three terms closeness, betweenness and degree
    df_centrality <- df_mapping %>%
      dplyr::mutate(
        closeness_norm = normalize_as_cytoscape(closeness),
        betweenness_norm = normalize_as_cytoscape(betweenness),
        degree_norm = normalize_as_cytoscape(degree)
      ) %>% 
      dplyr::mutate(
        avg_centrality = (closeness_norm+betweenness_norm+degree_norm)/3
      )
    
    # We store all the targets (searched) with the data found 
    # whereas in the other list we store only the data found
    network_targets_list[[idx]] <- sel_targets_df %>% 
      dplyr::left_join(df_centrality, by=c("genesymbol", "oil_eng", 
                                           "bacteria", "sum_contr")) %>% 
      dplyr::select(-STRING_id)
    
    names_to_export[idx] <- sprintf("%s %s%s", oil, 
                                    substring(strsplit(bact, " ")[[1]][1], 1, 1),
                                    strsplit(bact, " ")[[1]][2])
  }
}

# Here we have the number of targets found
df_combinations$perc_targets_found <- 100*df_combinations$targets_found/df_combinations$nr_targets



# 
# The sheet Debug is written only in the debug mode (set it in the config.R file)
# 
# Note: used the filter on zzscore_molpoc to limit the number of rows 
# (error message: the xlsx format does not support tables with 1M+ rows)
# 
if (config$debug) {
  writexl::write_xlsx(
    x = list(
      Debug = df_targets_final,
      FullData = df_full %>% dplyr::filter(zzscore_molpoc>0),
      AllOils = df_targets_filtered,
      Oregano = df_targets_filtered %>%dplyr::filter(oil_eng=="Oregano"),
      Thyme = df_targets_filtered %>%dplyr::filter(oil_eng=="Thyme"),
      Cinnamon = df_targets_filtered %>%dplyr::filter(oil_eng=="Cinnamon"),
      Targets = df_combinations,
      Centrality = dplyr::bind_rows(network_targets_list)
    ),
    path = file.path(output_dir, "Results_targetfishing.xlsx")
  )
} else {
  writexl::write_xlsx(
    x = list(
      FullData = df_full %>% dplyr::filter(zzscore_molpoc>0),
      AllOils = df_targets_filtered,
      Oregano = df_targets_filtered %>%dplyr::filter(oil_eng=="Oregano"),
      Thyme = df_targets_filtered %>%dplyr::filter(oil_eng=="Thyme"),
      Cinnamon = df_targets_filtered %>%dplyr::filter(oil_eng=="Cinnamon"),
      Targets = df_combinations,
      Centrality = dplyr::bind_rows(network_targets_list)
    ),
    path = file.path(output_dir, "Results_targetfishing.xlsx")
  )
}


#
# This is to check which top central nodes have alkso good sum_contr
# (as reported in Table 4)
# 
centrality_df <- dplyr::bind_rows(network_targets_list)

for (wanted_oil in c("Thyme", "Oregano", "Cinnamon")) {
  for (wanted_bacteria in c("Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus")) {
    
    df_tmp <- centrality_df %>% 
      dplyr::filter(oil_eng==wanted_oil, bacteria==wanted_bacteria) %>%
      dplyr::arrange(desc(sum_contr)) %>% 
      dplyr::group_by(oil_eng,bacteria) %>% 
      dplyr::slice_head(prop = 0.2) %>% 
      dplyr::ungroup()
    
    df_tmp_centr <- centrality_df %>% 
      dplyr::filter(oil_eng==wanted_oil, bacteria==wanted_bacteria) %>% 
      dplyr::arrange(desc(avg_centrality)) %>% 
      dplyr::group_by(oil_eng,bacteria) %>% 
      dplyr::slice_head(n = 10) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(check = genesymbol %in% df_tmp$genesymbol)
    
    cat("\n\nOIL:", wanted_oil, "BACTERIUM:", wanted_bacteria)
    df_tmp_centr %>% 
      dplyr::filter(check) %>% 
      dplyr::select(genesymbol) %>% 
      print()
    
  }
}


