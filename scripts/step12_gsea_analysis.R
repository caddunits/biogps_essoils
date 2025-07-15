# 
# This script run the first part of the analysis: pathways enrichment
# 
# Input are data from .RDS files
# Output are tabular data (as .XLSX) and a unique heatmap oils vs pathways
# 

library(dplyr)
library(magrittr)
library(tibble)
library(circlize)
library(writexl)
library(tidyr)
library(purrr)
library(igraph)
suppressPackageStartupMessages(library(ComplexHeatmap))

library(clusterProfiler)

.Random.seed <- 2025



# 
# Functions
# 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Function to extract the genes corresponding to the same pocket
# (for a given oil-bacteria pair)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
find_genes_same_pocket <- function(df, bacteria, oil) {
  
  # Step 1: Extract relevant data (wanted oil and bacteria)
  #         and focus on genes which are linked to the same pocket
  df_1 <- df %>% 
    dplyr::filter(bacteria==bacteria, oil==oil) %>% 
    dplyr::select(genesymbol, pocketname) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(pocketname) %>% 
    dplyr::add_count() %>% 
    dplyr::filter(n>1) %>% 
    dplyr::group_by(pocketname) %>%
    dplyr::mutate(
      concatenated_genesymbol = paste(genesymbol, collapse = "/")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(genesymbol, concatenated_genesymbol) %>%
    dplyr::rename(gene=genesymbol) %>%
    dplyr::distinct() %>%
    dplyr::group_by(gene) %>%
    dplyr::add_count() %>%
    dplyr::mutate(
      concatenated2 = paste(concatenated_genesymbol, collapse = "/")
    )
  
  # Step 2: find edges (link between nodes = corresponding genes)
  df_2 <- df_1 %>%
    dplyr::mutate(gene_list = strsplit(concatenated2, "/")) %>% # Split into individual genes
    dplyr::pull(gene_list) %>%                                  # Extract the list of gene groups
    purrr::map(~ combn(.x, 2, simplify = FALSE)) %>%            # Generate pairwise combinations (only for groups with >=2 genes)
    purrr::discard(is.null) %>%                                 # Remove NULL entries (groups with only one gene)
    unlist(recursive = FALSE)                                   # Flatten to a list of pairs
  
  if (is.null(df_2)) {
    final <- NULL
    
  } else {
    edges <- df_2 %>%
      do.call(rbind, .)                                         # Convert to a matrix of edges
    
    # Step 3: Build a graph from the edges
    graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
    
    # Step 4: Identify clusters (connected components)
    clusters <- igraph::components(graph)
    
    # Step 5: Create a dataframe with the cluster assignments
    result <- data.frame(
      gene = names(clusters$membership), # Individual gene names
      group = clusters$membership        # Cluster ID
    ) %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        all_genes = paste(sort(unique(gene)), collapse = "/") # Create group name
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        gene = strsplit(all_genes, "/") # Expand the genes in each cluster
      ) %>%
      tidyr::unnest(gene) # One row per gene
    
    # Step 6: Arrange final output
    final <- result %>%
      dplyr::select(concatenated_genesymbol = all_genes, genesymbol = gene)
  }
  
  return(final)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Function to extract the data for given oil-bacteria pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
extract_data <- function(biogps_data, bact_df, comp_df, 
                         sel_bact, thr_comp) {
  
  orig_df <- biogps_data$full
  df_zscore_mol <- biogps_data$mol
  df_zscore_poc <- biogps_data$poc
  
  cat("READ COMPOSITION DATA\n")
  oils <- unique(comp_df$oil_eng)
  
  data_df <- do.call('rbind', lapply(oils, function(selected_oil) {
    
    cat("\nWORK ON OIL", selected_oil, "\n")
    df_oil <- comp_df %>% 
      dplyr::filter(oil_eng==selected_oil) %>% 
      dplyr::select(molname, fullname, perc_comp, phytoclass) %>%
      dplyr::filter(perc_comp > thr_comp)
    
    cat("JOIN WITH BACTERIA DATA\n")
    # Refine data (only data for: given bacteria and molecules of the oil)
    df <- orig_df %>% 
      dplyr::filter(molname %in% df_oil$molname) %>% 
      dplyr::left_join(bact_df, by="organism") %>% 
      dplyr::filter(bacteria == sel_bact) %>% 
      dplyr::select(-id_pdb, -info,
                    -bact_details, -bact_ATCC, -KEGG_code) %>% 
      # dplyr::select(-Notes, -id_pdb, -info,
      #               -KEGG_Tnumber, -bact_details, -bact_ATCC, -KEGG_code) %>% 
      dplyr::distinct() %>%
      dplyr::arrange(genesymbol)
    
    cat("MERGE WITH BIOGPS DATA\n")
    filtered <- df %>% 
      dplyr::left_join(df_zscore_mol, by=c('id_mol', 'id_poc'), relationship = "many-to-many") %>%
      dplyr::left_join(df_zscore_poc, by=c('id_mol', 'id_poc', 'poc_refcode'), relationship = "many-to-many") %>% 
      # dplyr::mutate(zzscore_molpoc_ref = ifelse(zzscore_molpoc > config$threshold_zzscore, 
      #                                           zzscore_molpoc, 0)) %>% 
      dplyr::mutate(zzscore_molpoc_ref = ifelse(
        zzscore_molpoc < config$threshold_zzscore |
          zscore_mol < config$threshold_zscore_mol |
          zscore_poc < config$threshold_zscore_poc,
        0, zzscore_molpoc)) %>%
      dplyr::left_join(df_oil, by=c("molname")) %>% 
      dplyr::select(-id_poc, -id_mol, -poc_refcode) %>% 
      dplyr::mutate(perc_weight = ifelse(perc_comp>1, 1+log10(perc_comp), perc_comp)) %>%
      dplyr::mutate(contribution_ref = zzscore_molpoc_ref * perc_weight) %>%
      dplyr::mutate(contribution_orig = zzscore * perc_weight) %>% 
      dplyr::mutate(oil=selected_oil) %>% 
      dplyr::select(bacteria, organism, oil, 
                    genesymbol, uniprot, pdb, pocketname, nr_pockets,
                    molname, fullname, phytoclass, perc_comp, perc_weight,
                    zzscore, zscore_mol, zscore_poc, zzscore_molpoc_ref, 
                    contribution_ref, contribution_orig)
  }))
  
  return(data_df)
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Function to run gsea from the package clusterProfiler
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
safe_gsea <- function(geneList, pathway, perm, pvalue_threshold, verbose) {
  tryCatch(
    {
      result <- clusterProfiler::GSEA(
        geneList = geneList,
        TERM2GENE = pathway,
        seed = TRUE,
        by = "fgsea",
        maxGSSize= 100000000000,
        minGSSize = 1,
        nPermSimple = perm,
        pvalueCutoff = pvalue_threshold,
        pAdjustMethod = "none",
        verbose = verbose
      )
      return(result) # Return the result if successful
    },
    error = function(e) {
      message("Error encountered: ", e$message) # Print the error message
      return(NULL) # Return NULL in case of an error
    }
  )
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Function to draw a heatmap with colog_pvalues and text
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
create_complex_heatmap_pathways <- function(df, output_filename) {
  
  mycols <- circlize::colorRamp2(
    breaks = c(0, 4), 
    colors = c("white", "blue")
  )
  
  # Define colors for column labels
  column_label_colors <- c(
    Cinnamon = "#B66E34", 
    Oregano = "#244A0D", 
    Thyme = "#95C977"
  )
  # RGB 
  # Cinnamon: 182, 110,  52
  # Oregano:   36,  74,  13
  # Thyme:    149, 201, 119
  
  # The input dataframe has the following columns:
  # bacteria, oil, pathway, colog_pvalue, text
  # bacteria will be used to split the heatmaps on Y axis
  # oil will be used as columns
  # pathway will be used as rows
  # colog_pvalue will be used to colour the cells
  # text will be used to write content in the cells
  
  # As first step we need to create two matrix that have
  # the same order of colog_pvalue and text
  # extract data as dataframe in the right format (using pivot)
  df_with_values <- df %>% 
    dplyr::select(-text) %>%
    tidyr::pivot_wider(names_from = oil, values_from = colog_pvalue)
  # Transform data as matrix, set zero to NA values and assign the rownames
  # N.B. only numerical data; we suppose to have the first two columns as info
  xmatrix_values <- as.matrix(df_with_values[,c(3:ncol(df_with_values))])
  xmatrix_values[is.na(xmatrix_values)] <- 0
  row.names(xmatrix_values) <- df_with_values$pathway

  # Then create a corresponding matrix with text values 
  # Note we have the right format (using pivot) as the matrix with values
  df_with_text <- df %>% 
    dplyr::select(-colog_pvalue) %>%
    tidyr::pivot_wider(names_from = oil, values_from = text)
  # Transform data as matrix and assign the rownames
  # N.B. again, we suppose to have the first two columns as info
  xmatrix_text <- as.matrix(df_with_text[,c(3:ncol(df_with_text))])
  row.names(xmatrix_text) <- df_with_text$pathway
  

  # We also create a vector with bacteria data, in order to use as split  
  x_bacteria <- df_with_values$bacteria
  
  hm1 <- Heatmap(
    xmatrix_values, 
    
    # Legend
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = expression(-log[10]("p-value")),
      direction = "horizontal"
    ),
    
    # Rows
    show_row_names = TRUE, 
    row_names_gp = gpar(col = "black", fontsize = 8),
    
    # split = df$bacteria,
    split = x_bacteria,
    
    row_title_rot = 0,
    cluster_row_slices = FALSE,
    row_gap = unit(3, "mm"),
    
    # Columns
    show_column_names = TRUE,
    
    cluster_columns = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(
      fontsize = 12, 
      fontface = "bold",
      col = column_label_colors[colnames(xmatrix_values)]  # Set label colors
    ), 
    
    # Cells
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(xmatrix_text[i, j], x, y, gp = gpar(fontsize = 7, col="white"))
    },
    
    width = 12,
    # heatmap_width = unit(24, "cm"),
    border = TRUE,
    col = mycols
    
  )
  
  svg(output_filename)
  draw(hm1, heatmap_legend_side = "bottom")
  dev.off()
  
  return()
  
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


repo_dir <- getwd()

# Settings are set once, in a configuration file
source(file.path(repo_dir, "scripts", "config.R"))


# Input
data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")

# Output
output_dir <- file.path(repo_dir, "output", "gsea_pathways_enrichment")

# Create the ouput directory if it does not exist
dir.create(output_dir, showWarnings = FALSE)

# 
# First, read from .RDS files the data we need in the various steps
# 

biogps_data <- readRDS(file = file.path(data_dir, "Biogps_zscore_data", "biogps_data.RDS"))

df_kegg_bacteria_unified <- readRDS(file = file.path(mapp_dir, "kegg_bacteria_unified.RDS")) %>% 
  dplyr::filter(!grepl("^NA ", pathway_name)) # Note: this is to correct an error in input file pap

df_bacteria <- readRDS(file = file.path(mapp_dir, "bacteria_data.RDS"))

df_composition_new <- readRDS(file = file.path(mapp_dir, "composition_new.RDS"))

# Set a variable to decide whether to use data from manual curation
manual_curation <- TRUE
# manual_curation <- FALSE

# 
# Read data about multiple genes and manual curation
# 
fullpath_datacuration <- file.path(data_dir, "Curation_data", "manual_curation.xlsx")

df_multiple <- readxl::read_xlsx(path=fullpath_datacuration,
                                 sheet = "multiple_genes")

df_curation_pdb <- readxl::read_xlsx(path=fullpath_datacuration,
                                     sheet = "pdb_curation")

df_curation_pockets <- readxl::read_xlsx(path=fullpath_datacuration,
                                         sheet = "pocket_curation")

df_curation_targets <- readxl::read_xlsx(path=fullpath_datacuration,
                                         sheet = "target_curation")


# 
# Loop over all the oils and all the bacteria
# 

# Set an empty dataframe
df_pathways <- data.frame(
  bacteria = character(0),
  oil = character(0),
  pathway = character(0),
  nrP = numeric(0),
  set_size = numeric(0),
  enrichment_score = numeric(0),
  nes = numeric(0),
  pvalue = numeric(0),
  padjust= numeric(0),
  leading_edge = character(0),
  core_enrichment = character(0),
  flag_nrP = logical(0),
  flag_pvalue = logical(0),
  stringsAsFactors = FALSE
)

df_threshold_max_contrib_orig <- data.frame(
  oil = character(0),
  bacteria = character(0),
  value = numeric(0),
  stringsAsFactors = FALSE
)

# List of pathwys for debug 
# debug_pathways <- c(
#   "Fatty acid metabolism",
#   "Fatty acid biosynthesis",
#   "Tryptophan metabolism"
# )

# Create NULL dataframes where we will add information about: 
# overall contributions for each phytocomplex-bacteria pair
# single contributions (for each phytocomplex-bacteria/gene) from molecules
df_phytocomplexes <- NULL
df_moltargets_contributions <- NULL

pathways_and_genes <- list()
best_data <- list()
full_data <- list()

# Loop over the bacteria
for (wanted_bacteria in config$wanted_bacteria) {

  df_setK_gene <- df_kegg_bacteria_unified %>% 
    dplyr::filter(bacteria==wanted_bacteria) %>% 
    dplyr::select(bact_kegg, bacteria, gene) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(gene) %>% 
    dplyr::add_count() %>% 
    dplyr::arrange(gene, bact_kegg)
  
  genes_setK <- unique(df_setK_gene$gene)

  # Extract data 
  df_data_new <- extract_data(
    biogps_data = biogps_data, 
    bact_df = df_bacteria, 
    comp_df = df_composition_new, 
    sel_bact = wanted_bacteria, 
    thr_comp = config$threshold_comp
  )
  
  # 
  # Here we have to exclude data from manual curation
  # 
  if (manual_curation) {
    
    # Here we exclude rows for unwanted genesymbols
    df_curation_targets %<>% dplyr::filter(curation=="exclude")
    df_data_new <- df_data_new %>%
      dplyr::anti_join(df_curation_targets, by=c("bacteria", "genesymbol"))
    
    # Here we exclude rows for unwanted PDB entries or unwanted pocketnames
    df_curation_pdb_exclude_list <- df_curation_pdb %>% 
      dplyr::filter(curation=="exclude") %>% dplyr::select(pdb)
    
    df_curation_pockets_exclude_list <- df_curation_pockets %>% 
      dplyr::filter(curation=="exclude") %>% dplyr::select(pocketname)
    
    df_data_new <- df_data_new %>%
      dplyr::filter(!pdb %in% df_curation_pdb_exclude_list$pdb) %>%
      dplyr::filter(!pocketname %in% df_curation_pockets_exclude_list$pocketname)
  }
    
  # Store data in a list, in order to export as XLSX file 
  # Note that to avoid problems with the excel limits we only keep data with 
  # positive contribution_ref values (i.e. exclude all the rows with 0)
  best_data[[wanted_bacteria]] <- df_data_new %>% 
    dplyr::filter(contribution_ref > 0) %>% 
    dplyr::arrange(desc(contribution_ref)) %>%
    dplyr::group_by(bacteria, oil, genesymbol, molname) %>%
    dplyr::slice_head(n = 1)
  
  full_data[[wanted_bacteria]] <- df_data_new
    
  setB_genelist_all_oils <- unique(df_data_new$genesymbol)
  
  # Retain only genes that are in KEGG-pathways
  setB_genelist_in <- setB_genelist_all_oils[which(setB_genelist_all_oils %in% genes_setK)]


  # Pathways for the given bacteria
  df_setP <- df_kegg_bacteria_unified %>% 
    dplyr::filter(bacteria==wanted_bacteria) %>% 
    dplyr::select(pathway_name, gene, organism) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(pathway_name, gene) %>%
    dplyr::add_count() %>%
    dplyr::arrange(desc(n), gene)
  
  pathways_setP <- unique(df_setP$pathway_name)
  
  oils <- unique(df_data_new$oil)
  
  
  # Main loop on oils
  # oils = c()
  for (wanted_oil in oils) {
    cat("SELECTED OIL:", wanted_oil, "\n")
    
    multiplegenes_df <- find_genes_same_pocket(
      df = df_data_new,
      bacteria = wanted_bacteria,
      oil = wanted_oil
    )

    # 
    # We first rank the genes based on phytocomplex
    # 
    
    # Attention: special case of NULL results when searching for correspondences 
    if (is.null(multiplegenes_df)) {
      targets_phytocomplex_df <- df_data_new %>% 
        dplyr::filter(bacteria==wanted_bacteria, oil==wanted_oil) %>% 
        dplyr::select(genesymbol, molname, contribution_ref, contribution_orig) %>%
        dplyr::group_by(genesymbol, molname) %>%
        dplyr::summarise(max_contribution_ref = max(contribution_ref),
                         max_contribution_orig = max(contribution_orig)) %>% 
        dplyr::group_by(genesymbol) %>%
        dplyr::summarise(sum_contr_ref = sum(max_contribution_ref),
                         sum_contr_orig = sum(max_contribution_orig)) %>%
        dplyr::arrange(desc(sum_contr_ref), desc(sum_contr_orig)) %>%
        dplyr::distinct()
    } else {
      targets_phytocomplex_df <- df_data_new %>% 
        dplyr::filter(bacteria==wanted_bacteria, oil==wanted_oil) %>% 
        # dplyr::filter(!genesymbol %in% exclude_list$genesymbol) %>%                # QUESTA NON L'HO ANCORA SISTEMATO (PROBABILMENTE EFFETTO MINORE)
        dplyr::select(genesymbol, molname, contribution_ref, contribution_orig) %>%
        dplyr::left_join(multiplegenes_df, by="genesymbol") %>%
        dplyr::mutate(genesymbol = ifelse(is.na(concatenated_genesymbol),
                                          genesymbol, concatenated_genesymbol)) %>%
        dplyr::select(-concatenated_genesymbol) %>%
        dplyr::group_by(genesymbol, molname) %>%
        dplyr::summarise(max_contribution_ref = max(contribution_ref),
                         max_contribution_orig = max(contribution_orig)) %>% 
        dplyr::group_by(genesymbol) %>%
        dplyr::summarise(sum_contr_ref = sum(max_contribution_ref),
                         sum_contr_orig = sum(max_contribution_orig)) %>%
        dplyr::arrange(desc(sum_contr_ref), desc(sum_contr_orig)) %>%
        dplyr::distinct()
    }
    
    
    # We need a final ranked list, that will be obtained by using a new value
    # (first, identify the max value (for orig contribution) among the
    #  zeroes or negative values)
    max_contr_orig_on_zeros <- targets_phytocomplex_df %>% 
      dplyr::filter(sum_contr_ref <= 0) %>% 
      dplyr::mutate(x=max(sum_contr_orig)) %>% 
      dplyr::select(x) %>% 
      dplyr::distinct()

    addendum <- max_contr_orig_on_zeros$x
    cat("ADDENDUM DEFINED:", addendum, "\n")
    
    # Calculate the new contributions in a unique scale (by adding addendum) 
    df_new_targets_phytocomplex <- targets_phytocomplex_df %>% 
      dplyr::mutate(new_contr = ifelse(sum_contr_ref > 0, 
                                       sum_contr_ref + addendum, 
                                       sum_contr_orig)) %>% 
      dplyr::arrange(desc(new_contr)) %>% 
      dplyr::mutate(bacteria=wanted_bacteria,
                    oil=wanted_oil)
    
    # Add data with the contribution to a unique dataframe
    if (is.null(df_phytocomplexes)) {
      df_phytocomplexes <- df_new_targets_phytocomplex
    } else {
      df_phytocomplexes <- rbind(
        df_phytocomplexes,
        df_new_targets_phytocomplex
      )
    }
    
    # Extract the molecular contribution for each bacteria-oil-gene 
    df_moltargets <- df_data_new %>% 
      dplyr::filter(oil==wanted_oil, 
                    bacteria==wanted_bacteria,
                    contribution_ref>0) %>% 
      dplyr::select(bacteria, oil, 
                    genesymbol, uniprot, pdb, pocketname, nr_pockets,
                    molname, fullname, phytoclass,
                    perc_comp, perc_weight, zzscore, zscore_mol, zscore_poc,
                    zzscore_molpoc_ref, contribution_ref) %>% 
      dplyr::group_by(bacteria, oil, genesymbol, fullname) %>% 
      dplyr::summarise(max_contribution_ref = max(contribution_ref))
    
    
    # Add data with the mol/target contribution to a unique dataframe
    if (is.null(df_moltargets_contributions)) {
      df_moltargets_contributions <- df_moltargets
    } else {
      df_moltargets_contributions <- rbind(
        df_moltargets_contributions,
        df_moltargets
      )
    }
    
    # Add information on the addendum to a unique dataframe
    df_threshold_max_contrib_orig %<>% add_row(
      bacteria = wanted_bacteria, 
      oil = wanted_oil, 
      value = addendum
    )
    
    gene_list <- c(df_new_targets_phytocomplex$new_contr)
    names(gene_list) <- df_new_targets_phytocomplex$genesymbol
    
    # Loop over the pathways
    for (p in 1:length(pathways_setP)) {
      
      cat("\n--------------------------------------------------\n")
      cat("BACTERIA:", wanted_bacteria, "OIL:", wanted_oil, 
          "PATHWAY ID:", p, "/", length(pathways_setP), "PATHWAY NAME:", pathways_setP[p], "\n")

      # if (!pathways_setP[p] %in% debug_pathways) next
      
      # Genes of the pathway
      subsetP_genes <- unique(df_setP$gene[df_setP$pathway_name == pathways_setP[p]])
      nrP <- length(subsetP_genes)

      # setB_genelist_in IS THE LIST OF GENES THAT ARE IN THE KEGG PATHWAYS
      found_in_setB <- which(subsetP_genes %in% setB_genelist_in)

      if (length(found_in_setB)==0) {
        cat("\tPATHWAY'S GENES NOT FOUND\n")
        
      } else {
        # Table to be used in GSEA-Like analysis, with pathway and genes
        # (Also here we have to join correspondences)
        
        # Attention: special case of NULL results when searching for correspondences 
        if (is.null(multiplegenes_df)) {
          df_associated_genes <- df_setP %>% 
            dplyr::filter(pathway_name == pathways_setP[p]) %>%
            dplyr::select(pathway_name, gene) %>% 
            dplyr::distinct()
        } else {
          df_associated_genes <- df_setP %>% 
            dplyr::filter(pathway_name == pathways_setP[p]) %>%
            dplyr::select(pathway_name, gene) %>% 
            dplyr::left_join(multiplegenes_df, by=c("gene"="genesymbol")) %>% 
            dplyr::distinct() %>% 
            dplyr::mutate(gene = ifelse(is.na(concatenated_genesymbol),
                                        gene, concatenated_genesymbol)) %>%
            dplyr::select(-concatenated_genesymbol)
        }
          
        # Add information of pathways and associated genes for wanted bacteria
        # NOTE: we know this is repeated several times (as many times as oils 
        # we have). At the end we will apply the distinct() to avoid duplicates
        if (!wanted_bacteria %in% names(pathways_and_genes)) {
          pathways_and_genes[[wanted_bacteria]] <- df_associated_genes %>% 
            dplyr::mutate(bacteria = wanted_bacteria) %>% 
            dplyr::rename(genesymbol=gene, pathway=pathway_name)
        } else {
          pathways_and_genes[[wanted_bacteria]] <- rbind(
            pathways_and_genes[[wanted_bacteria]],
            df_associated_genes %>% 
              dplyr::mutate(bacteria = wanted_bacteria) %>% 
              dplyr::rename(genesymbol=gene, pathway=pathway_name)
          )
        }
        
        cat("NUMBER OF GENES OF THE PATHWAY:", nrow(df_associated_genes), "\n")
        if (config$verbose) print(df_associated_genes$gene)

        cat("- - - \n")
        cat("CHECK IF AVAILABLE IN THE GENE LIST (RANKED BY PHYTOCOMPLEX):\n")
        check_gene_list <- gene_list[names(gene_list) %in% subsetP_genes[found_in_setB]]
        nr_genes_found <- length(check_gene_list)
        print(check_gene_list)
        cat("NR OF GENES FOUND:", nr_genes_found, "\n")
        cat("- - - \n")

        if (nr_genes_found == 0) {
          cat("\tNO GENE FOUND\n")
          
        } else if (all(check_gene_list < addendum)) {
          cat("\tALL VALUES BELOW THRESHOLD (ADDENDUM)\n")
          
        } else {
          gse <- safe_gsea(
            geneList = gene_list, 
            pathway = df_associated_genes, 
            perm = config$nr_permutations, 
            pvalue_threshold = 1, # pvalue_threshold CAMBIATO SABATO 25/1 (era 0.05)
            verbose = config$verbose
          )

          if (is.null(gse)) {
            cat ("ERRORS\n")
            
          } else {
            # Extract results into a dataframe
            results_nrow <- nrow(gse@result)
            
            if (results_nrow == 0) {
              cat ("NO RESULTS\n")
              
            } else {
              cat("GSEA ANALYSIS DONE\n")
              print(gse@result$pvalue)
              
              df_single_pathway <- data.frame(
                bacteria = wanted_bacteria,
                oil = wanted_oil,
                pathway = pathways_setP[p],
                nrP = nrP,
                set_size = gse@result$setSize,
                enrichment_score = gse@result$enrichmentScore,
                nes = gse@result$NES,
                pvalue = gse@result$pvalue,
                padjust = gse@result$p.adjust,
                leading_edge = gse@result$leading_edge,
                core_enrichment = gse@result$core_enrichment,
                stringsAsFactors = FALSE
              ) %>%             
                dplyr::mutate(  
                  flag_nrP = nrP >= config$threshold_nrP,            # Add a flag for pathways with high number of genes
                  flag_pvalue = pvalue <= config$threshold_pvalue    # Add a flag for pvalue (true if below the threshold)
                )
              
              df_pathways <- rbind(
                df_pathways,
                df_single_pathway
              )
            }
          
          }

        }

      }

    }
    
  }
  
}

# Convert data available as lists into dataframes (to be easily handled by join)
df_pathways_and_genes <- dplyr::bind_rows(pathways_and_genes, .id = "bacteria") %>% 
  distinct()

df_best_data <- dplyr::bind_rows(best_data, .id = "bacteria")
df_full_data <- dplyr::bind_rows(full_data, .id = "bacteria")

# Filter the pathways (only with nrP below threshold and significant pvalue)
df_pathways_filtered <- df_pathways %>% 
  dplyr::filter(flag_nrP==FALSE, flag_pvalue==TRUE)


# Create a dataframe with both columns for genesymbol, one
# with original values and one with split
df_phytocomplexes_split <- df_phytocomplexes %>% 
  dplyr::mutate(genesymbol_split=genesymbol) %>% 
  tidyr::separate_rows(genesymbol_split, sep = "/")

# For each significant pathway we add the corresponding genes (targets) 
df_pathways_and_targets <- df_pathways_filtered %>% 
  dplyr::select(-flag_nrP, -flag_pvalue, -enrichment_score, -nes, -padjust) %>% 
  dplyr::left_join(df_pathways_and_genes, 
                   by=c("bacteria", "pathway"),
                   relationship = "many-to-many")


# Add information about phyocomplexes, i.e. three columns 
# for sum contributions of each target. Note that we join
# by using genesymbol where complexes are included (ex: ACCA/ACCD)
# and corresponding rows are duplicated
df_pathways_targets_contrib <- df_pathways_and_targets %>% 
  dplyr::left_join(df_phytocomplexes_split,
                   by=c("bacteria", "oil", "genesymbol"),
                   relationship = "many-to-many") %>% 
  dplyr::mutate(gene_in_core = stringr::str_detect(core_enrichment,
                paste0("\\b", genesymbol_split, "\\b"))) %>% 
  dplyr::arrange(bacteria, oil, pathway, desc(new_contr), genesymbol)


# However, in order to count how many targets of different classes 
# (STRONG/WEAK/NONE) are associated to the pathway we need to group again
# in order that complexes are counted only once. First we skip genes not in core
df_pathways_targets_interaction <- df_pathways_targets_contrib %>% 
  dplyr::filter(gene_in_core==TRUE) %>% 
  dplyr::arrange(desc(new_contr)) %>% 
  dplyr::select(bacteria, oil, pathway, genesymbol, 
                sum_contr_ref, sum_contr_orig, new_contr) %>%
  dplyr::distinct() %>% 
  dplyr::mutate(class_interaction = ifelse(sum_contr_ref > 0, "STRONG",
                                          ifelse(sum_contr_orig > 0,
                                          "WEAK", "NONE")))

# Finally we count how many targets have interaction of the different classes
# We will have three columns for STRONG, WEAK and NONE
df_pathways_targets_counters <- df_pathways_targets_interaction %>% 
  dplyr::group_by(bacteria, oil, pathway) %>%
  dplyr::summarise(
    nr_strong = sum(class_interaction == "STRONG"),
    nr_weak = sum(class_interaction == "WEAK"),
    nr_none = sum(class_interaction == "NONE")
  ) 

# Add a flag to check whether pathways of class NONE are more
# than pathways of class STRONG or WEAK
df_pathways_flag_strength <- df_pathways_targets_counters %>% 
  dplyr::mutate(flag_strength = nr_none<=nr_strong&nr_none<=nr_weak)

# Focus on the molecules that interact with the selected targets
df_pathways_molecules_interaction <- df_pathways_targets_contrib %>% 
  dplyr::filter(sum_contr_ref > 0) %>% 
  dplyr::left_join(df_best_data,
                   by=c("bacteria", "oil", "genesymbol_split"="genesymbol"),
                   relationship = "many-to-many") %>%
  dplyr::filter(!is.na(molname)) %>% 
  dplyr::arrange(desc(new_contr), desc(contribution_ref)) %>% 
  dplyr::select(bacteria, oil, pathway, genesymbol, genesymbol_split, 
                uniprot, pdb, pocketname, nr_pockets, 
                molname, fullname, phytoclass, perc_comp, perc_weight,
                zzscore, zscore_mol, zscore_poc, zzscore_molpoc_ref,
                contribution_ref, contribution_orig, 
                sum_contr_ref, sum_contr_orig, new_contr)  
  
# Finally we count how many molecules have interactions, focusing only
# on the class STRONG. First we extract data
df_pathways_molecules_targets <- df_pathways_molecules_interaction %>% 
  dplyr::left_join(df_pathways_targets_interaction, 
                   by=c("bacteria", "oil", "pathway","genesymbol"))

# and the we derive the counter
df_pathways_molecules_counters <- df_pathways_molecules_targets %>% 
  dplyr::group_by(bacteria, oil, pathway) %>% 
  dplyr::summarise(nr_molecules = dplyr::n_distinct(fullname), .groups = "drop")

# Put together the counters for molecules and targets
df_pathways_counters <- df_pathways_targets_counters %>% 
  dplyr::left_join(df_pathways_molecules_counters, 
                   by=c("bacteria", "oil", "pathway"))

# Create a dataframe with the content for the cells
df_pathways_celltext <- df_pathways_counters %>% 
  dplyr::mutate(text = sprintf("%d / %d", nr_strong, nr_molecules)) %>% 
  dplyr::select(bacteria, oil, pathway, text)

# 
# NOTE: we decided to write the number of targets (STRONG) and of molecules
# but we could also write something else 
# nr_strong / nr_weak / nr_none
# 


# 
# Heatmaps bacteria vs pathways
# 
# A unique heatmap includes all the bacteria (this is because heatmaps for 
# single rows or single columns are not created by heatmaply)
# 

ceiling <- 10

df_pathways_final <- df_pathways %>%
  dplyr::left_join(df_pathways_flag_strength, 
                   by=c("bacteria", "oil", "pathway")) %>% 
  dplyr::filter(flag_pvalue==TRUE, flag_nrP==FALSE, flag_strength==TRUE) %>%
  # dplyr::filter(flag_pvalue==TRUE, flag_nrP==FALSE) %>%
  dplyr::mutate(colog_pvalue = -1*log10(pvalue)) %>%
  dplyr::mutate(colog_pvalue = ifelse(colog_pvalue>ceiling, ceiling, colog_pvalue) )


# New dataframe that has not only the colog_pvalue but also the text
# that should be reported in each cell
df_pathways_heatmap_text <- df_pathways_final %>% 
  dplyr::select(-nrP, -set_size, -enrichment_score, -nes, -pvalue,
                -padjust, -leading_edge, -core_enrichment, 
                -flag_nrP, -flag_pvalue, -flag_strength,
                -nr_strong, -nr_weak, -nr_none) %>% 
  dplyr::left_join(df_pathways_celltext,
                   by=c("bacteria", "oil", "pathway"))

# Draw and save the heatmap
heatmap_filename <- file.path(output_dir, "Heatmap_3oils_allbacts_gsea.svg")
create_complex_heatmap_pathways(df = df_pathways_heatmap_text,
                                output_filename = heatmap_filename)


# 
# Write the results in an excel file and a .RDS file
# 

df_pathways_addedstrength <- df_pathways %>% 
  dplyr::left_join(df_pathways_flag_strength, 
                 by=c("bacteria", "oil", "pathway"))

writexl::write_xlsx(
  x = list(
    Pathways = df_pathways_addedstrength,         # Results from GSEA analysis
    Pathways_filtered = df_pathways_final,        # Filtered pathways, by flag_pvalue(TRUE), 
                                                  # flag_nrP(FALSE) and flag_strength(TRUE)
    Heatmap = df_pathways_heatmap_text,           # Same data reported in the heatmap
    Pathways_and_genes = df_pathways_and_genes,   # All the associations gene-pathway (for genes included in the study) 
    Addendum = df_threshold_max_contrib_orig,     
    Phytocomplexes = df_phytocomplexes,
    Targets_interaction = df_pathways_targets_interaction,
    Pathways_details = df_pathways_molecules_interaction,
    Pathways_counters = df_pathways_counters,   # Added 2025-06-09
    Moltargets = df_moltargets_contributions,
    Bestdata = df_best_data
  ),
  path = file.path(output_dir, "Results_pathwaysenrichment_gsea.xlsx")
)

# Save relevant data for further analysis (alluvial plots)
saveRDS(object=df_pathways_addedstrength, 
        file=file.path(mapp_dir, "all_pathways.RDS"))
# saveRDS(object=df_pathways, 
#         file=file.path(mapp_dir, "all_pathways.RDS"))
saveRDS(object=df_pathways_final, 
        file=file.path(mapp_dir, "sign_pathways.RDS"))
saveRDS(object=df_pathways_molecules_interaction, 
        file=file.path(mapp_dir, "molecules_interaction.RDS"))

# 
# AGGIUNGERE FILTRO SUI PATHWAYS BASATO SU STRONG,WEAK e NONE 
# (nuova variabile: flag_strength)
# 
