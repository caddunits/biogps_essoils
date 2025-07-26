library(readxl)
library(dplyr)
library(magrittr)
library(writexl)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)


# 
# FUNCTIONS --------------------------------------------------------------------
#
# Read a dataframe and extract only relevant data
extract_target_data <- function(df) {
  
  df1 <- df %>% 
    dplyr::select(gene, target_overall) %>% 
    dplyr::arrange(desc(target_overall)) %>% 
    dplyr::distinct()
  
  df2 <- df %>% 
    dplyr::select(gene, target_overall, centrality) %>% 
    dplyr::arrange(desc(target_overall)) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise(mean_centrality=mean(centrality)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct()
  
  df3 <- df1 %>% 
    dplyr::right_join(df2, by="gene")
  
  return(df3)
}

extract_perc_values <- function(df_comp, oil) {
  df_comp %>% 
    dplyr::filter(oil_eng == !!oil) %>%  # Attention used the symbol !! for using oil as variable
    dplyr::arrange(desc(perc_comp)) %>%
    dplyr::select(molecule=fullname, perc_comp)
}

create_empty_matrix <- function(df_mol, df_tar) {
  # Create an empty matrix with targets as rows and molecules as columns
  m <- matrix(0, nrow = nrow(df_tar), ncol = nrow(df_mol))
  colnames(m) <- df_mol$molecule
  rownames(m) <- df_tar$gene
  return(m)
}

fill_matrix_with_data <- function(bact, oil, m, df) {
  for (r in seq(1, nrow(df))) {
    if(df$bacteria[r] == bact &
       df$oil_eng[r] == oil) {
      loop_target <- df$gene[r]
      loop_molecule <- df$fullname[r]
      zzscore <- df$zzscore[r]
      if (config$debug) {
        cat("INDEX r", r, "TARGET", loop_target, "MOL",
            loop_molecule, "ZZSCORE", zzscore, "\n")
      }
      m[ which(rownames(m)==loop_target), 
         which(colnames(m)==loop_molecule) ] <- zzscore
    }
  }
  return(m)
}


# 
# Settings and input data ------------------------------------------------------
#
# From previous analysis we have found that data are unbalanced.
# We will produce heatmaps only for bacteria with the most abundant data
# We take the list of bacteria from config$major_bacteria

# 
# Setting Parameters
# 
# Number of top genes to be considered for each oil
# n_top_genes <- 15     # For the figures of the paper
n_top_genes <- 20     # For the figures of the paper
# n_top_genes <- 30   # For the figures of Supporting Information

# NOTE: The first iteration, with n_top_genes=30, provides a list of 
#       targets that are manually curated. 
#       The second iteration, with n_top_genes=15, works with a refined 
#       list of targets, after manual curation
manual_curation <- TRUE
# manual_curation <- FALSE

# Define manual colors for each oil
oil_colors <- c(
  "Cinnamon" = "#B66E34", 
  "Thyme" = "#95C977", 
  "Oregano" = "#244A0D"
)  


# Change this line according to the path of your repository
repo_dir <- getwd()

data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")

# Source the configuration file
source(file.path(repo_dir, "scripts", "config.R"))

# 
# NOTE: Below are checks about the existance of files from previous steps
#
# Read data from .RDS file (composition data)
if (file.exists(file.path(mapp_dir, "composition_new.RDS"))) {
  df_composition_new <- readRDS(file = file.path(mapp_dir, "composition_new.RDS"))
} else {
  stop("\nPlease check previous steps: problems with file composition_new.RDS")
}

# Read data from the previous step (xlsx file)
output_dir <- file.path(repo_dir, "output", "target_fishing")
fullpath_res_targetfishing <- file.path(output_dir, "Results_targetfishing.xlsx")

if (!file.exists(fullpath_res_targetfishing)) {
  stop("\nPlease check previous steps: problems with file Results_targetfishing.xlsx")
  
}

# 
# These targets are here manually modified 
# (based on lists of targets which correspond to the same multiple-target) 
# (Potential improvement in the future is to create this dataframe by scripting)
# 
# Escherichia coli
# MCBA/MCBB/MCBC/MCBD
# CARA/CARB
# 
# Pseudomonas aeruginosa
# AMIC/AMIR
# GSPI/GSPJ/GSPK
# NORB/NORC
# 
# Staphylococcus aureus
# ACCA/ACCD
# 
# fullpath_multiplegenes <- file.path(data_dir, "multiple_genes.xlsx")

# 
# Read data about multiple genes
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

df_curation_proteins <- readxl::read_xlsx(path=fullpath_datacuration,
                                          sheet = "proteins") %>% 
  dplyr::select(genesymbol, protein) %>% 
  dplyr::distinct()
# Note: given that targets name is the same for different bacteria
#       (for the same target) we can operate the given simplification above

# 
# Read data from excel about centrality: this is an average of different scores
# obtained in a previous step
# 
df_targets_contr_and_centrality <- readxl::read_xlsx(
  path=fullpath_res_targetfishing, sheet="Centrality") %>% 
  dplyr::filter(!is.na(degree)) %>% 
  dplyr::select(oil_eng, bacteria, genesymbol, sum_contr, avg_centrality) %>% 
  dplyr::rename(sum_contr_all = sum_contr)


# 
# Read data from excel about scores, only for major bacteria
# In addition we filter by zzscore and zscore values, considering the 
# thresholds defined in config.R (that are not set to 0, but should be 
# 2 for zzscore and 1 for zscore_mol and zscore_poc)
# 
df_targets_filtered <- readxl::read_xlsx(
  path=fullpath_res_targetfishing, sheet="AllOils") %>% 
  dplyr::filter(zzscore_molpoc > config$threshold_zzscore_barplot,
                zscore_mol > config$threshold_zscore_mol_barplot,
                zscore_poc > config$threshold_zscore_poc_barplot) %>%
  dplyr::filter(!is.na(bacteria)) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(bacteria, oil_eng, genesymbol, fullname) %>%                  # Added to solve a problem of multiple data for a few cases
  dplyr::slice_max(order_by = zzscore_molpoc, n = 1, with_ties = FALSE) %>%     # Added to solve a problem of multiple data for a few cases
  dplyr::ungroup() %>%                                                          # Added to solve a problem of multiple data for a few cases
  dplyr::rename(zzscore = zzscore_molpoc) %>% 
  dplyr::select(bacteria, oil_eng, fullname, dplyr::everything()) %>% 
  dplyr::arrange(oil_eng, bacteria, desc(zzscore), fullname)


# 
# 2025-05-26
# Here we can apply some corrections based on manual data curation
# 
if (manual_curation) {
  
  # Here we exclude rows for unwanted genesymbols
  df_curation_targets %<>% dplyr::filter(curation=="exclude")
  df_targets_filtered <- df_targets_filtered %>% 
    dplyr::anti_join(df_curation_targets, by=c("bacteria", "genesymbol"))
  
  # Here we exclude rows for unwanted PDB entries or unwanted pocketnames
  df_curation_pdb_exclude_list <- df_curation_pdb %>% 
    dplyr::filter(curation=="exclude") %>% dplyr::select(pdb)
  
  df_curation_pockets_exclude_list <- df_curation_pockets %>% 
    dplyr::filter(curation=="exclude") %>% dplyr::select(pocketname)
  
  df_targets_filtered <- df_targets_filtered %>%
    dplyr::filter(!pdb %in% df_curation_pdb_exclude_list$pdb) %>% 
    dplyr::filter(!pocketname %in% df_curation_pockets_exclude_list$pocketname)
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# FOCUS ON MOLECULES: HOW MANY BACTERIAL TARGETS DO THEY HIT?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# Calculate how many times each molecule hits a target (for bacteria-oil pair)
# 
df_counts_major <- df_targets_filtered %>%
  dplyr::filter(bacteria %in% config$major_bacteria) %>% 
  dplyr::group_by(oil_eng, bacteria, fullname) %>%
  dplyr::summarise(count = n(), .groups = "drop")

df_counts_minor <- df_targets_filtered %>%
  dplyr::filter(bacteria %in% config$minor_bacteria) %>% 
  dplyr::group_by(oil_eng, bacteria, fullname) %>%
  dplyr::summarise(count = n(), .groups = "drop")


# Extract a fixed set of molecule names across all datasets
# Note: we sort by descendant order by using the cumulative number of targets
# over all the bacteria
fixed_levels_major <- df_counts_major %>%
  dplyr::group_by(fullname) %>% 
  dplyr::summarise(count_t = sum(count), .groups = "drop") %>%
  dplyr::arrange(count_t) %>% 
  dplyr::pull(fullname)

fixed_levels_minor <- df_counts_minor %>%
  dplyr::group_by(fullname) %>% 
  dplyr::summarise(count_t = sum(count), .groups = "drop") %>%
  dplyr::arrange(count_t) %>% 
  dplyr::pull(fullname)


# 
# Barplots, one for each bacterium, in the same plot
# (separately for major and minor bacteria)
# 
combined_plot_mol_major <- ggplot(df_counts_major, 
                                  aes(x = count, y = factor(fullname, levels = fixed_levels_major))) +
  geom_col(fill = "blue") +
  facet_wrap(~ bacteria, ncol = 3) +  # Keeps all panels the same height
  theme_minimal() +
  labs(x = "# Targets", y = "Molecules") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.line = element_line(color = "black", linewidth = 1.2)
  ) +
  scale_y_discrete(limits = fixed_levels_major, drop = FALSE)  # Keeps all molecules in every plot


combined_plot_mol_minor <- ggplot(df_counts_minor, 
                                  aes(x = count, y = factor(fullname, levels = fixed_levels_minor))) +
  geom_col(fill = "blue") +
  facet_wrap(~ bacteria, ncol = 3) +  # Keeps all panels the same height
  theme_minimal() +
  labs(x = "# Targets", y = "Molecules") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.line = element_line(color = "black", linewidth = 1.2)
  ) +
  scale_y_discrete(limits = fixed_levels_minor, drop = FALSE)  # Keeps all molecules in every plot

# 
# Save the plots as an SVG files
# 
ggsave(file.path(output_dir, "barplot_mol_nrtar_majorbacteria.svg"), 
       plot = combined_plot_mol_major, 
       width = 12, height = 8, units = "in", device = "svg")

ggsave(file.path(output_dir, "barplot_mol_nrtar_minorbacteria.svg"), 
       plot = combined_plot_mol_minor, 
       width = 12, height = 8, units = "in", device = "svg")


#
# Phytocomplexes ---------------------------------------------------------------
# 
# FOCUS ON PHYTOCOMPLEXES: HOW DO THE PHYTOCOMPLEXES INTERACT WITH EACH TARGET?  
# First, we calculate a contribution for each target, by considering the 
# zzscore and the percentage (composition)
# 
df_targets_contributions <- df_targets_filtered %>%
  dplyr::select(bacteria, oil_eng, fullname, perc_comp, genesymbol, zzscore) %>%
  dplyr::mutate(perc_weight = ifelse(perc_comp>1, 1+log10(perc_comp), perc_comp)) %>%
  dplyr::mutate(contribution = zzscore * perc_weight)

# 
# Then, we calculate the sum of contributions over all the molecules of 
# phytocomplexes and add the data about centrality
# 
df_targets_phytocomplex_centr <- df_targets_contributions %>%
  dplyr::group_by(bacteria, genesymbol, oil_eng) %>%
  dplyr::summarise(sum_contr = sum(contribution)) %>% 
  dplyr::left_join(df_targets_contr_and_centrality, 
                   by=c("oil_eng", "bacteria", "genesymbol"))



# 
# We store the data of zzscore from molecule-target of oil-bacterium pairs
# 
df_bactdata_zzscores <- df_targets_contributions %>%
  dplyr::left_join(df_multiple, by=c("bacteria", "genesymbol")) %>%
  dplyr::mutate(gene=if_else(!is.na(genesymbol_group),
                             genesymbol_group, genesymbol)) %>%
  dplyr::group_by(bacteria, oil_eng, fullname, gene) %>%
  slice_max(order_by = zzscore, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(bacteria, oil_eng, fullname, zzscore, gene)


# 
# Extract top n genes for each bacteria/oil pair
# We can set the variable n_top_genes
# 
df_top_genes <- df_targets_phytocomplex_centr %>%
  group_by(oil_eng, bacteria) %>%
  slice_max(order_by = sum_contr, n = n_top_genes, with_ties = FALSE) %>%
  ungroup()


# Create a unique list of genes per bacteria (i.e., each bacteria's top n genes)
unique_genes_per_bacteria <- df_top_genes %>%
  group_by(bacteria) %>%
  distinct(genesymbol) %>%
  ungroup()


# Use the shortlist of genes (in dataframe format) to filter the previous data
df_filtered <- df_targets_phytocomplex_centr %>%
  dplyr::right_join(unique_genes_per_bacteria, by=c("bacteria", "genesymbol"))


# Here we filter data by applying some exclusion rules based on manual curation
# data (if PDB correspond to "exclude" or pocket correspond to "exclude")
# if (manual_curation) {
#   df_filtered <- df_filtered
# }

# Calculate the sum of `sum_contr` over the three oils for each gene
# This will be used to rank genes in the Y axis of the heatmaps and scatterplot
df_ranked_genes <- df_filtered %>%
  group_by(bacteria, genesymbol) %>%
  summarize(total_sum_contr = sum(sum_contr, na.rm = TRUE), .groups = "drop") %>%
  arrange(bacteria, desc(total_sum_contr))  # Rank genes by descending sum_contr


# Merge the ranked genes back with the filtered dataset (which is now ranked)
df_filtered_ranked <- df_filtered %>%
  left_join(df_ranked_genes, by = c("bacteria", "genesymbol"))


# We create a new dataframe in which we add a new column, that is named 'gene'
# This column is the name of grouped genes whenever available, otherwise it is
# the old genesymbol name
# Then we use this (column gene) with bacteria and oil_eng to group rows and 
# take only one row for group (with highest centrality)
df_filtered_groupedgenes <- df_filtered_ranked %>% 
  dplyr::left_join(df_multiple, by=c("bacteria", "genesymbol")) %>% 
  dplyr::mutate(gene=if_else(!is.na(genesymbol_group), 
                             genesymbol_group, genesymbol)) %>% 
  dplyr::group_by(bacteria, oil_eng, gene) %>% 
  slice_max(order_by = avg_centrality, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(defined_centrality = ifelse(is.na(avg_centrality),  
                                            "Undefined", "Defined")) %>% 
  dplyr::mutate(defined_centrality = as.factor(defined_centrality))


# 
# Loop over the bacteria
# 

for (bact in config$major_bacteria) {

  # 
  # extract from df_filtered_groupedgenes data of the selected bacterium
  # 
  df_bactdata_tar <- df_filtered_groupedgenes %>% 
    dplyr::filter(bacteria==bact) %>% 
    dplyr::select(oil=oil_eng, gene, defined_centrality,
                  centrality=avg_centrality,
                  target_overall=total_sum_contr,
                  target_zzscore_sum=sum_contr)
  
  # 
  # This is a dataframe where we have the selected targets for the given bacterium
  # 
  df_bactdata <- extract_target_data(df_bactdata_tar)
  
  # ------------------------------------------------------------------------------
  # Three steps (three functions) to create the matrix for each bacterium: 
  #   i) extract percent data; 
  #   ii) create empty matrix; 
  #   iii) fill the matrix with data
  # ------------------------------------------------------------------------------
  
  moldata_df1 <- extract_perc_values(df_comp=df_composition_new, oil="Cinnamon")
  moldata_df2 <- extract_perc_values(df_comp=df_composition_new, oil="Oregano")
  moldata_df3 <- extract_perc_values(df_comp=df_composition_new, oil="Thyme")
  
  matrix1 <- create_empty_matrix(df_mol=moldata_df1, df_tar=df_bactdata)
  matrix2 <- create_empty_matrix(df_mol=moldata_df2, df_tar=df_bactdata)
  matrix3 <- create_empty_matrix(df_mol=moldata_df3, df_tar=df_bactdata)
  
  # Note: df_bactdata_zzscores contains data for all the bacteria, 
  #       we filter inside the function 'fill_matrix_with_data'
  matrix1 <- fill_matrix_with_data(bact=bact, oil="Cinnamon", m=matrix1, 
                                   df=df_bactdata_zzscores)
  matrix2 <- fill_matrix_with_data(bact=bact, oil="Oregano", m=matrix2, 
                                   df=df_bactdata_zzscores)
  matrix3 <- fill_matrix_with_data(bact=bact, oil="Thyme", m=matrix3, 
                                   df=df_bactdata_zzscores)
  
  nr_targets <- nrow(matrix1)
  
  # 
  # Update the names of the proteins (according to manual curation file)
  # 
  for (n in 1:nrow(matrix1)) {
    orig_protein_name <- rownames(matrix1)[n]
    if (orig_protein_name %in% df_curation_proteins$genesymbol) {
      new_protein_name <- df_curation_proteins$protein[df_curation_proteins$genesymbol==orig_protein_name]
      rownames(matrix1)[n] <- new_protein_name
      if (config$debug) {
        cat("Orig protein name", orig_protein_name, "\tnew protein name", new_protein_name, "\n")
      }
    }
  }
  
  
  # 
  # Define the bottom annotations
  #
  column_ha1 <- HeatmapAnnotation(
    "Composition (%)" = anno_barplot(
      moldata_df1$perc_comp, 
      gp = gpar(fill = "orange"), ylim = c(0, 100),
      axis_param = list(side = "left", at = c(0, 50, 100), labels = c("0", "50", "100"))
    ), 
    height = unit(1, "cm"),
    show_annotation_name = c("Composition (%)" = TRUE),
    annotation_name_rot = c("Composition (%)" = 0),
    annotation_name_offset = c("Composition (%)" = "8mm"),
    annotation_name_side = "left"
  )
  column_ha2 <- HeatmapAnnotation(
    "Composition (%)" = anno_barplot(
      moldata_df2$perc_comp, 
      gp = gpar(fill = "orange"), ylim = c(0, 100),
      axis_param = list(side = "left", at = c(0, 50, 100), labels = c("0", "50", "100"))
    ), 
    height = unit(1, "cm"),
    show_annotation_name = c("Composition (%)" = FALSE),
    annotation_name_side = "left"
  )
  column_ha3 <- HeatmapAnnotation(
    "Composition (%)" = anno_barplot(
      moldata_df3$perc_comp, 
      gp = gpar(fill = "orange"), ylim = c(0, 100),
      axis_param = list(side = "left", at = c(0, 50, 100), labels = c("0", "50", "100"))
    ), 
    height = unit(1, "cm"),
    show_annotation_name = c("Composition (%)" = FALSE),
    annotation_name_side = "left"
  )
  
  
  # 
  # Define the top annotations
  #
  top_ha1 = HeatmapAnnotation(
    title_anno = anno_block(gp = gpar(fill = oil_colors["Cinnamon"]),
      labels = c("Cinnamon"), 
      labels_gp = gpar(col = "Black", fontface = "bold", fontsize = 12))
  )
  top_ha2 = HeatmapAnnotation(
    title_anno = anno_block(gp = gpar(fill = oil_colors["Oregano"]),
      labels = c("Oregano"), 
      labels_gp = gpar(col = "white", fontface = "bold", fontsize = 12))
  )
  top_ha3 = HeatmapAnnotation(
    title_anno = anno_block(gp = gpar(fill = oil_colors["Thyme"]),
      labels = c("Thyme"), 
      labels_gp = gpar(col = "Black", fontface = "bold", fontsize = 12))
  )
  
  
  # 
  # Parameters for the three heatmaps
  # 
  row_fontsize <- 8
  col_fontsize <- 6
  col_heatmap <- colorRamp2(c(0, 5), c("white", "red"))
  col_centrality <- colorRamp2(c(0, 1), c("white", "purple"))
  
  # 
  # Heatmaps for the three oils
  # 
  # ht1=Cinnamon
  # ht2=Oregano
  # ht3=Thyme
  # 
  ht1 <- Heatmap(
    matrix1,
    name = "ZZ-score",
    rect_gp = gpar(col = "#E5E4E2", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_fontsize),
    column_names_gp = gpar(fontsize = col_fontsize),
    column_names_rot = 60,
    col = col_heatmap,
    top_annotation = top_ha1,
    bottom_annotation = column_ha1,
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  ht2 <- Heatmap(
    matrix2,
    name = "ZZ-score",
    rect_gp = gpar(col = "#E5E4E2", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_fontsize),
    column_names_gp = gpar(fontsize = col_fontsize),
    column_names_rot = 60,
    col = col_heatmap,
    top_annotation = top_ha2,
    bottom_annotation = column_ha2,
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  ht3 <- Heatmap(
    matrix3,
    name = "ZZ-score",
    rect_gp = gpar(col = "#E5E4E2", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_fontsize),
    column_names_gp = gpar(fontsize = col_fontsize),
    column_names_rot = 60,
    col = col_heatmap,
    top_annotation = top_ha3,
    bottom_annotation = column_ha3,
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  
  # 
  # Define the centrality as annotation
  # 
  ha_centrality <- rowAnnotation(
    Centrality = df_bactdata$mean_centrality, 
    col = list(Centrality = col_centrality),
    na_col = "#E5E4E2",
    annotation_name_offset = c(Centrality = "3mm"),
    annotation_legend_param = list(Centrality = list(direction = "horizontal"))
  )

  
  # 
  # Define the scatterplot as annotation
  # 
  # First, we create an empty matrix with rows sorted as the other heatmaps
  # Then, we insert data
  matrix_scatterplot <- matrix(NA, nrow = nrow(df_bactdata), ncol = 3)
  colnames(matrix_scatterplot) <- c("Cinnamon", "Oregano", "Thyme")
  rownames(matrix_scatterplot) <- df_bactdata$gene
  
  for (r in seq(1:nrow(df_bactdata_tar))) {
    if (df_bactdata_tar$oil[r] %in% colnames(matrix_scatterplot) &
        df_bactdata_tar$gene[r] %in% rownames(matrix_scatterplot)) {
      matrix_scatterplot[ which(rownames(matrix_scatterplot)==df_bactdata_tar$gene[r]),
                          which(colnames(matrix_scatterplot)==df_bactdata_tar$oil[r]) ] <- df_bactdata_tar$target_zzscore_sum[r]
    }
  }
  
  
  # Finally, we create the scatterplot annotation
  ha_scatterplot = rowAnnotation(
    "WSumZZscore_refined" = anno_points(matrix_scatterplot, ylim = c(0, 40),
                                        pch = 19, gp = gpar(col = oil_colors)),
    width = unit(5, "cm")
  )
  
  # Define the heatmap list (including two row annotations)
  ht_list <- ht1 + ht2 + ht3 + ha_centrality + ha_scatterplot

  # Open SVG device
  output_file <- file.path(output_dir, sprintf("heatmaps_%s_best%d.svg", bact, n_top_genes))
  svg(output_file, width = 12, height = 1.5+nr_targets/5)
  
  draw(
    ht_list,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends = TRUE,
    ht_gap = unit(8, "mm"),
    padding = unit(c(5, 5, 15, 5), "mm")  # Add padding for space around the plot
  )
  
  # Add title
  grid.text(bact, x = 0.5, y = 0.97, gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Close device
  dev.off()
  
} 



# 
#     
# 
