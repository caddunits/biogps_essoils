# 
# This script runs the analysis with alluvial plots 
# 
# Input are data from .RDS files:
#  sign_pathways.RDS
#  molecules_interaction.RDS
#  
# Output are alluvial plots (as .SVG)
# 
# Note: given that pathway names vary in length, their style is curated with 
#       specific rules to correctly split names in two lines when too long
#       

library(dplyr)
library(magrittr)
library(writexl)
library(ggplot2)
library(ggalluvial)


# 
# Functions
# 
# Assign colors based on the presence of the oils in the dataframe
# Note that the alphabetic order is guaranteed by the arrange function  
assign_flow_colors <- function(oils) {
  # Define line-colors based on the presence of different oils in the plot
  if (all(c("Cinnamon", "Oregano", "Thyme") %in% oils)) {
    scale_fill_manual_colors <- c("#B66E34", "#244A0D", "#95C977")
  } else if (all(c("Cinnamon", "Oregano") %in% oils)) {
    scale_fill_manual_colors <- c("#B66E34", "#244A0D")
  } else if (all(c("Cinnamon", "Thyme") %in% oils)) {
    scale_fill_manual_colors <- c("#B66E34", "#95C977")
  } else if (all(c("Oregano", "Thyme") %in% oils)) {
    scale_fill_manual_colors <- c("#244A0D", "#95C977")
  } else if (all(c("Oregano") %in% oils)) {
    scale_fill_manual_colors <- c("#244A0D")
  } else if (all(c("Thyme") %in% oils)) {
    scale_fill_manual_colors <- c("#95C977")
  } else if (all(c("Cinnamon") %in% oils)) {
    scale_fill_manual_colors <- c("#B66E34")
  } else {
    scale_fill_manual_colors <- c("yellow")
  }
  return(scale_fill_manual_colors)
}


# Change this line according to the path of your repository
repo_dir <- getwd()

# Settings are set once, in a configuration file
source(file.path(repo_dir, "scripts", "config.R"))


# Input
data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")

# 
# Manual curation to deal with genesymbol names
# 
fullpath_datacuration <- file.path(data_dir, "Curation_data", "manual_curation.xlsx")

# df_multiple <- readxl::read_xlsx(path=fullpath_datacuration,
#                                  sheet = "multiple_genes")

df_curation_targets <- readxl::read_xlsx(path=fullpath_datacuration,
                                          sheet = "proteins")

# if (orig_protein_name %in% df_curation_proteins$genesymbol) {
#   new_protein_name <- df_curation_proteins$protein[df_curation_proteins$genesymbol==orig_protein_name]
#   rownames(matrix1)[n] <- new_protein_name
  
# Output
output_dir <- file.path(repo_dir, "output", "gsea_alluvialplots")

# Create the ouput directory if it does not exist
dir.create(output_dir, showWarnings = FALSE)



# 
# First, read from .RDS files the data we need
# 

# df_gsea_pathways <- readRDS(file = file.path(mapp_dir, "sign_pathways.RDS"))
# df_gsea_molecules_interaction <- readRDS(file = file.path(mapp_dir, "molecules_interaction.RDS"))

# 
# modi 2025-06-09
# We read data from the XLSX file from gsea_analysis instead from .RDS
# so we can take different sheets (including the info about STRONG interactions) 
# 
pathwaysenrichment_gsea <- file.path(repo_dir, "output", 
                                     "gsea_pathways_enrichment",
                                     "Results_pathwaysenrichment_gsea.xlsx")


# df_gsea_pathways <- readRDS(file = file.path(mapp_dir, "sign_pathways.RDS"))
df_gsea_pathways <- readxl::read_xlsx(path=pathwaysenrichment_gsea,
                                 sheet = "Pathways_filtered")


# df_gsea_molecules_interaction <- readRDS(file = file.path(mapp_dir, 
#                                                           "molecules_interaction.RDS"))
df_gsea_molecules_interaction <- readxl::read_xlsx(path=pathwaysenrichment_gsea,
                                                    sheet = "Pathways_details")

# This is new
df_gsea_targets_interaction <- readxl::read_xlsx(path=pathwaysenrichment_gsea,
                                                    sheet = "Targets_interaction")

alplot_params <- list(
  size_oil = 5,
  size_mol = 4,
  size_gene = 4,
  size_pathway = 5,
  title_size = 18,
  size_xlabel = 14,
  threshold_pathname = 12
)


processed_pairs <- c()


# Loop over all the pathways
# for (p in 57:57) {
for (p in 1:nrow(df_gsea_pathways)) {
  
  wanted_bacteria <- df_gsea_pathways$bacteria[p]
  wanted_pathway <- df_gsea_pathways$pathway[p]
  cat("\n * * * ID:", p, "\nWorking on bacteria =", wanted_bacteria, 
      "pathway =", wanted_pathway, "\n")
  pair <- sprintf("%s_%s", wanted_bacteria, wanted_pathway)
  if (pair %in% processed_pairs) {
    cat("SKIP PAIR")
    next()
  }
  processed_pairs <- c(processed_pairs, pair)
  
  # Extract relevant data (genes with strong interactions)
  strong_targets <- df_gsea_targets_interaction %>% 
    dplyr::filter(bacteria == wanted_bacteria,
                  pathway == wanted_pathway) %>% 
    dplyr::filter(class_interaction == "STRONG") %>% 
    dplyr::select(genesymbol) %>% 
    dplyr::distinct()
  
  df_gsea_single_pathway <- df_gsea_molecules_interaction %>% 
    dplyr::filter(bacteria == wanted_bacteria,
                  pathway == wanted_pathway) %>% 
    dplyr::filter(genesymbol %in% strong_targets$genesymbol) %>% 
    dplyr::select(bacteria, pathway, oil, fullname, phytoclass, perc_comp, genesymbol) %>% 
    dplyr::rename(molecule = fullname) %>% 
    dplyr::distinct() %>% 
    dplyr::left_join(df_curation_targets, by=c("bacteria", "genesymbol"))
  
  missing_curation <- df_gsea_single_pathway$genesymbol[is.na(df_gsea_single_pathway$protein)]
  if (length(missing_curation)==0) {
    cat("All targets are manually curated")
    # print(unique(df_gsea_single_pathway$genesymbol))
  } else {
    cat("Please check and curate the name of the following targets:\n")
    print(unique(missing_curation))
  }
  
  # Define whether use one line or two lines for the pathway name
  nr_flows <- nrow(df_gsea_single_pathway)
  length_pathname <- nchar(wanted_pathway)
  
  if (length_pathname/nr_flows > alplot_params$threshold_pathname) {
    cat("FLOWS:", nr_flows, "PATH-LENGTH:", length_pathname, "\tRATIO", length_pathname/nr_flows, "\n")
    cat("NEED SPLIT THE PATHWAY NAME IN TWO PARTS\n")
    
    # Identify the positions of the empty spaces in the pathway name
    emptyspace_positions <- stringr::str_locate_all(string = wanted_pathway, pattern=" ")[[1]][,2]

    # We need to identify which is the closest to the middle of the name
    # NOTE: This procedure is valid for one replacement (so, two-lines names)
    #       but it might be adapted also to three-lines names
    length_half_pathway <- ceiling(length_pathname/2)
    length_possible_replacers <- length_pathname - emptyspace_positions
    distance_possible_replacers <- abs(length_possible_replacers - length_half_pathway)
    best_replacer_id <- which(distance_possible_replacers == min(distance_possible_replacers))
    best_replacer_position <- emptyspace_positions[best_replacer_id]
    
    # Split the name into character, replace the space and join characters again
    wanted_pathway_strsplit <- strsplit(wanted_pathway, "")[[1]]
    wanted_pathway_strsplit[best_replacer_position] <- "\n"
    twolined_pathwayname <- paste(wanted_pathway_strsplit, collapse = "")
    
    # Finally insert the new name in the dataframe
    df_gsea_single_pathway$pathway <- twolined_pathwayname
  }
  
  
  # Use more lines for the genesymbol name when there is the symbol / 
  genesymbol_newline <- gsub("/", "\n", df_gsea_single_pathway$genesymbol)
  df_gsea_single_pathway$genesymbol <- genesymbol_newline

  assigned_colors <- assign_flow_colors(oils=df_gsea_single_pathway$oil)
  
  # Stratified data is created based on the dataframe df_gsea_single_pathway
  # Two additional columns are created: 
  # object contains all the objects of the plot
  # categories contains the category each object belongs to
  strata <- rbind(
    df_gsea_single_pathway %>% 
      dplyr::select(object=oil) %>% 
      dplyr::mutate(category="oil", additional=NA) %>% 
      dplyr::arrange(dplyr::desc(object)),
    
    df_gsea_single_pathway %>% 
      dplyr::select(object=molecule, additional=phytoclass) %>%
      dplyr::mutate(category="molecule") %>% 
      dplyr::arrange(dplyr::desc(object)),
    
    df_gsea_single_pathway %>% 
      dplyr::select(object=protein) %>%
      dplyr::mutate(category="genesymbol", additional=NA) %>%
      # dplyr::select(object=genesymbol) %>%
      # dplyr::mutate(category="genesymbol", additional=NA) %>%
      dplyr::arrange(dplyr::desc(object)),
    
    df_gsea_single_pathway %>% 
      dplyr::select(object=pathway) %>%
      dplyr::mutate(category="pathway", additional=NA) %>%
      dplyr::arrange(dplyr::desc(object))
    
  ) %>% 
    dplyr::distinct()
  
  strata_width <- 10
  strata_height <- 0.8 * nrow(strata)
  
  strata_addedinfo <- strata %>% mutate(
    width = case_when(
      category == "oil" ~ 4/12,
      category == "molecule" ~ 6/12,
      category == "genesymbol" ~ 4/12,
      category == "pathway" ~ 6/12
    ),
    angle = case_when(
      category == "oil" ~ 90,
      category == "molecule" ~ 0,
      category == "genesymbol" ~ 0,
      category == "pathway" ~ 90
    ),
    textcolor = case_when(
      object == "Cinnamon"  ~ "black",
      object == "Oregano"  ~ "white",
      object == "Thyme"  ~ "black",
      category == "molecule" ~ "darkgreen",
      category == "genesymbol" ~ "darkblue",
      category == "pathway" ~ "white"
    ),
    textsize = case_when(
      category == "oil" ~ alplot_params$size_oil,
      category == "molecule" ~ alplot_params$size_mol,
      category == "genesymbol" ~ alplot_params$size_gene,
      category == "pathway" ~ alplot_params$size_pathway
    ),
    fontface = case_when(
      category == "oil" ~ "bold",
      category == "molecule" ~ "italic",
      category == "genesymbol" ~ "plain",
      category == "pathway" ~ "bold"
    ),
    fill = case_when(
      object == "Cinnamon"  ~ "#B66E34",
      object == "Oregano"  ~ "#244A0D",
      object == "Thyme"  ~ "#95C977",
      # We could chenge color according to phytoclasses
      # additional == "Oxygenated monoterpenoids" ~ colors_phytoclasses$oxygenated_monoterpenoids,
      # additional == "Monoterpenoids" ~ colors_phytoclasses$monoterpenoids,
      # additional == "Sesquiterpenes" ~ colors_phytoclasses$sesquiterpenes,
      # additional == "Phenylpropanoids" ~ colors_phytoclasses$phenylpropanoids,
      # additional == "Other organic compounds" ~ colors_phytoclasses$other_organic_compounds,
      category == "molecule" ~ "beige",
      category == "genesymbol" ~ "#AEE7F8",
      category == "pathway" ~ "blue"
    )
  ) %>% 
    dplyr::distinct()
  
  alplot <- ggplot(data = df_gsea_single_pathway,
                   aes(axis1 = oil,
                       axis2 = molecule,
                       axis3 = protein,
                       # axis3 = genesymbol,
                       axis4 = pathway)) +
    geom_alluvium(aes(fill = oil), show.legend = FALSE) +
    scale_fill_manual(name = "", values = assigned_colors) +
    geom_stratum(width = strata_addedinfo$width,
                 fill = strata_addedinfo$fill) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),
              color = strata_addedinfo$textcolor,
              fontface = strata_addedinfo$fontface,
              size = strata_addedinfo$textsize,
              angle = strata_addedinfo$angle) +
    theme_minimal() +
    scale_x_discrete(limits = c("Oil", "Molecule", "Target", "Pathway"),
    # scale_x_discrete(limits = c("Oil", "Molecule", "Gene", "Pathway"),
                     expand = c(0, 0)) + 
    ggtitle(sprintf("%s: Composition and Targets", wanted_bacteria),
            subtitle = sprintf("Pathway = %s", wanted_pathway)) +
    theme(legend.position = "none",
          plot.title = element_text(size = alplot_params$title_size, face = "bold"),
          plot.subtitle = element_text(size = alplot_params$title_size - 2, face = "italic"),
          plot.background = element_rect(fill="white"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = alplot_params$size_xlabel, face = "bold")
    )
  
  # alplot_code <- sprintf("pathway%03d", p)
  # alplot_filename <- file.path(output_dir, paste0("Alluvial_", alplot_code, ".svg", collapse = ""))
  
  alplot_filename <- sprintf("%s %s.svg", wanted_bacteria, wanted_pathway)
  alplot_fullpath <- file.path(output_dir, alplot_filename)
  ggsave(filename = alplot_fullpath, device = "svg", 
         plot = alplot, width = strata_width, height = strata_height)
  
}

