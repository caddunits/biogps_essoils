# NOTE: If the script, for some reasons, crashes with the following message
#       Error in .getUrl(url, .flatFileParser) : Forbidden (HTTP 403)
#       adjust the content of this line of code: Sys.sleep(0.3)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# It is preferable to run the script step-by-step (for example, on Rstudio)
# 

library(KEGGREST)
library(dplyr)
library(magrittr)

repo_dir <- getwd()

# Input
data_dir <- file.path(repo_dir, "data")

# Output
kegg_dir <- file.path(data_dir, "KEGG_data")


# 
# For each KEGG code we should obtain a corresponding .RDS file
# 

verbose <- FALSE
# 
# NOTE: Change the row below according to the org_code which is in progress
# Last run was with this org_code 
# 
org_code <- "edh"


org_data <- keggList(org_code)

cat("Number of entries:", length(org_data), "\n")

org_df <- data.frame(
  k_name = character(0),
  k_entry = character(0),
  uniprot = character(0),
  NCBI = character(0),
  pw_name = character(0),
  pw_value = character(0),
  stringsAsFactors = FALSE
)

org_strange <- c()
org_unknown <- c()

# 
# This loop might crash, likely because of the number of connections, or else. 
# Slightly better when adding the line Sys.sleep
# 
# Anyway, in order to bypass the problem, the script checks the data already
# in the dataframe, so it is enough to re-run the loop several times
# until it correctly ends.
# 

wanted_idx <- seq(1:length(org_data))

# VERY IMPORTANT NOTE: In case the loop ends with an error, run again it.
# In some cases it may require several runs.
for (k in wanted_idx) {
  cat("index:", k, "\n")
  
  k_info <- org_data[k]
  k_name <- names(org_data[k])
  
  # Pass further in case data are already in the df
  if (k_name %in% org_df$k_name) next()
  if (k_name %in% org_strange) next()
  if (k_name %in% org_unknown) next()
  
  # 
  # Pay attention, this line is repeated below without the tryCatch
  # It is correct, and the use of tryCatch could lead to wrong assignments
  # in case of network problems due to the multiple request to the server
  # that after a while cause the server to fail in the response.
  # 
  
  Sys.sleep(0.3)
  k_data <- tryCatch(keggGet(k_name), error=function(e) NULL)
  k_entry <- k_data[[1]]$ENTRY
  
  if (k_info == "gene") {
    uniprot_code <- NA
    pathway_names <- NA
    pathway_values <- NA
    
  } else if (k_info == "CDS") {
    if (verbose) cat("\tK_NAME:", k_name, "\n")
    
    # See the comment above
    k_data <- keggGet(k_name)
    k_entry <- k_data[[1]]$ENTRY
    
    if (length(k_data) == 1) {
      
      if (!"PATHWAY" %in% names(k_data[[1]])) {
        
        uniprot_code <- NA
        NCBI_code <- NA
        pathway_names <- NA
        pathway_values <- NA
        
      } else {
        
        # Extract data for pathways
        pathway_names <- names(k_data[[1]]$PATHWAY)
        pathway_values <- k_data[[1]]$PATHWAY %>% unlist()
        
        if (verbose) {
          cat("CHECK PATHWAYS DATA...\n")
          print(pathway_names)
          print(pathway_values)
        }
        
        # Extract data for uniprot
        is_uniprot <- which(grepl("UniProt: ", k_data[[1]]$DBLINKS))
        #print(length(is_uniprot))
        
        is_NCBI <- which(grepl("NCBI-ProteinID: ", k_data[[1]]$DBLINKS))
        #print(length(is_NCBI))
        
        # Strange data (neither Uniprot nor NCBI)
        if (length(is_uniprot) == 0 & length(is_NCBI) == 0) {
          if (verbose) {
            cat("STRANGE UNIPROT & NCBI\n")
            print(k_data)
          }
          org_strange <- c(org_strange, k_name)
          next()
          
        } else {
          uniprot_code <- NA
          NCBI_code <- NA
          # Either Uniprot or NCBI is fine (or both)
          
          if (length(is_uniprot) == 1) {
            uniprot_code <- gsub(pattern="UniProt: ",
                                 replacement = "",
                                 x = k_data[[1]]$DBLINKS[is_uniprot])
          }
          if (length(is_NCBI) == 1) {
            NCBI_code <- gsub(pattern="NCBI-ProteinID: ",
                                 replacement = "",
                                 x = k_data[[1]]$DBLINKS[is_NCBI])
          } 
          if (is.null(pathway_names)) pathway_names <- NA
          if (is.null(pathway_values)) pathway_values <- NA
          
        }
        
      }
      
    } else {
      if (verbose) {
        cat("STRANGE LENGTH\n")
        print(k_data)
      }
      org_strange <- c(org_strange, k_name)
      next()
    }
    
  } else {
    if (verbose) {
      cat("\nUNKNOWN:\n")
      print(k_info)
    }
    org_unknown <- c(org_unknown, k_name)
    next()
  }
  
  # Finally create a new dataframe to bind to the data
  single_df <- data.frame(
    k_name = k_name,
    k_entry = k_entry,
    uniprot = uniprot_code,
    NCBI = NCBI_code,
    pw_name = pathway_names,
    pw_value = pathway_values,
    stringsAsFactors = FALSE
  )
  org_df <- rbind(org_df, single_df)
  
}


# To check the state of calculations (percentage)
success_rate <- length(unique(org_df$k_entry))/length(org_data)
print(success_rate)

row.names(org_df) <- NULL

# Save dato to a .RDS file
saveRDS(object = org_df, file=file.path(kegg_dir, paste0(org_code, ".RDS")))

