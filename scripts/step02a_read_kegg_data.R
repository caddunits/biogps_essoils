# NOTE: If the script, for some reasons, crashes with the following message
#       Error in .getUrl(url, .flatFileParser) : Forbidden (HTTP 403)
#       adjust the content of this line of code: Sys.sleep(0.3)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# It is preferable to run the script step-by-step (for example, on Rstudio)
# 

library(KEGGREST)
library(dplyr)
library(magrittr)
library(purrr)

repo_dir <- getwd()

# Input
data_dir <- file.path(repo_dir, "data")

# Output
kegg_dir <- file.path(data_dir, "KEGG_data")

if (!file.exists(kegg_dir)){
  dir.create(kegg_dir)
}
# 
# For each KEGG code we should obtain a corresponding .RDS file
# 

verbose <- FALSE
# 
# NOTE: Change the row below according to the org_code which is in progress
# 

orgs <- c("adh","eab","ebd","ebe","ebl","ebr","ebw","ecc","ecd","ece","ecf","ecg",
"eci","ecj","eck","ecl","ecm","eco","ecoa","ecob","ecoc","ecoh","ecoi","ecoj",
"ecok","ecol","ecoo","ecos","ecp","ecq","ecr","ecs","ect","ecv","ecw","ecx",
"ecy","ecz","edj","efa","efd","efi","efl","efn","efq","efs","eih","ekf","eko",
"elc","eld","elf","elh","ell","eln","elo","elp","elr","elu","elw","elx","ena",
"ene","eoc","eoh","eoi","eoj","eok","ese","esl","esm","eso","etw","eum","eun",
"kpa","kpb","kpc","kpg","kph","kpi","kpj","kpm","kpn","kpne","kpnk","kpnu","kpo",
"kpp","kpq","kpr","kps","kpt","kpu","kpv","kpw","kpx","kpy","kpz","pae","paeb",
"paec","paeg","paei","pael","paem","paep","paer","paes","paeu","paev","paf","pag",
"pap","pau","pdk","pnc","prp","psg","saa","sab","sac","sad","sae","sah","sam",
"sams","sao","sar","sas","sau","saua","saub","sauc","saud","saue","sauf","saug",
"saui","sauj","sauk","saum","saun","sauq","saur","saus","saut","sauu","sauv",
"sauw","saux","sauy","sauz","sav","saw","sax","sech","sep","sepp","seps","ser",
"sud","sue","suf","sug","suj","suk","suq","sut","suu","suv","suw","sux","suy","suz")

for (org_code in orgs){

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
  
  # VERY IMPORTANT NOTE: In case the loop ends with an error, run it again
  # In some cases it may require several runs.
  for (k in wanted_idx) {
    cat("index:", k, "\n")
    
    k_info <- org_data[k]
    k_name <- names(org_data[k])
    
    # Pass further in case data are already in the df
    if (k_name %in% org_df$k_name) next()
    if (k_name %in% org_strange) next()
    if (k_name %in% org_unknown) next()
    
    if (k_info == "gene") {
      uniprot_code <- NA
      pathway_names <- NA
      pathway_values <- NA
      
    } else if (k_info == "CDS") {
      if (verbose) cat("\tK_NAME:", k_name, "\n")
      
      # 
      # Repeat the loop untill it safely reaches the end
      # 
      repeat{
         query <- safely(keggGet)(k_name)
         if(is.null(query$error))
           break
         else if(query$error$message != "Forbidden (HTTP 403).")
           stop(query$error$message)
         else
           Sys.sleep(5)
      }
      
      k_data <- query$result
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
}