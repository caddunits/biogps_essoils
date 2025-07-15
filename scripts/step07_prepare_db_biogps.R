library(DBI)
library(RSQLite)
library(dplyr)
library(magrittr)
library(tibble)

library(BBmisc)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("bio3d")
# BiocManager::install("reactome.db")
# BiocManager::install("biomaRt")
library(bio3d)
library(reactome.db)
library(biomaRt)


# 
# Functions --------------------------------------------------------------------
# 

# This function read only one line, with the release date of the pdb2uniprot
pdbtosp_read_releasedate <- function(pdbtosp_filename) {
  tmp_info <- data.table::fread(
    file=pdbtosp_filename, header=FALSE, sep="\t", skip=9, nrows = 1,
    col.names=c("origdata"), strip.white=FALSE, blank.lines.skip=TRUE
  ) %>% as.character()
  info <- stringr::str_replace(string=tmp_info, 
                               pattern="Release: +", replacement="")
  return(info)
}


# This function reads the information of the pdb - uniprot conversion
pdbtosp_read_content <- function(pdbtosp_filename) {
  
  # skip_lines_content <- "^Distributed|^Copyrighted|^-----"
  # skip_lines_content <- skip_lines_content+"|^UniProt - Swiss-Prot"
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


#
# Functions to read data using biomaRt
#
get_uniprot_gene_relationships_from_mart <- function(biomart, dataset) {
  # Create biomart Mart (needs internet connection)
  mart <- biomaRt::useMart(biomart = biomart, dataset = dataset)
  
  # Extract info for wanted uniprot codes (it needs internet connection)
  df <- biomaRt::getBM(
    attributes = c("uniprotswissprot", "entrezgene_id"), 
    bmHeader = T, 
    mart = mart
  ) %>% 
    dplyr::rename("uniprot"="UniProtKB/Swiss-Prot ID") %>% 
    dplyr::mutate(entrezgene_id = as.character(`NCBI gene (formerly Entrezgene) ID`)) %>% 
    dplyr::filter(!is.na(entrezgene_id)) %>% 
    dplyr::filter(uniprot != "") %>% 
    dplyr::select(entrezgene_id, uniprot)
  return(df)  
}


# Reads pathways for all organisms
get_reactome_pathways <- function() {
  # Import all pathways from REACTOME DB
  df <- as.data.frame(reactome.db::reactomePATHID2NAME) %>% 
    dplyr::rename(path_code = DB_ID, path_desc = path_name) %>% 
    dplyr::select(path_code, path_desc)
  
  df2 <- df %>% 
    dplyr::mutate(path_orgcode = substr(path_code, 3, 5)) %>% 
    dplyr::mutate(path_organism = dplyr::case_when(
      path_orgcode == "BTA" ~ "Bos taurus",
      path_orgcode == "CEL" ~ "Caenorhabditis elegans",
      path_orgcode == "CFA" ~ "Canis familiaris",
      path_orgcode == "DRE" ~ "Danio rerio",
      path_orgcode == "DDI" ~ "Dictyostelium discoideum",
      path_orgcode == "DME" ~ "Drosophila melanogaster",
      path_orgcode == "GGA" ~ "Gallus gallus",
      path_orgcode == "HSA" ~ "Homo sapiens",
      path_orgcode == "MMU" ~ "Mus musculus",
      path_orgcode == "MTU" ~ "Mycobacterium tuberculosis",
      path_orgcode == "PFA" ~ "Plasmodium falciparum",
      path_orgcode == "RNO" ~ "Rattus norvegicus",
      path_orgcode == "SCE" ~ "Saccharomyces cerevisiae",
      path_orgcode == "SPO" ~ "Schizosaccharomyces pombe",
      path_orgcode == "SSC" ~ "Sus scrofa",
      path_orgcode == "XTR" ~ "Xenopus tropicalis",
      .default = "UNDEFINED"
    )) %>% 
    dplyr::mutate(path_description = dplyr::case_when(
      grepl("Bos taurus:", path_desc) ~ gsub("Bos taurus: ", "", path_desc),
      grepl("Caenorhabditis elegans:", path_desc) ~ gsub("Caenorhabditis elegans: ", "", path_desc),
      grepl("Canis familiaris:", path_desc) ~ gsub("Canis familiaris: ", "", path_desc),
      grepl("Danio rerio:", path_desc) ~ gsub("Danio rerio: ", "", path_desc),
      grepl("Dictyostelium discoideum:", path_desc) ~ gsub("Dictyostelium discoideum: ", "", path_desc),
      grepl("Drosophila melanogaster:", path_desc) ~ gsub("Drosophila melanogaster: ", "", path_desc),
      grepl("Gallus gallus:", path_desc) ~ gsub("Gallus gallus: ", "", path_desc),
      grepl("Homo sapiens:", path_desc) ~ gsub("Homo sapiens: ", "", path_desc),
      grepl("Mus musculus:", path_desc) ~ gsub("Mus musculus: ", "", path_desc),
      grepl("Mycobacterium tuberculosis:", path_desc) ~ gsub("Mycobacterium tuberculosis: ", "", path_desc),
      grepl("Plasmodium falciparum:", path_desc) ~ gsub("Plasmodium falciparum: ", "", path_desc),
      grepl("Rattus norvegicus:", path_desc) ~ gsub("Rattus norvegicus: ", "", path_desc),
      grepl("Saccharomyces cerevisiae:", path_desc) ~ gsub("Saccharomyces cerevisiae: ", "", path_desc),
      grepl("Schizosaccharomyces pombe:", path_desc) ~ gsub("Schizosaccharomyces pombe: ", "", path_desc),
      grepl("Sus scrofa:", path_desc) ~ gsub("Sus scrofa: ", "", path_desc),
      grepl("Xenopus tropicalis:", path_desc) ~ gsub("Xenopus tropicalis: ", "", path_desc),
      .default = path_desc
    )) %>% 
    dplyr::select(path_code, path_organism, path_description)
  
  return(df2)
}


# Reads id for entrezgenes
get_reactome_entrezgenes <- function() {
  ls <- as.list(reactome.db::reactomeEXTID2PATHID)
  return(names(ls))
}


# Reads relationships between genes and pathways
get_reactome_relationships <- function() {
  
  reactome_data_aslist <- as.list(reactome.db::reactomeEXTID2PATHID)
  # reactomeEXTID2PATHID 
  # An annotation data object that maps Entrez Gene identifiers to Reactome 
  # pathway identifiers.
  # R object containing key and value pairs. Keys are Entrez Gene identifiers 
  # and values are the corresponding Reactome pathway identifiers.
  # Values are vectors of length 1 or greater depending on whether a given 
  # external identifier can be mapped to only one or more Reactome pathway
  # identifiers.
  
  list_idx <- seq(1:length(reactome_data_aslist))
  df <- do.call('rbind', lapply(list_idx, function(r) {
    data.frame(
      entrezgene_id = as.character(names(reactome_data_aslist)[r]), 
      pathway = reactome_data_aslist[[r]],
      stringsAsFactors = FALSE
    )
  }))
  return(df)
}


get_uniprot_gene_relationships_from_mart <- function(biomart, dataset) {
  # Create biomart Mart (needs internet connection)
  mart <- biomaRt::useMart(biomart = biomart, dataset = dataset)
  
  # Extract info for wanted uniprot codes (it needs internet connection)
  df <- biomaRt::getBM(
    attributes = c("uniprotswissprot", "entrezgene_id"), 
    bmHeader = T, 
    mart = mart
  ) %>% 
    dplyr::rename("uniprot"="UniProtKB/Swiss-Prot ID") %>% 
    dplyr::mutate(entrezgene_id = as.character(`NCBI gene (formerly Entrezgene) ID`)) %>% 
    dplyr::filter(!is.na(entrezgene_id)) %>% 
    dplyr::filter(uniprot != "") %>% 
    dplyr::select(entrezgene_id, uniprot)
  return(df)  
}


read_GS_file <- function(gsfile) {
  
  tmp1_df <- read.delim2(gsfile, header = TRUE) %>%
    dplyr::mutate(mol = gsub(" +$", "", Molecule_Name)) %>%
    dplyr::select(-c("InChIKey", "X", "Molecule_Name"))
  
  tmp2_df <- reshape2::melt(tmp1_df, value.name = "value_str", id = "mol") %>%
    dplyr::mutate(value_num=as.numeric(value_str)) %>%
    dplyr::select(-c("value_str"))
  
  colnames(tmp2_df) <- c("molname", "poc", "value_num")
  
  tmp3_df <- tmp2_df %>%
    dplyr::mutate(pocketname = stringr::str_replace(poc, "^X", ""))
  
  tmp4_df <- tmp3_df %>% dplyr::select(molname, pocketname, value_num) %>% 
    dplyr::mutate(pocketname = substr(pocketname, 1, 9))
  
  return(tmp4_df)
  
}


# This creates a dataframe with info for the molecules,
# by reading a dataframe extracted from .txt file from function read_GS_file
extract_molecules <- function(scores_df) {
  
  # Check the input df
  if (!"molname" %in% colnames(scores_df)) return(NULL)
  
  molecules_df <- scores_df %>% 
    dplyr::select(molname) %>% 
    dplyr::distinct(molname) %>% 
    dplyr::arrange(molname)
  
  return(molecules_df)
  
}

# This creates a dataframe with info for the pockets,
# by reading a dataframe extracted from .txt file from function read_GS_file
# Attention: in some examples chain has value different from ______________
extract_pockets <- function(scores_df) {
  
  # Check the input df
  if (!"pocketname" %in% colnames(scores_df)) return(NULL)
  
  pockets_df <- scores_df %>% 
    dplyr::select(pocketname) %>% 
    dplyr::distinct(pocketname) %>% 
    dplyr::arrange(pocketname) %>%
    dplyr::mutate(
      pdb = toupper(substr(pocketname, 1, 4)),
      intid = substr(pocketname, 7, 9),
      chain = gsub("_$", "", substr(pocketname, 11, 100))
    ) %>% 
    dplyr::mutate(pocketname = substr(pocketname, 1, 9))
  
  return(pockets_df)
  
}


# This access the PDB data through the package bio3d, and extract information
# PAY ATTENTION: This function needs internet connection to download annotations
get_pdb_annotation <- function(wanted_pdb_idx) {
  
  # Requested annotations
  anno_terms <- c("structureId", "experimentalTechnique", "resolution")
  
  tryCatch(
    {
      annotation_df <- bio3d::pdb.annotate(wanted_pdb_idx, anno.terms=anno_terms) %>% 
        dplyr::distinct() %>% 
        dplyr::rename(pdb=structureId) %>% 
        dplyr::mutate(
          exp_res = round(resolution, digits = 2), 
          exp_method = dplyr::case_when(
            experimentalTechnique == "X-ray" ~ "XRAY",
            experimentalTechnique == "NMR" ~ "NMR",
            experimentalTechnique == "EM" ~ "EM",
            experimentalTechnique == "IR" ~ "IR",
            experimentalTechnique == "Fiber" ~ "FIBER",
            experimentalTechnique == "Other" ~ "OTHER",
            experimentalTechnique == "Neutron" ~ "NEUTRON"
          )) %>% 
        dplyr::select(pdb, exp_method, exp_res)
      return(annotation_df)
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
      # },
      # warning=function(w) {
      # message('A Warning Occurred')
      # print(w)
      # return(NA)
    }
  )
  
}


# This add data to the table pdb (only for data not present)
eventually_add_pdb <- function(db_con, single_pockets_df) {
  
  # Extract again data for the PDB indeces in the table
  index_pdb_df <- dbGetQuery(db_con, "SELECT id AS id_pdb, pdb FROM pdb;")
  last_pdb_id <- max(as.numeric(index_pdb_df$id_pdb))
  
  # If a pdb is not present in the database, we should add it
  # Thus, extract the ids and the information from the web using package bio3d
  wanted_pdb_ids <- single_pockets_df %>% 
    dplyr::select(pdb) %>% 
    dplyr::distinct(pdb) %>% 
    dplyr::filter(!pdb %in% index_pdb_df$pdb)
  
  # In case there are some data to be added
  if (length(wanted_pdb_ids$pdb) > 0) {
    
    cat("\tNEW PDB DATA FOR ", nrow(wanted_pdb_ids), "PDB\n")
    
    # Access PDB to extract annotation data
    single_annotation_df <- get_pdb_annotation(wanted_pdb_idx = wanted_pdb_ids$pdb)
    # cat("XXXX\n")
    # print(single_annotation_df)
    # Check if the results is an error 
    # (catched through TryCatch, error caused by bio3d::pdb.annotate)
    # 
    if(!BBmisc::is.error(single_annotation_df)) {
      
      if (nrow(single_annotation_df) > 0) {
        
        # INSERT IN DATABASE DATA FOR PDB ENTRIES
        pdb4db_new <- single_annotation_df %>% 
          dplyr::mutate(id = seq(last_pdb_id+1, last_pdb_id+nrow(single_annotation_df))) %>% 
          dplyr::select(id, dplyr::everything())
        
        cat("\tDATA READY FOR DATABASE:pdb", nrow(pdb4db_new), "\n")
        
        # INSERT IN DATABASE DATA FOR POCKETS
        # Note: the table pdb already exists, thus it is with append=TRUE 
        dbWriteTable(db_con, "pdb", pdb4db_new, append=TRUE)
        
      } else {
        
        cat("\tNO NEW DATA SUITABLE FOR DATABASE:pdb\n")
        
      }
      
    } else {
      
      cat("\tNO NEW DATA SUITABLE FOR DATABASE:pdb - ERROR ON bio3d::pdb.annotate\n")
      
    }
    
  } else {
    
    cat("\tNO NEW DATA SUITABLE FOR DATABASE:pdb\n")
    
  }
  
  return()
  
}


# This add data to the table pockets (only for data not present)
eventually_add_pocket <- function(db_con, single_pockets_df) {
  
  # Extract again data for the PDB indeces in the table
  index_pdb_df <- dbGetQuery(db_con, "SELECT id AS id_pdb, pdb FROM pdb;")
  
  # Extract the ids for pockets. If no pockets (empty table), use 0
  index_pockets_df <- dbGetQuery(db_con, "SELECT id AS id_poc, name FROM pockets;")
  # cat("???\n")
  # print(nrow(index_pockets_df))
  # print(head(index_pockets_df$id_poc))
  # print(max(as.numeric(index_pockets_df$id_poc)))
  
  if (nrow(index_pockets_df) > 0) {
    last_pocket_id <- max(as.numeric(index_pockets_df$id_poc))
  } else {
    last_pocket_id <- 0
  }
  
  wanted_pocket_ids <- single_pockets_df %>% 
    dplyr::select(pocketname) %>% 
    dplyr::distinct(pocketname) %>% 
    dplyr::filter(!pocketname %in% index_pockets_df$name)
  
  wanted_pockets_df <- single_pockets_df %>% 
    dplyr::filter(pocketname %in% wanted_pocket_ids$pocketname) 
  
  # In case there are some data to be added
  if (nrow(wanted_pockets_df) > 0) {
    
    cat("\tNEW POCKETS DATA FOR ", nrow(wanted_pockets_df), "POCKETS\n")
    
    # These indeces are written in the table of the pockets
    pockets4db <- wanted_pockets_df %>%
      dplyr::left_join(index_pdb_df, by="pdb") %>%
      dplyr::select(-pdb) %>%
      dplyr::mutate(id = seq(last_pocket_id+1, last_pocket_id+nrow(wanted_pockets_df))) %>% 
      dplyr::select(id, dplyr::everything()) %>% 
      dplyr::rename(name=pocketname)
    
    cat("\tDATA READY FOR DATABASE:pockets", nrow(pockets4db), "\n")
    # print(head(pockets4db))
    
    dbWriteTable(db_con, "pockets", pockets4db, append=TRUE)
    
  } else {
    
    cat("\tNO NEW DATA SUITABLE FOR DATABASE:pockets\n")
    
  }
  
  return()
  
}


# This add data to the table molecules (only for data not present)
eventually_add_molecule <- function(db_con, scores_df) {
  
  # Create a dataframe with info for the molecules
  single_molecules_df <- extract_molecules(scores_df = scores_df)
  # print(colnames(single_molecules_df))
  
  # Extract the ids for moleculess. If no molecules (empty table), use 0
  index_molecules_df <- dbGetQuery(db_con, "SELECT id AS id_mol, code AS molname FROM molecules;")
  if (nrow(index_molecules_df) > 0) {
    last_mol_id <- max(as.numeric(index_molecules_df$id_mol))
  } else {
    last_mol_id <- 0
  }
  # print(last_mol_id)
  wanted_molecules_ids <- single_molecules_df %>% 
    dplyr::select(molname) %>% 
    dplyr::distinct(molname) %>% 
    dplyr::filter(!molname %in% index_molecules_df$molname)
  
  wanted_molecules_df <- single_molecules_df %>% 
    dplyr::filter(molname %in% wanted_molecules_ids$molname)
  # print(wanted_molecules_df)
  # In case there are some data to be added
  if (nrow(wanted_molecules_df) > 0) {
    
    cat("\tNEW POCKETS DATA FOR ", nrow(wanted_molecules_df), "MOLECULES\n")
    
    # cat("\tADD TO DATABASE", nrow(wanted_molecules_df), "MOLECULES\n")
    
    # These indeces are written in the table of the pockets
    molecules4db <- wanted_molecules_df %>%
      dplyr::mutate(id = seq(last_mol_id+1, last_mol_id+nrow(wanted_molecules_df))) %>% 
      dplyr::select(id, dplyr::everything()) %>% 
      dplyr::mutate(code=as.character(molname),
                    name=as.character(molname)) %>% 
      dplyr::select(-molname)
    
    cat("\tDATA READY FOR DATABASE:molecules", nrow(molecules4db), "\n")
    
    dbWriteTable(db_con, "molecules", molecules4db, append=TRUE)
    
  } else {
    
    cat("\tNO NEW DATA SUITABLE FOR DATABASE:molecules\n")
    
  }
  
  return()
  
}


# This adds the values of the GS score
add_score_values <- function(db_con, scores_df) {
  
  # Finally, add numeric values
  index_scores_df <- dbGetQuery(db_con, "SELECT id AS id_scorevalue FROM scores;")
  if (nrow(index_scores_df) > 0) {
    last_score_id <- max(as.numeric(index_scores_df$id_scorevalue))
  } else {
    last_score_id <- 0
  }
  
  index_molecule_df <- dbGetQuery(db_con, "SELECT id AS id_mol, code AS molname FROM molecules;")
  index_pocket_df <- dbGetQuery(db_con, "SELECT id AS id_poc, name AS pocketname FROM pockets;")
  
  # print(colnames(scores_df))
  # print(head(scores_df))
  
  # Read the (GS) score values and change names with id, 
  # then add a the identifier for the table
  scores4db <- scores_df %>% 
    dplyr::left_join(index_pocket_df, by="pocketname") %>% 
    dplyr::left_join(index_molecule_df, by="molname") %>% 
    dplyr::mutate(id = seq(last_score_id+1, last_score_id+nrow(scores_df))) %>% 
    dplyr::select(-molname, -pocketname) %>% 
    dplyr::rename(score = value_num)
  
  cat("\tDATA READY FOR DATABASE:scores", nrow(scores4db), "\n")
  
  # Note: these are only GS values
  dbWriteTable(db_con, "scores", scores4db, append=TRUE)
  
  return()
}


calculate_zscore_pockets <- function(db_con, wanted_refcode, wanted_classcode, poc_n_threshold) {
  query <- sprintf(("SELECT poc.id AS id_poc
, mol.id AS id_mol
, poc.name AS pocname
, sco.score
, mol.code
, mol.name AS molname
, mcl.molclass
, '%s' AS poc_refcode
, stp.poc_n
, stp.poc_mean
, stp.poc_sd
, (sco.score-stp.poc_mean)/stp.poc_sd AS zscore_poc
FROM scores sco
LEFT JOIN molecules mol ON mol.id = sco.id_mol
LEFT JOIN molclass_assignment mca ON mca.id_mol = mol.id
LEFT JOIN molclasses mcl ON mcl.id = mca.id_class
LEFT JOIN pockets poc ON poc.id = sco.id_poc
LEFT JOIN stats_poc stp ON stp.id_poc = poc.id AND stp.poc_refcode == '%s'
WHERE score > 0 AND poc_n > %d
AND molclass == '%s'"), wanted_refcode, wanted_refcode, poc_n_threshold, wanted_classcode)
  
  df <- dbGetQuery(db_con, query)
  return(df)
}



# 
# Define directories -----------------------------------------------------------
# 

repo_dir <- getwd()
data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")

#
# Data PDB<->Uniprot -----------------------------------------------------------
#


url <- "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pdbtosp.txt"
destfile <- file.path(mapp_dir, "pdbtosp.txt")

# Download in case file does not exist
if (!file.exists(destfile)) download.file(url, destfile, method = "curl")

# Read the release date (info to be stored in the database)
pdbtosp_releasedate <- pdbtosp_read_releasedate(destfile)

# Read the file (the procedure takes some minutes)
pdbtosp_df <- pdbtosp_read_content(destfile)

selected_methods <- c("XRAY", "NMR", "EM", "NEUTRON", "IR", "FIBER", "OTHER")
# Use this notation if want only one method: selected_methods <- c("FIBER")

pdb2uniprot_filename <- file.path(mapp_dir, "pdb2uniprot.txt")
# write.csv2(x=pdb2uniprot_df, file=pdb2uniprot_filename,row.names=FALSE)
if (!file.exists(pdb2uniprot_filename)) {
  pdb2uniprot_df <- pdbtosp_extract_content(df=pdbtosp_df, methods=selected_methods)
} else {
  pdb2uniprot_df <- read.csv2(file=pdb2uniprot_filename)
}




#
# DB connection ----------------------------------------------------------------
#

db_filename <- file.path(repo_dir, "db", "biogps.db")
db_con <- dbConnect(RSQLite::SQLite(), db_filename)



#
# Tables creation --------------------------------------------------------------
#
# Queries for tables creation
#

query_create_organisms <- "
CREATE TABLE organisms (
  id INTEGER PRIMARY KEY,
  organism TEXT NOT NULL
)"


query_create_genesymbols <- "
CREATE TABLE genesymbols (
  id INTEGER PRIMARY KEY,
  genesymbol TEXT NOT NULL
)"


query_create_uniprot <- "
CREATE TABLE uniprot (
  id INTEGER PRIMARY KEY,
  uniprot TEXT NOT NULL,
  gene TEXT NOT NULL,
  id_organism INTEGER,
  id_genesymbol INTEGER,
  FOREIGN KEY (id_organism) 
    REFERENCES organisms (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION,
  FOREIGN KEY (id_genesymbol) 
    REFERENCES genesymbols (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


query_create_pdb <- "
CREATE TABLE pdb (
  id INTEGER PRIMARY KEY,
  pdb TEXT NOT NULL,
  exp_method TEXT,
  exp_res REAL
)"


query_create_mapping_pdb_uniprot <- "
CREATE TABLE mapping_pdb_uniprot (
  id INTEGER PRIMARY KEY,
  id_pdb INTEGER,
  id_uniprot INTEGER,
  FOREIGN KEY (id_pdb) 
    REFERENCES pdb (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION,
  FOREIGN KEY (id_uniprot) 
    REFERENCES uniprot (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


# 
# Create a table for details of calculations, which will be filled in later
# Note that the details for the release date of pdb-uniprot mapping are set to
# default based on the file that is read here
# 
query_create_calcdetails <- sprintf("
CREATE TABLE calcdetails (
  projname TEXT,
  probes TEXT,
  protonation TEXT,
  pdbuniprot TEXT DEFAULT \"%s\"
)", pdbtosp_releasedate)



dbExecute(db_con, query_create_organisms)
dbExecute(db_con, query_create_genesymbols)
dbExecute(db_con, query_create_uniprot)
dbExecute(db_con, query_create_pdb)
dbExecute(db_con, query_create_mapping_pdb_uniprot)
dbExecute(db_con, query_create_calcdetails)


# This should be an empty dataframe
remaining_df <- pdbtosp_df %>% 
  dplyr::anti_join(pdb2uniprot_df, by="origdata")
cat("CHECK: in case all methods have been selected, this would be 0:", 
    nrow(remaining_df), "\n")


# Final dataframe with all the information
mapping_pdb2uniprot_df <- pdb2uniprot_df %>% 
  tidyr::separate_longer_delim(gene_uniprot, delim = ",") %>% 
  dplyr::mutate(info_split = gsub("\\(|\\)", "", gene_uniprot)) %>% 
  tidyr::separate_wider_delim(cols=info_split, delim=" ", 
                              names=c("gene", "uniprot")) %>% 
  tidyr::separate_wider_delim(cols=gene, delim="_", 
                              names=c("genesymbol", "organism"),cols_remove = FALSE) %>% 
  dplyr::select(-gene_uniprot) %>% 
  dplyr::rename(pdb=pdbid)


# 
# Prepare data to be inserted into database
#  

pdb4db <- mapping_pdb2uniprot_df %>% 
  dplyr::select(pdb, exp_method, exp_res) %>% 
  dplyr::arrange(pdb) %>% 
  distinct() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename(id=rowname)


organisms4db <- mapping_pdb2uniprot_df %>% 
  dplyr::select(organism) %>% 
  dplyr::arrange(organism) %>% 
  distinct() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename(id=rowname)


genesymbols4db <- mapping_pdb2uniprot_df %>% 
  dplyr::select(genesymbol) %>% 
  dplyr::arrange(genesymbol) %>% 
  distinct() %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename(id=rowname)


uniprot4db <- mapping_pdb2uniprot_df %>% 
  dplyr::select(uniprot, gene, genesymbol, organism) %>% 
  dplyr::arrange(uniprot) %>% 
  distinct() %>% 
  dplyr::left_join(organisms4db, by="organism") %>% 
  dplyr::rename(id_organism = id) %>% 
  dplyr::select(-organism) %>% 
  dplyr::left_join(genesymbols4db, by="genesymbol") %>% 
  dplyr::rename(id_genesymbol = id) %>% 
  dplyr::select(-genesymbol) %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename(id=rowname)


mapping4db <- mapping_pdb2uniprot_df %>% 
  dplyr::select(pdb, uniprot) %>% 
  dplyr::arrange(pdb, uniprot) %>% 
  distinct() %>% 
  dplyr::left_join(pdb4db, by="pdb") %>% 
  dplyr::rename(id_pdb = id) %>% 
  dplyr::left_join(uniprot4db, by="uniprot") %>% 
  dplyr::rename(id_uniprot = id) %>% 
  dplyr::select(id_pdb, id_uniprot) %>%
  tibble::rownames_to_column() %>% 
  dplyr::rename(id=rowname)



# Write general data into db
dbWriteTable(db_con, "pdb", pdb4db, append=TRUE)
dbWriteTable(db_con, "organisms", organisms4db, append=TRUE)
dbWriteTable(db_con, "genesymbols", genesymbols4db, append=TRUE)
dbWriteTable(db_con, "uniprot", uniprot4db, append=TRUE)
dbWriteTable(db_con, "mapping_pdb_uniprot", mapping4db, append=TRUE)


# First we need to extract some data previously stored (necessary for join)
df_uniprot <- dbGetQuery(db_con, "SELECT * FROM uniprot")

query_create_pathway_databases <- "
CREATE TABLE pathway_databases (
  id INTEGER PRIMARY KEY,
  name TEXT
)"


query_create_pathways <- "
CREATE TABLE pathways (
  id INTEGER PRIMARY KEY,
  code TEXT,
  organism TEXT,
  description TEXT,
  id_pathwaydb INTEGER,
  FOREIGN KEY (id_pathwaydb) 
    REFERENCES pathway_databases (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


query_create_mapping_uniprot_pathway <- "
CREATE TABLE mapping_uniprot_pathway (
  id INTEGER PRIMARY KEY,
  id_uniprot INTEGER,
  id_pathway INTEGER,
  annotation_type TEXT,
  FOREIGN KEY (id_uniprot) 
    REFERENCES uniprot (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION,
  FOREIGN KEY (id_pathway) 
    REFERENCES pathways (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


dbExecute(db_con, query_create_pathway_databases)
dbExecute(db_con, query_create_pathways)
dbExecute(db_con, query_create_mapping_uniprot_pathway)



# 
# Data from Reactome (Pathways and mapping to uniprot)
# 
pathway_databases4db <- data.frame(
  id = c(1, 2),
  name = c("REACTOME", "KEGG"),
  stringsAsFactors = FALSE
)

# 
# Important: 
# The following code is for reactome pathways. 
# These are used on projects where the organism is 'human' 
# We insert here to have a unique database for all the data
# 
# We should add to the database the pathways from KEGG, using the same schema
# for this reason we keep here the code.
# 

# Extract data from Reactome
reactome_pathways_df <- get_reactome_pathways()
reactome_entrezgenes_ls <- get_reactome_entrezgenes()
reactome_mapping_df <- get_reactome_relationships()


# Extract uniprot-gene assignments using biomart 
# Note: needs internet
mapping_uniprot_entrezgenes_df <- get_uniprot_gene_relationships_from_mart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)

# Check number of unique entries
length(unique(mapping_uniprot_entrezgenes_df$uniprot))
length(unique(mapping_uniprot_entrezgenes_df$entrezgene_id))

# How many genes have NOT a correspondence in reactome
sum(!mapping_uniprot_entrezgenes_df$entrezgene_id %in% reactome_entrezgenes_ls)
# How many genes have a correspondence in reactome
sum(mapping_uniprot_entrezgenes_df$entrezgene_id %in% reactome_entrezgenes_ls)

# How many gene entries in reactome
length(reactome_entrezgenes_ls)
# How many have NOT a correspondence extracted by using biomaRt
sum(!reactome_entrezgenes_ls %in% mapping_uniprot_entrezgenes_df$entrezgene_id)
# How many have a correspondence extracted by using biomaRt
sum(reactome_entrezgenes_ls %in% mapping_uniprot_entrezgenes_df$entrezgene_id)


# We should clean and remove the cases in which more uniprot codes correspond 
# to a given gene id
tmp_mapping_df <- mapping_uniprot_entrezgenes_df %>%
  dplyr::group_by(entrezgene_id) %>%
  dplyr::add_count()

clean_mapping_uniprot_entrezgenes_df <- tmp_mapping_df %>%
  dplyr::filter(n==1) %>%
  dplyr::select(-n)

wrong_mapping_uniprot_entrezgenes_df <- tmp_mapping_df %>%
  dplyr::filter(n>1)

# These are UNIPROT codes with multiple mapping. Should they be checked?
cat("These are UNIPROT codes with multiple mapping: should they be checked?\n")
unique(wrong_mapping_uniprot_entrezgenes_df$uniprot)



pathways4db <- reactome_pathways_df %>%
  dplyr::arrange(path_code) %>%
  dplyr::mutate(id_pathwaydb = "1") %>% 
  dplyr::rename(code=path_code,
                organism=path_organism,
                description=path_description) %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename(id=rowname)
# NOTE: this has been implemented with foreign key to pathways databases
#       (Up to now: REACTOME=1, KEGG=2)


mapping_uniprot2reactome4db <- reactome_mapping_df %>% 
  dplyr::left_join(clean_mapping_uniprot_entrezgenes_df, by="entrezgene_id") %>% 
  dplyr::left_join(df_uniprot, by="uniprot") %>%
  # dplyr::left_join(uniprot4db, by="uniprot") %>%
  dplyr::rename(id_uniprot = id, code = pathway) %>% 
  dplyr::left_join(pathways4db, by="code") %>%
  dplyr::rename(id_pathway = id) %>% 
  dplyr::filter(!is.na(id_pathway), !is.na(id_uniprot)) %>% 
  dplyr::select(id_uniprot, id_pathway) %>% 
  dplyr::arrange(id_uniprot, id_pathway) %>% 
  dplyr::distinct() %>%
  tibble::rownames_to_column() %>% 
  dplyr::rename(id=rowname)


# Write general data into db
dbWriteTable(db_con, "pathway_databases", pathway_databases4db, append=TRUE)
dbWriteTable(db_con, "pathways", pathways4db, append=TRUE)
dbWriteTable(db_con, "mapping_uniprot_pathway", mapping_uniprot2reactome4db, append=TRUE)


query_check_data <- "
SELECT 
  pdb.pdb,
  pdb.exp_method,
  pdb.exp_res,
  sym.genesymbol,
  org.organism,
  upr.gene,
  upr.uniprot
  
FROM mapping_pdb_uniprot map
  LEFT JOIN pdb ON pdb.id = map.id_pdb
  LEFT JOIN uniprot upr ON upr.id = map.id_uniprot
  LEFT JOIN organisms org ON org.id = upr.id_organism
  LEFT JOIN genesymbols sym ON sym.id = upr.id_genesymbol
ORDER BY pdb, uniprot;
"

# 
# Create tables for pockets, molecules, molclasses and scores
# 
query_create_pockets <- "
CREATE TABLE pockets (
  id INTEGER PRIMARY KEY,
  name TEXT,
  intid TEXT,
  chain TEXT,
  id_pdb INTEGER,
  FOREIGN KEY (id_pdb) 
    REFERENCES pdb (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


query_create_molecules <- "
CREATE TABLE molecules (
  id INTEGER PRIMARY KEY,
  code TEXT NOT NULL,
  name TEXT NOT NULL,
  info TEXT
)"


query_create_molclasses <- "
CREATE TABLE molclasses (
  id INTEGER PRIMARY KEY,
  molclass TEXT
)"


query_create_molclass_assignment <- "
CREATE TABLE molclass_assignment (
  id INTEGER PRIMARY KEY,
  id_mol INTEGER,
  id_class INTEGER,
  FOREIGN KEY (id_mol)
    REFERENCES molecules (id)
      ON DELETE CASCADE
      ON UPDATE NO ACTION,
  FOREIGN KEY (id_class)
    REFERENCES molclasses (id)
      ON DELETE CASCADE
      ON UPDATE NO ACTION
)"


query_create_scores <- "
CREATE TABLE scores (
  id INTEGER PRIMARY KEY,
  id_poc INTEGER,
  id_mol INTEGER,
  score REAL,
  FOREIGN KEY (id_poc) 
    REFERENCES pockets (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION,
  FOREIGN KEY (id_mol) 
    REFERENCES molecules (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


query_create_stats_mol <- "
CREATE TABLE stats_mol (
  id_mol INTEGER,
  mol_mean REAL,
  mol_sd REAL,
  mol_n INTEGER,
  mol_refvalue REAL,
  FOREIGN KEY (id_mol)
    REFERENCES molecules (id)
      ON DELETE CASCADE
      ON UPDATE NO ACTION
)"


query_create_stats_poc <- "
CREATE TABLE stats_poc (
  id_poc INTEGER,
  poc_mean REAL,
  poc_sd REAL,
  poc_n INTEGER,
  poc_refvalue REAL,
  poc_refcode TEXT,
  FOREIGN KEY (id_poc)
    REFERENCES pockets (id)
      ON DELETE CASCADE
      ON UPDATE NO ACTION
)"


query_create_zscores_mol <- "
CREATE TABLE zscores_mol (
  id_poc INTEGER,
  id_mol INTEGER,
  zscore_mol REAL,
  FOREIGN KEY (id_poc) 
    REFERENCES pockets (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION,
  FOREIGN KEY (id_mol) 
    REFERENCES molecules (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


query_create_zscores_poc <- "
CREATE TABLE zscores_poc (
  id_poc INTEGER,
  id_mol INTEGER,
  zscore_poc REAL,
  poc_refcode TEXT,
  FOREIGN KEY (id_poc) 
    REFERENCES pockets (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION,
  FOREIGN KEY (id_mol) 
    REFERENCES molecules (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


query_create_zzscores_molpoc <- "
CREATE TABLE zzscores_molpoc (
  id_mol INTEGER,
  id_poc INTEGER,
  zzscore_molpoc REAL,
  poc_refcode TEXT,
  FOREIGN KEY (id_mol) 
    REFERENCES molecules (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION,
  FOREIGN KEY (id_poc) 
    REFERENCES pockets (id) 
      ON DELETE CASCADE 
      ON UPDATE NO ACTION
)"


dbExecute(db_con, query_create_pockets)
dbExecute(db_con, query_create_molecules)
dbExecute(db_con, query_create_molclasses)
dbExecute(db_con, query_create_molclass_assignment)
dbExecute(db_con, query_create_scores)
dbExecute(db_con, query_create_stats_mol)
dbExecute(db_con, query_create_stats_poc)
dbExecute(db_con, query_create_zscores_mol)
dbExecute(db_con, query_create_zscores_poc)
dbExecute(db_con, query_create_zzscores_molpoc)


# 
# Add information about molecules ----------------------------------------------
# 


index_molecule_df <- dbGetQuery(db_con, "SELECT id AS id_mol, code AS molname FROM molecules;")
if (nrow(index_molecule_df) > 0) {
  last_mol_id <- max(as.numeric(index_molecule_df$id_mol))
} else {
  last_mol_id <- 0
}

ligands_filename <- file.path(data_dir, "Ligands_data", "LIGANDS.xlsx")
ligands_df <- readxl::read_xlsx(path=ligands_filename) %>%
  dplyr::mutate(molname=as.character(molname)) %>% 
  dplyr::filter(!molname %in% index_molecule_df$molname)

molecules4db <- ligands_df %>%
  dplyr::mutate(id = seq(last_mol_id+1, last_mol_id+nrow(ligands_df))) %>%
  dplyr::select(id, dplyr::everything()) %>%
  dplyr::mutate(code=as.character(molname),
                name=as.character(molname)) %>%
  dplyr::select(-molname, -molclass)

cat("\tDATA READY FOR DATABASE:molecules", nrow(molecules4db), "\n")

dbWriteTable(db_con, "molecules", molecules4db, append=TRUE)

molclasses4db <- ligands_df %>%
  dplyr::select(molclass) %>%
  dplyr::distinct(molclass) %>%
  dplyr::arrange(molclass) %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(id=rowname)

dbWriteTable(db_con, "molclasses", molclasses4db, append=TRUE)


#
# Assign molclass  -------------------------------------------------------------
# 


# From the database extract the index (to be used as FK)
query_index_molecule <- "SELECT id AS id_mol, name AS molname FROM molecules;"
query_index_class <- "SELECT id AS id_class, molclass FROM molclasses;"

index_molecule_df <- dbGetQuery(db_con, query_index_molecule)
index_class_df <- dbGetQuery(db_con, query_index_class)

molclass_assignment4db <- ligands_df %>%
  dplyr::left_join(index_molecule_df, by="molname") %>%
  dplyr::left_join(index_class_df, by="molclass") %>%
  dplyr::filter(!is.na(id_mol) & !is.na(id_class)) %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(id=rowname) %>%
  dplyr::select(id, id_mol, id_class)

dbWriteTable(db_con, "molclass_assignment", molclass_assignment4db, append=TRUE)



#
# Read data from GS.txt files for drugcentral
#
dir_screening_data <- file.path(data_dir, "Biogps_Screening_data")

dir_drucentral <- file.path(dir_screening_data, "drugcentral")
# dir_drucentral <- file.path(dir_screening_data, "drugcentral_tmp")
output_gs_listfiles <- list.files(dir_drucentral, pattern = "*.GS.txt", recursive = TRUE)
output_gs_drucentral <- file.path(dir_drucentral, output_gs_listfiles)

dir_natural <- file.path(dir_screening_data, "natural")
# dir_natural <- file.path(dir_screening_data, "natural_tmp")
output_gs_listfiles <- list.files(dir_natural, pattern = "*.GS.txt", recursive = TRUE)
output_gs_natural <- file.path(dir_natural, output_gs_listfiles)

output_gs_filenames <- c(output_gs_drucentral, output_gs_natural)

# 
# Loop over OUTPUT_GS.txt files ------------------------------------------------
#
for (output_gs_filename in output_gs_filenames) {
  cat("\nWORKING ON FILENAME: ", output_gs_filename, "\n")
  # Extract data (expected columns molname, pocketname and value_num)
  single_scores_df <- read_GS_file(gsfile = output_gs_filename)
  
  # Create a dataframe with info for the pockets
  single_pockets_df <- extract_pockets(scores_df = single_scores_df)
  
  # Add pdb to the database (in case of data are not present)
  eventually_add_pdb(db_con=db_con, single_pockets_df=single_pockets_df)
  
  # Add pockets to the database (in case of data are not present)
  eventually_add_pocket(db_con=db_con, single_pockets_df=single_pockets_df)
  
  # Add molecules to the database (in case of data are not present)
  eventually_add_molecule(db_con=db_con, scores_df=single_scores_df)
  
  # Add score values to the database (do not check for data if already present)
  add_score_values(db_con=db_con, scores_df=single_scores_df)
  
}




# Check added data into table scores -------------------------------------------


# At the end we should check that all the molecules have (about) the same
# number of calculated pockets: the query below has this aim.

query_extract_scores <- "
SELECT mc.molclass,
  sc.id_poc,
  sc.id_mol,
  mol.code,
  mol.name,
  sc.score
FROM scores sc
  LEFT JOIN molclass_assignment ma ON ma.id_mol = sc.id_mol
  LEFT JOIN molecules mol ON mol.id = sc.id_mol
  LEFT JOIN molclasses mc ON mc.id = ma.id_class;"

# Extract scores
dbscores_df <- dbGetQuery(db_con, query_extract_scores)

# Report the number of pockets of each molecule
cat("Below all the molecules should have the same number of pockets:\n")
tabulate(dbscores_df$id_mol) %>% print()


# Print all the molecules and molclasses available in the database
dbscores_df %>% 
  dplyr::select(id_mol, name, molclass) %>% 
  dplyr::distinct() %>% 
  print()


# Parameters for molecules (calculated for all the molecules in the database)
stats_mol4db <- dbscores_df %>%
  dplyr::filter(score > 0) %>%
  dplyr::group_by(id_mol) %>%
  dplyr::summarise(mol_mean=mean(score), mol_sd=sd(score), mol_n=n()) %>%
  dplyr::mutate(mol_refvalue = mol_mean + mol_sd)

dbWriteTable(db_con, "stats_mol", stats_mol4db, append=TRUE)


stats_poc4db <- dbscores_df %>%
  dplyr::filter(score > 0) %>%
  dplyr::filter(molclass == "drugcentral") %>%
  dplyr::group_by(id_poc) %>%
  dplyr::summarise(poc_mean=mean(score), poc_sd=sd(score), poc_n=n()) %>%
  dplyr::mutate(poc_refvalue = poc_mean + poc_sd,
                poc_refcode = "drugcentral")

dbWriteTable(db_con, "stats_poc", stats_poc4db, append=TRUE)


query_calculate_zscore_mol <- "SELECT mol.id AS id_mol
, mol.code AS molcode
, mol.name AS molname
, mcl.molclass
, poc.id AS id_poc
, poc.name AS pocketname
, sco.score
, stm.mol_n
, stm.mol_mean
, stm.mol_sd
, (sco.score-stm.mol_mean)/stm.mol_sd AS zscore_mol
FROM scores sco
LEFT JOIN molecules mol ON mol.id = sco.id_mol
LEFT JOIN molclass_assignment mca ON mca.id_mol = mol.id
LEFT JOIN molclasses mcl ON mcl.id = mca.id_class
LEFT JOIN pockets poc ON poc.id = sco.id_poc
LEFT JOIN stats_mol stm ON stm.id_mol = mol.id
WHERE score > 0"

zscores_mol_df <- dbGetQuery(db_con, query_calculate_zscore_mol)
zscores_mol4db <- zscores_mol_df %>%
  dplyr::select(id_mol, id_poc, zscore_mol)

dbWriteTable(db_con, "zscores_mol", zscores_mol4db, append=TRUE)


# We calculate zscore for natural molecules compared to drugcentral
# We apply two filters: 
# - score should be > 0 
# - number of molecules for the given pocket above a given threshold
poc_n_threshold <- 20
zscores_poc_df <- calculate_zscore_pockets(db_con, wanted_refcode='drugcentral',
                                           wanted_classcode='natural', 
                                           poc_n_threshold=poc_n_threshold)


zscores_poc4db <- zscores_poc_df %>%
  dplyr::select(id_poc, id_mol, poc_refcode, zscore_poc)

# Check the ten most significant zscore_poc values
zscores_poc4db %>% 
  arrange(-zscore_poc) %>% head(10) %>% 
  print()

dbWriteTable(db_con, "zscores_poc", zscores_poc4db, append=TRUE)


query_zzscore <- "
SELECT DISTINCT sco.id_poc
, sco.id_mol
, mol.name AS molname
, mcl.molclass
, sco.score
, zsm.zscore_mol
, zsp.zscore_poc
, zsp.poc_refcode
, zsm.zscore_mol+zsp.zscore_poc AS zzscore_molpoc
FROM scores sco
LEFT JOIN molecules mol ON mol.id = sco.id_mol
LEFT JOIN molclass_assignment mca ON mca.id_mol = sco.id_mol
LEFT JOIN molclasses mcl ON mcl.id = mca.id_class
LEFT JOIN zscores_mol zsm ON zsm.id_poc = sco.id_poc AND zsm.id_mol = sco.id_mol
LEFT JOIN zscores_poc zsp ON zsp.id_poc = sco.id_poc AND zsp.id_mol = sco.id_mol
WHERE (mcl.molclass == 'drugcentral' OR mcl.molclass == 'natural')
AND zscore_mol IS NOT NULL
AND zscore_poc IS NOT NULL
"

zzscore_df <- dbGetQuery(db_con, query_zzscore)

zzscores_molpoc4db <- zzscore_df %>%
  dplyr::filter(score > 0) %>%
  dplyr::select(id_mol, id_poc, poc_refcode, zzscore_molpoc)

dbWriteTable(db_con, "zzscores_molpoc", zzscores_molpoc4db, append=TRUE)




# Da qua è chiaro che per drugcentral non abbiamo tutte le pocket
# dbscores_df %>% 
#   dplyr::group_by(id_mol) %>% 
#   dplyr::add_count() %>%
#   dplyr::select(id_mol, code,n) %>% 
#   dplyr::distinct() %>% print(n=200)




#
# disconnect the DB ------------------------------------------------------------
#

dbDisconnect(db_con) # Disconnect at the end
