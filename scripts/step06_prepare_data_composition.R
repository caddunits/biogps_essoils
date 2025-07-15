#
# Script to prepare data on the composition of the phytocomplex
# 


library(dplyr)
library(magrittr)

repo_dir <- getwd()
data_dir <- file.path(repo_dir, "data")
mapp_dir <- file.path(data_dir, "Mapping_data")

# Input
composition_filename <- file.path(data_dir, "Composition_data", "OLI_COMPDATA.xlsx")

# Output
composition <- file.path(mapp_dir, "composition.RDS")
composition_new <- file.path(mapp_dir, "composition_new.RDS")


# Composition Data from XLS file
df_composition <- readxl::read_excel(path=composition_filename) %>% 
  dplyr::mutate(molname = as.character(molname),
                phytoclass = as.character(phytoclass))

saveRDS(object=df_composition, file=composition)


df_composition_new <- rbind(
  df_composition %>%
    dplyr::filter(cannella > 0) %>%
    dplyr::select(-origano, -timo) %>%
    dplyr::mutate(oil="cannella", oil_eng="Cinnamon") %>%
    dplyr::mutate(perc_comp=as.numeric(cannella)) %>%
    dplyr::select(-cannella),
  df_composition %>%
    dplyr::filter(origano > 0) %>%
    dplyr::select(-cannella, -timo) %>%
    dplyr::mutate(oil="origano", oil_eng="Oregano") %>%
    dplyr::mutate(perc_comp=as.numeric(origano)) %>%
    dplyr::select(-origano),
  df_composition %>%
    dplyr::filter(timo > 0) %>%
    dplyr::select(-cannella, -origano) %>%
    dplyr::mutate(oil="timo", oil_eng="Thyme") %>%
    dplyr::mutate(perc_comp=as.numeric(timo)) %>%
    dplyr::select(-timo)
)

saveRDS(object=df_composition_new, file=composition_new)
