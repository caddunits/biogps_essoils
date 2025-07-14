# biogps_essoils
BioGPS application for phytocomplexes from essential oils

### Prepare Data

#### Data mapping Uniprot <-> ATCC
1. Run the script step1_prepare_data_atcc_mapping.R

This script is used to create a file named "atcc_bacteria_uniprot.RDS"
which is stored in the directory data/Mapping_data

NOTE: according to ATCC website, the download of CSV files is available only to 
"ATCC Genome Portal supporting members or those who have purchased a 
 corresponding physical product".
Therefore, in the official repository of the current project we cannot put
CSV files, but only a modified version, as .RDS file with only the necessary 
columns for running the other scripts. 

However, We keep in the repository the script to create the file RDS from the
set of CSV files.


#### Data mapping Uniprot <-> KEGG pathways
2. Run the script step2a_read_kegg_data.R (to be run with several codes)
This script was used to extract KEGG data corresponding to different values of
org_code (retrieved in the website https://www.kegg.jp/brite/br08611)

For some values of org_code, the script may fail, because of some problems 
with the access, with the following error message
"Errore in .getUrl(url, .flatFileParser) : Forbidden (HTTP 403)."
The problem was solved by adding Sys.sleep(x), with x=0.3

The script saves data into a RDS file (one for each org_code). We cannot put 
in the repository the RDS data, but only the script to create the file RDS by
accessing the KEGG webserver through the R package KEGGREST keggList(org_code)

For each org_code, at the end of the script we report the success_rate, 
a percentage value which estimates how many objects have been retrieved, 
by running the command keggList(org_code)
This data is reported as supplementary information in the paper. 

In order to collect all the information from different .RDS files,
use the additional script step2b_check_kegg_data.R 
(it produces an xlsx file, with the information - numer of uniprot, 
number of pathways, ... available as supplementary information).


3. Run the script step3_prepare_data_kegg_mapping.R
This script creates a file named "kegg_bacteria_uniprot.RDS"
which is stored in the directory data/Mapping_data
and contains the information from all the .RDS files obtained for the different 
bacteria (previous point).

# ATTENZIONE: NELLA PARTE FINALE PROPOSTA MODIFICA PER SEMPLIFICARE SEQUENZA, 
# DA VALUTARE ED AGGIORNARE CODICE


#### Data mapping Uniprot <-> Genes
4. Run the script step4_prepare_data_genes_mapping.R
This script creates a file named "gene_uniprot.RDS"
The mapping is limited to the Uniprot entries which have a correspondence in
the PDB archive. In fact, data are taken from the web (file pdbtosp.txt from Uniprot web site)
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pdbtosp.txt
A data mining procedure elaborates the data by including only some methods:
XRAY, NMR, EM, NEUTRON, IR, FIBER, OTHER.
The data contains the mapping for all the organisms, not only those studied
in the present work.


#### Data from bacteria
5. Run the script step5_prepare_data_bacteria.R
This script creates a file named "bacteria_data.RDS" which contains
info on bacteria, including codes on different databases (ATCC, KEGG)
Input data is a manually curated file (BACTERIADATA.xlsx) starting from the 
five-letter organism codes extracted from the dataframe previously obtained; 
information for data curation is obtained from Uniprot.

