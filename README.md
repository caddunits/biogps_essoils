# biogps_essoils
BioGPS application for phytocomplexes from essential oils

### Prepare Data

#### Data mapping Uniprot <-> ATCC
1. Run the script step01_prepare_data_atcc_mapping.R

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
2. Run the script step02a_read_kegg_data.R (to be run with several codes)
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
use the additional script step02b_check_kegg_data.R 
(it produces an xlsx file, with the information - numer of uniprot, 
number of pathways, ... available as supplementary information).


3. Run the script step03_prepare_data_kegg_mapping.R
This script creates a file named "kegg_bacteria_uniprot.RDS"
which is stored in the directory data/Mapping_data
and contains the information from all the .RDS files obtained for the different 
bacteria (previous point).

# ATTENZIONE: NELLA PARTE FINALE PROPOSTA MODIFICA PER SEMPLIFICARE SEQUENZA, 
# DA VALUTARE ED AGGIORNARE CODICE


#### Data mapping Uniprot <-> Genes
4. Run the script step04_prepare_data_genes_mapping.R
This script creates a file named "gene_uniprot.RDS"
The mapping is limited to the Uniprot entries which have a correspondence in
the PDB archive. In fact, data are taken from the web (file pdbtosp.txt from Uniprot web site)
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pdbtosp.txt
A data mining procedure elaborates the data by including only some methods:
XRAY, NMR, EM, NEUTRON, IR, FIBER, OTHER.
The data contains the mapping for all the organisms, not only those studied
in the present work.


#### Data from bacteria
5. Run the script step05_prepare_data_bacteria.R
This script creates a file named "bacteria_data.RDS" which contains
info on bacteria, including codes on different databases (ATCC, KEGG)
Input data is a manually curated file (BACTERIADATA.xlsx) starting from the 
five-letter organism codes extracted from the dataframe previously obtained; 
information for data curation is obtained from Uniprot.


#### Data from composition
6. Run the script step06_prepare_data_composition.R
This script creates a file named "composition.RDS"
and another file names "composition_new.RDS"
It extracts data about composition from an excel file. There are columns related
to the presence (as percentage) in the three oils, with the value 0.0333 we mean 
'traces'.
The file composition_new is in a wider format, with one column to define the
numerical value and another column to define the oil; instead, the file 
composition has one column for each oil.
In both files there is also the phytochemical classification.


#### Data from BioGPS calculations
Run BioGPS commandline for the list of pockets and the database of molecules. 
It can be convenient to split the job in several pieces, because at the end 
of the job it takes several time to save data and exit.
    
7. Run the script step07_prepare_db_biogps.R
The script first creates an empty database and then populates with some data
for the data linked to the pockets (pdb, uniprot, genes, pathways).

8. Run the script step08_prepare_data_biogps.R
The script reads data from the database (zzscores_molpoc, zscores_mol and zscores_poc)
and saves an object in the format .RDS to be used in the analysis.


### Check Data

#### Data for the table in the paper with all the bacteria numeric information 
9. Run the script step09_check_tablesummary.R
Note that for the part related to biogps it is necessary the pocketome database
(kindly provided by Molecular Discovery Ltd).
However, this is not necessary for the analysis



### Analysis

#### Target fishing and network analysis
10. Run the script step10_targetfishing_analysis.R
(needs internet connection to download data from stringdb)

Data are read from the following .RDS files:
biogps_data.RDS
bacteria_data.RDS
composition_new.RDS

Output are tabular data, including network centrality, saved in the file
Results_targetfishing.xlsx (in the directory output/targetfishing)


#### Results as barplots and heatmaps
11. Run the script step11_targetfishing_barplots_and_heatmaps.R

Data are read from the following files:
composition_new.RDS
Results_targetfishing.xlsx (from step10)

Additional input data is in the file manual_curation.xlsx 
(available in the directory data)

Several sheets are present, including:
multiple_genes
proteins (correct characters -uppercase/lowercase- given that all 
the genesymbols are stored in the databsase as upercase)
target_curation
pdb_curation
pubmed
pocket_curation
The sheet Non_ecoli reports data tha were removed 
(originally present for some errors)

In the sheet multiple_genes, for each bacteria we report genes as groups 
whenever they correspond to the same complex or should be considered as
duplicates for any other reason. 
Some genes from the analysis that have this peculiarity are the following:
Staphylococcus aureus	ACCA	ACCA/ACCD
Staphylococcus aureus	ACCD	ACCA/ACCD
Pseudomonas aeruginosa	AMIC	AMIC/AMIR
Pseudomonas aeruginosa	AMIR	AMIC/AMIR
Pseudomonas aeruginosa	PQSB	PQSB/PQSC
Pseudomonas aeruginosa	PQSC	PQSB/PQSC
Pseudomonas aeruginosa	GSPI	GSPI/GSPJ/GSPK
Pseudomonas aeruginosa	GSPJ	GSPI/GSPJ/GSPK
Pseudomonas aeruginosa	GSPK	GSPI/GSPJ/GSPK
Pseudomonas aeruginosa	NORB	NORB/NORC
Pseudomonas aeruginosa	NORC	NORB/NORC
Escherichia coli	MCBA	MCBA/MCBB/MCBC/MCBD
Escherichia coli	MCBB	MCBA/MCBB/MCBC/MCBD
Escherichia coli	MCBC	MCBA/MCBB/MCBC/MCBD
Escherichia coli	MCBD	MCBA/MCBB/MCBC/MCBD
Escherichia coli	CARA	CARA/CARB
Escherichia coli	CARB	CARA/CARB
Escherichia coli	NARG	NARG/NARH/NARI
Escherichia coli	NARH	NARG/NARH/NARI
Escherichia coli	NARI	NARG/NARH/NARI
...

In the script these are then considered as the same genegroup and
the higher centrality value along the genegroup is assigned
In a future development, this part should be automathized

The script produces the two barplots (for major and minor bacteria) and
the heatmaps for the major bacteria

Data are read from the XLSX file of the previous step
First, data are considered focusing on molecules; the question is:
How many bacterial targets does each molecule hit?
Second, data are considered focusing on phytocomplexes; the question is:
How do the phytocomplexes interact with each target (through which molecules)?

We calculate a contribution for each target, by considering the zzscore and the
percentage (composition). Then, we calculate the sum of contributions over all
the molecules of phytocomplexes, and add the data about centrality
We store the data of zzscore from molecule-target of oil-bacterium pairs and
then extract the top n genes for each bacteria/oil pair (can adjust this value)

So, after having created a unique list of genes per bacteria (i.e., each 
bacteria's top n genes) we use these shortlist of genes (in dataframe format) 
to filter the previous data. 

In order to get a shared list among the bacteria, we calculate the sum of
`sum_contr` over the three oils for each gene; this will be used to rank genes 
in the Y axis of the heatmaps and scatterplot.

Thus, we finally merge the ranked genes back with the filtered dataset
(which is now ranked).

In order to solve the problem of genes that should be considered as a goup, we
create a new dataframe in which we add a new column, that is named 'gene'.
This column is the name of grouped genes (whenever available), otherwise it is
the old genesymbol name. Then we use this (column gene) with bacteria and
oil_eng to group rows and take only one row for group (with highest centrality).

