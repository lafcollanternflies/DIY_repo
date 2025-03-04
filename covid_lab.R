# We only need a single package
library(tidyverse)
# Begin by reading in study design that includes all info for both human and ferret samples
targets <- read_tsv("covid_metadata.txt")
# then read in the human covid data and convert to a matrix with gene symbols as rownames
human_covid_data <- read_tsv("GSE147507_RawReadCounts_Human.tsv")
human_covid_data <- as.matrix(column_to_rownames(human_covid_data, "...1")) %>% as.data.frame(.)
# repeat for the ferret covid data
ferret_covid_data <- read_tsv("GSE147507_RawReadCounts_Ferret.tsv")
ferret_covid_data <- as.matrix(column_to_rownames(ferret_covid_data, "...1")) %>% as.data.frame(.)

# Now proceed with your exploration and analysis of the data!
library(ensembldb) 
library(EnsDb.Hsapiens.v86) 
library(biomaRt)

####Question: 
#Comparing the Immune Response (CCL20, IL6, and TNF) of SARS-CoV-2 infected A549 cells (Human lung cell line) vs SARS-CoV-2 infected 4 month old Ferret
genes_of_interest <- c("CCL20","IL6","TNF","IL6R","TLR3","IL17A","ACE2","IRF3","IFNB1","IFNAR1") #The only ones I found in ferret someone highlighted in paper were first 3, added ones that are known from research to have some sort pf resp in covid
fer.samples_of_interest <- c("Series10_FerretNW_SARS-CoV-2_d1_1","Series10_FerretNW_SARS-CoV-2_d1_2","Series11_FerretNW_SARS-CoV-2_d3_1","Series11_FerretNW_SARS-CoV-2_d3_2") #To make simplier, I onlt chose a few, being the ones with the most initial respose (1 day, 3 day) since human data was taken 24 hrs
human.samples_of_interest <- c("Series5_A549_SARS-CoV-2_1","Series5_A549_SARS-CoV-2_2","Series5_A549_SARS-CoV-2_3","Series2_A549_SARS-CoV-2_3")

#I had fer.anno and fer.attributes previously from lab
fer.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mpfuro_gene_ensembl")
fer.attributes <- listAttributes(fer.anno)
#need to map gene ID to the transcript in the rownames of the ferret covid data
Tx.fer <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id','external_gene_name',
                             'description', 'entrezgene_id', 'pfam'),
                mart = fer.anno)

Tx.fer <- as_tibble(Tx.fer)
Tx.want.fer <- Tx.fer %>% filter(external_gene_name%in%genes_of_interest) #Getting the Gene IDs so we can filter down the large experssion data

fer.expression <- ferret_covid_data %>% filter(rownames(.)%in%Tx.want.fer$ensembl_gene_id) %>% dplyr::select(any_of(fer.samples_of_interest))
human.expression <- human_covid_data %>% filter(rownames(.)%in%genes_of_interest) %>% dplyr::select(any_of(human.samples_of_interest))

#rename ferret rows to gene name
gene_name_mapping <- Tx.want.fer %>%
    distinct(ensembl_gene_id, external_gene_name) %>%  # Ensure unique mapping
    filter(ensembl_gene_id %in% rownames(fer.expression)) %>%  # Keep only relevant genes
    deframe() 
rownames(fer.expression) <- gene_name_mapping[rownames(fer.expression)]

#need to normalize I think to compare between species... not sure how
