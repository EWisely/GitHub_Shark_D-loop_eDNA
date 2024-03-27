### Create Metabarlist object from Dada2 abundance table, obi3 results output in metabar format, and sample sheet and sequencing metadata files ####
#Author: Eldridge Wisely
#Date: 10-16-2023

#Metabar and Phyloseq section
library("metabaR")
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("tibble")       # Needed for converting column to row names
# Load requested package for plotting
library("reshape2")
library("tidyverse")
library("dplyr")        # filter and reformat data frames


#reads (abundance table)
metabar_reads_table<-read.csv("data/21_Climb2/04_Usearch_annotated/21_Climb_dada2_metabar_reads_table.csv", header = TRUE, strip.white = TRUE)
metabar_reads_table<-column_to_rownames(metabar_reads_table, var="id")
rownames_to_remove_from_reads<- c("ClusterID","ClusterID.1","asv_seqs")

metabar_reads_table<- metabar_reads_table[!(row.names(metabar_reads_table) %in% rownames_to_remove_from_reads),]
reads_colnames<-colnames(metabar_reads_table)
reads_rownames<-rownames(metabar_reads_table)

#fix the accidental replacement of - with . 
reads_colnames <- sapply(reads_colnames, gsub, pattern = "ASV", replacement = "Climb2_ASV", fixed = TRUE)
reads_rownames <- sapply(reads_rownames, gsub, pattern = ".", replacement = "-", fixed =TRUE)

#put it in a numerical matrix then put the rownames and colnames back on
Metabar_formatted_reads<-base::as.matrix(metabar_reads_table, 
                                         ncol = ncol(metabar_reads_table))

Metabar_formatted_reads<-as.numeric(Metabar_formatted_reads)
Metabar_formatted_reads<-base::matrix(Metabar_formatted_reads, 
                                      ncol = ncol(metabar_reads_table))
colnames(Metabar_formatted_reads)<-reads_colnames
rownames(Metabar_formatted_reads)<-reads_rownames



#ASVs/motus
Dada2_ecotagged_usearched_results<-read.csv("data/21_Climb2/04_Usearch_annotated/annotated_23bthaps_2021_results.tab", sep="\t")
Dada2_ecotagged_usearched_results<- column_to_rownames(Dada2_ecotagged_usearched_results, var = "ID")
Dada2_ecotagged_usearched_results<-as.data.frame(Dada2_ecotagged_usearched_results)

write.csv(Dada2_ecotagged_usearched_results, "data/21_Climb2/04_Usearch_annotated/annotated_23bthaps_2021_results.csv")

#aggregate annotations into three new columns, then re-import.
Dada2_ecotagged_usearched_results1<-read.csv("data/21_Climb2/04_Usearch_annotated/annotated_23bthaps_2021_results_aggregated_annots.csv")
Dada2_ecotagged_usearched_results1<- column_to_rownames(Dada2_ecotagged_usearched_results1, var = "X")
Dada2_ecotagged_usearched_results1<-as.data.frame(Dada2_ecotagged_usearched_results1)


Metabar_Dada2_ecotagged_usearched_results <- Dada2_ecotagged_usearched_results1 %>% 
  dplyr::rename(sequence = NUC_SEQ)
Metabar_Dada2_ecotagged_usearched_results = subset(Metabar_Dada2_ecotagged_usearched_results, select = -c(COUNT, DEFINITION) )
Metabar_Dada2_ecotagged_usearched_results$seq_len<-nchar(Metabar_Dada2_ecotagged_usearched_results$sequence)
#rename motus: 
Metabar_formatted_motus<-Metabar_Dada2_ecotagged_usearched_results
#put motus in numerical order like they are in the reads object.
Metabar_formatted_motus$ASV_order = sapply(strsplit(rownames(Metabar_Dada2_ecotagged_usearched_results), "ASV_"), "[[", 2)
#order motus same as reads
Metabar_formatted_motus<- Metabar_formatted_motus[colnames(Metabar_formatted_reads),]


#samples

#edited the sample names from the Shark Project Sample sheet.xlxs to make it match the PCRs file (by removing spaces and shortening names as was necessary for the sequencing core and saved as "data/metadata/22-23_Climb2_sample_metadata.xlsx" and saved as a csv as well for import.
Metabar_formatted_samples<-read.csv("data/21_Climb2/04_Usearch_annotated/Climb2_21_samples.csv", na.strings = c("NA",""))
Metabar_formatted_samples<-column_to_rownames(Metabar_formatted_samples, var="Sample")



#pcrs

pcrs <- data.frame(
  sample_id = sapply(strsplit(rownames(Metabar_formatted_reads), "_"), "[[", 1),
  #rep = sapply(strsplit(rownames(Metabar_formatted_reads), "Slew"), "[[", 2),
  type = "sample",
  control_type = "NA",
  row.names = rownames(Metabar_formatted_reads)
)
pcrs_to_edit<-rownames_to_column(pcrs, var = "id")

write_csv(pcrs_to_edit, "data/21_Climb2/04_Usearch_annotated/21_Climb2_pcrs_to_edit.csv")
#edit in excel to make sure that type and control type are set correctly! save as tab-delimited text file with the following name: Metabar_formatted_PCRs_edited.txt
pcrs<- read_tsv("data/21_Climb2/04_Usearch_annotated/21_Climb2_pcrs_edited.txt", col_names = TRUE)
Metabar_formatted_pcrs<-column_to_rownames(pcrs, var="id")
Metabar_formatted_pcrs<- as.data.frame(Metabar_formatted_pcrs)


#Write intermediate files

write.csv(Metabar_formatted_reads, "data/21_Climb2/05_MetabaR_filtered/Metabar_formatted_reads.csv")
write.csv(Metabar_formatted_motus, "data/21_Climb2/05_MetabaR_filtered/Metabar_formatted_motus.csv")
write.csv(Metabar_formatted_pcrs, "data/21_Climb2/05_MetabaR_filtered/Metabar_formatted_pcrs_Climb_2021.csv")
write.csv(Metabar_formatted_samples, "data/21_Climb2/05_MetabaR_filtered/Metabar_formatted_samples.csv")


Climb_21<-metabarlist_generator(reads=Metabar_formatted_reads, Metabar_formatted_motus, pcrs=Metabar_formatted_pcrs, samples=Metabar_formatted_samples)


summary_metabarlist(Climb_21)

Climb_21$motus$count = colSums(Climb_21$reads)

#Arguments
#reads	
#MOTU abundance table. Rows and rownames of the table should correspond to PCRs and their names respectively. Columns and colnames should correspond to MOTUs and their names. Rownames in this table should correspond to PCR names respectively.

#motus	
#MOTU characteristics table (e.g. taxonomy, sequence, etc.). Rows and rownames of the table should correspond to MOTUs and their names respectively, and the columns to their characteristics. Mandatory fields: 'sequence', i.e. the sequence representative of the MOTU.

#pcrs	
#PCR characteristics table (e.g. tags, primers, plate wells, etc.). Rows and rownames of the table should correspond to PCRs and their names respectively, and the columns to their characteristics. Mandatory fields: (i) 'sample_id', i.e. the name of each biological sample. (ii) 'type', i.e. the type of PCR; can be 'sample' or 'control'. (iii) 'control_type', i.e. the type of control if applicable. Should be either: 'NA' for samples, 'extraction' for extraction negative controls, 'pcr' for PCR negative controls, 'sequencing' for sequencing negative controls (e.g. unused tag combinations), or 'positive' for positive controls.

#samples	
#Sample characteristics table. Rows and rownames of the table should correspond to biological samples and their names respectively, and the columns to their environnemental characteristics.


# Compute the number of reads per pcr
Climb_21$pcrs$nb_reads <- rowSums(Climb_21$reads)

# Compute the number of motus per pcr
Climb_21$pcrs$nb_motus <- rowSums(Climb_21$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_21 <- subset_metabarlist(Climb_21, table = "pcrs",
                                  indices = Climb_21$pcrs$nb_reads>0)

summary_metabarlist(Climb_21)
save(Climb_21, file = "data/21_Climb2/05_MetabaR_filtered/Climb_21_Metabarlist.RData")

