#export final dataset (just the annotated ones, disregarding obi3 assignment status for phyloseq (ASVs/haplotype can be merged in phyloseq) 

library(dplyr)

#reads
df.FINAL.reads_Climb_21_hap_matched<- as.data.frame(Climb_21_hap_matched_clean$reads)
write.csv(df.FINAL.reads_Climb_21_hap_matched, "data/21_Climb2/05_MetabaR_filtered/post_MetabaR_reads_Climb_21_hap_matched_ASVs.csv")
#make the sample names match those in the samples file by removing the primer name in front
df.FINAL.reads_Climb_21_hap_matched<-as.data.frame(read.csv("data/21_Climb2/05_MetabaR_filtered/post_MetabaR_reads_Climb_21_hap_matched_ASVs1.csv"))

#motus
df.FINAL.motus_Climb_21_hap_matched<- as.data.frame(Climb_21_hap_matched_clean$motus)
write.csv(df.FINAL.motus_Climb_21_hap_matched, "data/21_Climb2/05_MetabaR_filtered/post_MetabaR_ASVs_Climb_21_hap_matched_ASVs.csv")


#samples
df.FINAL.samples_Climb_21_hap_matched<- as.data.frame(Climb_21_hap_matched_clean$samples)
write.csv(df.FINAL.samples_Climb_21_hap_matched, "data/21_Climb2/05_MetabaR_filtered/post_MetabaR_samples_Climb_21_hap_matched_ASVs.csv")

#pcrs (the metadata in here should be unnecessary at this point)
df.FINAL.pcrs_Climb_21_hap_matched<- as.data.frame(Climb_21_hap_matched_clean$pcrs)
write.csv(df.FINAL.pcrs_Climb_21_hap_matched, "data/21_Climb2/05_MetabaR_filtered/post_MetabaR_pcrs_Climb_21_amplified_hap_matched_ASVs.csv")


#Save R Global Environment
save.image(file = "data/21_Climb2/05_MetabaR_filtered/post_MetabaR_Climb_21_hap_matched_ASVs.RData")

#Full Dataset Phyloseq starting from existing R objects made by the above code:



#############

#Use View(GAL21_Diet_Combined_agg_sum1[["pcrs"]]) and View(GAL21_Diet_Combined[["pcrs"]]) to figure out the rarefaction curve (sortof) from 1 rep to 2 reps to all three (if the three were sequenced separately)

###############

#Phyloseq!
#otu_table - Works on any numeric matrix. You must also specify if the species are rows or columns
#sample_data - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object
#tax_table - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
#phyloseq - Takes as argument an otu_table and any unordered list of valid phyloseq components: sample_data, tax_table, phylo, or XStringSet. The tip labels of a phylo-object (tree) must match the OTU names of the otu_table, and similarly, the sequence names of an XStringSet object must match the OTU names of the otu_table.
#merge_phyloseq - Can take any number of phyloseq objects and/or phyloseq components, and attempts to combine them into one larger phyloseq object. This is most-useful for adding separately-imported components to an already-created phyloseq object.


load("data/21_Climb2/05_MetabaR_filtered/post_MetabaR_Climb_21_hap_matched_ASVs.RData")
#otu table

#clean.otu.df<-as.data.frame(t(GAL21_Diet_Combined1$reads))
#change this to be the aggregated_mean of pcrs into samples dataset

clean.otu.df<-as.data.frame(t(Climb_21_hap_matched_clean$reads))
clean.otu.df<-rownames_to_column(clean.otu.df, var="id")

#taxa table

library(taxonomizr)
prepareDatabase(getAccessions=FALSE)
#make a taxa table for all obi3 result ASVs
taxaId<-df.FINAL.motus_Climb_21_hap_matched$TAXID
hapID<-df.FINAL.motus_Climb_21_hap_matched$annot_hap
#uncomment the 2 lines below to make a taxa table for just the final aggregated data
#df.FINAL.motus_Climb_21_hap_matched<-rownames_to_column(df.FINAL.motus_GAL21_Diet_Combined_aggregated, var = "id")
#taxaId<-df.FINAL.motus_GAL21_Diet_Combined_aggregated1$TAXID

taxa<-getTaxonomy(taxaId,'nameNode.sqlite') #error here.  Good thing I don't need this for the C. limbatus haplotypes.
print(taxa)
class(taxa)
taxa.df<-as.data.frame.array(taxa)

#rownames(taxa.df)<-df.21GAL_Diet_results_merged_motus$id
rownames(taxa.df)<-df.FINAL.motus_Climb_21_hap_matched$id

#backup plan for haplotypes instead of full biodiversity dataset:
haplo.df<- df.FINAL.motus_Climb_21_hap_matched %>% select(TAXID, SCIENTIFIC_NAME, BEST_IDENTITY, annot_hap, annot_qual, annot_pctid)


#samples table

df.FINAL.samples_Climb_21_hap_matched<-rownames_to_column(df.FINAL.samples_Climb_21_hap_matched, var="sample_id")

samples.df<-full_join(df.FINAL.pcrs_Climb_21_hap_matched,df.FINAL.samples_Climb_21_hap_matched, by="sample_id")


samples.df <- samples.df %>% 
  tibble::column_to_rownames("sample_id") 

#adjust these values below to select just the needed columns
#samples.df<-samples.df[c(1,10,32:61)]

#uncomment the line below for the rep_aggregated dataset

#samples.df<-df.FINAL.samples_GAL21_Diet_Combined_aggregated



#make otu matrix

clean.otu.df<-column_to_rownames(clean.otu.df, var="id")
clean.otu.mat<-as.matrix(clean.otu.df)

#make taxa matrix #changed taxa.df below to haplo.df
clean.taxa.mat<-as.matrix(haplo.df)


#Save and Load R Global Environment
save.image(file = "data/21_Climb2/06_Phyloseq_viz/MetabaR_output_ready_for_Phyloseq")
#load("data/21_Climb2/06_Phyloseq_viz/MetabaR_output_ready_for_Phyloseq")


