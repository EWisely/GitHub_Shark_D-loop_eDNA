#use Phyloseq to visualize Climb eDNA data
#Author Eldridge Wisely
#Date: 26-03-24

#Load R Global Environment
load("data/21_Climb2/06_Phyloseq_viz/MetabaR_output_ready_for_Phyloseq")

#install.packages("ggraph") # for taxatree_plots()
#install.packages("DT") # for tax_fix_interactive()
#install.packages("corncob") # for example datasets and beta binomial models
library(microViz)
library(ggraph)
library(DT)
library(corncob)
library(metagMisc)
library(phyloseq)


#Starting Phlyoseq!
# Climb_21_haplo.ps<- phyloseq(otu_table(clean.otu.mat, taxa_are_rows = TRUE), 
#                          sample_data(samples.df), 
#                          tax_table(clean.taxa.mat))

Climb_21_haplo.ps<-readRDS("data/21_Climb2/06_Phyloseq_viz/Climb_21_haplo.ps.RDS")
Climb_21_haplo.ps #(done with usearch haplotype identified data not considering the obitools3 taxonomic identification as well.)

Climb_21_species.ps<-readRDS("data/21_Climb2/06_Phyloseq_viz/Climb_21_species.ps.RDS")
Climb_21_species.ps

Climb_21_species_merge.ps<-tax_glom(Climb_21_species.ps, "species", NArm = TRUE)

#make Site_Date column

df.samples.Climb_21_species_merge.ps<-sample_data(Climb_21_species_merge.ps)

df.samples.Climb_21_species_merge.ps.1<-df.samples.Climb_21_species_merge.ps%>%
  mutate(Site_Date= paste(Site,Date))
df.samples.Climb_21_species_merge.ps.1$Site_Date

sample_data(Climb_21_species_merge.ps)<-df.samples.Climb_21_species_merge.ps.1

Climb1_species_merge_samples.df<-as.data.frame(cbind(sample_data(Climb_21_species_merge.ps)))
write.csv(Climb1_species_merge_samples.df,"data/21_Climb2/06_Phyloseq_viz/2021_Climb_species_merged_samples.csv")

#samples

Climb21_species_merge_samples.df<-as.data.frame(cbind(sample_data(Climb_21_species_merge.ps)))
write.csv(Climb21_species_merge_samples.df,"data/21_Climb2/06_Phyloseq_viz/2021_Climb_species_merged_samples.csv")

#reads

Climb21_species_merge_reads.df<-as.data.frame(t(cbind(otu_table(Climb_21_species_merge.ps))))
write.csv(Climb21_species_merge_reads.df,"data/21_Climb2/06_Phyloseq_viz/2021_Climb_species_merged_reads.csv")

# Transform data 

#species detection object

Climb_21_species_merge.ps.bin <- phyloseq_standardize_otu_abundance(Climb_21_species_merge.ps,"pa")

Climb_21_species_merge.eDNAindex.ps<-phyloseq_standardize_otu_abundance(Climb_21_species_merge.ps, "wisconsin")

Climb_21_species_merge.ps.prop <- phyloseq_standardize_otu_abundance(Climb_21_species_merge.ps, method = "total")

Climb_21_species_merge.hell.ps<- phyloseq_standardize_otu_abundance(Climb_21_species_merge.ps, method = "hellinger")

Climb_21_species_merge.norm.ps<- phyloseq_standardize_otu_abundance(Climb_21_species_merge.ps, method = "normalize")


#plot proportions of haplotype reads by number of ASVs by Site
plot_bar(Climb_21_species_merge.eDNAindex.ps, x="Site_Date",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Species_detections_eDNAindex_of_ASVs_by_Site_Date.jpg")

#plot proportions of haplotype reads by number of ASVs by Site
plot_bar(Climb_21_species_merge.ps.bin, x="Site_Date",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Species_detections_proportions_of_ASVs_by_Site_Date.jpg")

plot_bar(Climb_21_species_merge.ps.bin, x="Associated_Blood_Samples",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Species_detections_proportions_of_ASVs_by_Blood_Sampling.jpg")




#haplotypes
Climb_21_haplo.ps.bin <- phyloseq_standardize_otu_abundance(Climb_21_haplo.ps,"pa")

Climb_21_haplo.eDNAindex.ps<-phyloseq_standardize_otu_abundance(Climb_21_haplo.ps, "wisconsin")

Climb_21_haplo.ps.prop <- phyloseq_standardize_otu_abundance(Climb_21_haplo.ps, method = "total")

Climb_21_haplo.hell.ps<- phyloseq_standardize_otu_abundance(Climb_21_haplo.ps, method = "hellinger")

Climb_21_haplo.norm.ps<- phyloseq_standardize_otu_abundance(Climb_21_haplo.ps, method = "normalize")



#plot proportions of haplotype reads by number of ASVs by Site
plot_bar(Climb_21_haplo.ps.prop, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_proportions_of_ASVs_by_Site.jpg")

#plot haplotype reads by number of ASVs by Site
plot_bar(Climb_21_haplo.ps, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_reads_of ASVs_by_Site.jpg")

#plot eDNAindex of haplotype reads by number of ASVs by Site
plot_bar(Climb_21_haplo.eDNAindex.ps, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_eDNAindex_of_ASVs_by_Site.jpg")

#plot hellinger transformation of haplotype reads by number of ASVs by Site
plot_bar(Climb_21_haplo.hell.ps, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_hellinger_of_ASVs_by_Site.jpg")


#plot haplotype ASV presence by Site (not aggregated by haplotype)
plot_bar(Climb_21_haplo.ps.bin, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_ASV_presence_by_Site.jpg")


#######Now do the same but sorted by Date #####

#plot proportions of haplotype reads by number of ASVs by Date
plot_bar(Climb_21_haplo.ps.prop, x="Date",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_proportions_of_ASVs_by_Date.jpg")

#plot haplotype reads by number of ASVs by Date
plot_bar(Climb_21_haplo.ps, x="Date",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_reads_of ASVs_by_Date.jpg")

#plot eDNAindex of haplotype reads by number of ASVs by Date
plot_bar(Climb_21_haplo.eDNAindex.ps, x="Date",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_eDNAindex_of_ASVs_by_Date.jpg")

#plot Hellinger transformation of haplotype reads by number of ASVs by Date
plot_bar(Climb_21_haplo.hell.ps, x="Date",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_Hellinger_trans_of_ASVs_by_Date.jpg")



#plot haplotype ASV presence by Date (not aggregated by haplotype)
plot_bar(Climb_21_haplo.ps.bin, x="Date",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_ASV_presence_by_Date.jpg")

###############

#make a site-date column for the haplo phyloseq object (proportions)

df.samples.Climb_21_haplo.ps.prop<-sample_data(Climb_21_haplo.ps.prop)

df.samples.Climb_21_haplo.ps.prop1<-df.samples.Climb_21_haplo.ps.prop%>%
  mutate(Site_Date= paste(Site,Date))
df.samples.Climb_21_haplo.ps.prop1$Site_Date

sample_data(Climb_21_haplo.ps.prop)<-df.samples.Climb_21_haplo.ps.prop1
sample_data(Climb_21_haplo.eDNAindex.ps)<-df.samples.Climb_21_haplo.ps.prop1


#plot proportions of haplotype reads by number of ASVs by Site_Date
plot_bar(Climb_21_haplo.ps.prop, x="Site_Date",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_proportions_of_ASVs_by_Site-Date.jpg")


#plot eDNAindex of haplotype reads by number of ASVs by Site_Date
plot_bar(Climb_21_haplo.eDNAindex.ps, x="Site_Date",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/21_Climb2/06_Phyloseq_viz/2021_Bar_Chart_of_Climb_eDNA_Haplotype_eDNAindex_of_ASVs_by_Site-Date.jpg")


library(RColorBrewer)


#Normalize number of reads in each sample using median sequencing depth.

#I thought this would help, but now everything is really close to zero.

total = median(sample_sums(Climb_21_haplo.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
Climb_21_haplo.ps.med.norm = transform_sample_counts(Climb_21_haplo.ps, standf)

plot_bar(Climb_21_haplo.ps.med.norm, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")




#Faceting!


plot_bar(Climb_21_haplo.ps.prop, fill = "annot_hap", facet_grid = "Site")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")

ggsave("data/21_Climb2/06_Phyloseq_viz/Faceted_sample_Site_proportions_annot_hap.jpg")

plot_bar(Climb_21_haplo.eDNAindex.ps, fill = "annot_hap", facet_grid = "Site")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")

ggsave("data/21_Climb2/06_Phyloseq_viz/Faceted_sample_Site_eDNAindex_annot_hap.jpg")

plot_bar(Climb_21_haplo.ps, fill = "annot_hap", facet_grid = "Site")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")

ggsave("data/21_Climb2/06_Phyloseq_viz/Faceted_sample_Site_Reads_of_ASVs_annot_hap.jpg")

plot_bar(Climb_21_haplo.ps.bin, fill = "annot_hap", facet_grid = "Site")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")

ggsave("data/21_Climb2/06_Phyloseq_viz/Faceted_sample_Site_freq_of_ASVs_annot_hap.jpg")

#This is the most important part to me to make the otu table readable.
#Fist I'll have to put the annot_hap column last
#merge motus by annot_hap  
Climb_21_hap_merge.ps <- tax_glom(Climb_21_haplo.ps, taxrank = "annot_hap", NArm = FALSE)
plot_bar(Climb_21_hap_merge.ps, fill = "annot_hap") +
geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
# 
ggsave("data/21_Climb2/06_Phyloseq_viz/Sample_ASV_read_counts_by_annot_hap.jpg")

# 
# Climb_21_hap_merge.ps
# sample_variables(Climb_21_hap_merge.ps)
# rank_names(Climb_21_hap_merge.ps)
# tax_table(Climb_21_hap_merge.ps)

#calculate Bray-Curtis or Sørensen dissimilarity between blood and associated waters.


