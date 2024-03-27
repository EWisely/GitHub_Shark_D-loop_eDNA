#use Phyloseq to visualize Climb eDNA data
#Author Eldridge Wisely
#Date: 26-03-24

#Load R Global Environment
load("data/22-23_Climb2_amplified/06_Phyloseq_viz/MetabaR_output_ready_for_Phyloseq")

#install.packages("ggraph") # for taxatree_plots()
#install.packages("DT") # for tax_fix_interactive()
#install.packages("corncob") # for example datasets and beta binomial models
library(microViz)
library(ggraph)
library(DT)
library(corncob)



#Starting Phlyoseq!
Climb_22_23.ps<- phyloseq(otu_table(clean.otu.mat, taxa_are_rows = TRUE), 
                         sample_data(samples.df), 
                         tax_table(clean.taxa.mat))

Climb_22_23.ps #(done with usearch haplotype identified data not considering the obitools3 taxonomic identification as well.)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 8 taxa and 24 samples ]
#sample_data() Sample Data:       [ 24 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 8 taxa by 6 taxonomic ranks ]



#write.csv(clean.otu.df, "aggregated_combined_GAL21_Diet_ready_for_phyloseq_otus_no_empty_pcrs.csv")
#write.csv(samples.df, "aggregated_combined_GAL21_Diet_ready_for_phyloseq_samples.csv")
#write.csv(taxa.df, "aggregated_combined_GAL21_Diet_ready_for_phyloseq_taxa.csv")



plot_richness(Climb_22_23.ps, x="Date","Site", measures=c("Shannon", "Simpson"), color="Site")

#ggsave("diet_alpha_diversity_by_species_and_site.jpg")




plot_richness(Climb_22_23.ps, x="Site","Associated_Blood_Samples", measures=c("Shannon", "Simpson"), color="Site")


# Transform data to proportions as appropriate for Bray-Curtis distances
Climb_22_23.ps.prop <- transform_sample_counts(Climb_22_23.ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(Climb_22_23.ps.prop, method="NMDS", distance="bray")
#Error in if (autotransform && xam > 50) { : 
#missing value where TRUE/FALSE needed

plot_ordination(Climb_22_23.ps.prop, ord.nmds.bray, color="annot_hap"
                , title="Bray NMDS")

Climb_22_23.ps.rclr<- tax_transform(Climb_22_23.ps, trans = "rclr")
Climb_22_23.ps.bin <- tax_transform(Climb_22_23.ps, trans = "binary")



#plot proportions of haplotype reads by number of ASVs by Site
plot_bar(Climb_22_23.ps.prop, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/22-23_Climb2_amplified/06_Phyloseq_viz/Bar_Chart_of_22-23_eDNA_Haplotype_proportions_of_ASVs_by_Site.jpg")

#plot haplotype reads by number of ASVs by Site
plot_bar(Climb_22_23.ps, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/22-23_Climb2_amplified/06_Phyloseq_viz/Bar_Chart_of_22-23_eDNA_Haplotype_reads_of ASVs_by_Site.jpg")


#plot haplotype ASV presence by Site (not aggregated by haplotype)
plot_bar(Climb_22_23.ps.bin, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")
ggsave("data/22-23_Climb2_amplified/06_Phyloseq_viz/Bar_Chart_of_22-23_eDNA_Haplotype_ASV_presence_by_Site.jpg")

###############

library(RColorBrewer)

plot_bar(ps.top20, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")+
  scale_fill_brewer(palette = "Paired")+
  scale_color_brewer(palette = "Paired")

#Normalize number of reads in each sample using median sequencing depth.

#I thought this would help, but now everything is really close to zero.

total = median(sample_sums(Climb_22_23.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
Climb_22_23.ps.med.norm = transform_sample_counts(Climb_22_23.ps, standf)

plot_bar(Climb_22_23.ps.med.norm, x="Site",fill="annot_hap")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")

#ggsave("Bar_Chart_of_Species_by_PreyClass.jpg")
#######Now do the same for the MF_GENUS phyloseq object


plot_bar(Climb_22_23.ps.prop, fill = "annot_hap", facet_grid = "Site")+ 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")

#ggsave("Faceted_sample_Site_proportions_annot_hap.jpg")



#This is the most important part to me to make the otu table readable.
#merge motus by annot_hap  
Climb_22_23_hap_merge.ps <- tax_glom(Climb_22_23.ps, taxrank = "annot_hap", NArm = FALSE)
plot_bar(Climb_22_23_hap_merge.ps, fill = "annot_hap") + 
  geom_bar(aes(color=annot_hap, fill=annot_hap), stat="identity", position="stack")

#DIDN'T WORK


Climb_22_23_hap_merge.ps
sample_variables(Climb_22_23_hap_merge.ps)
rank_names(Climb_22_23_hap_merge.ps)
tax_table(Climb_22_23_hap_merge.ps)
