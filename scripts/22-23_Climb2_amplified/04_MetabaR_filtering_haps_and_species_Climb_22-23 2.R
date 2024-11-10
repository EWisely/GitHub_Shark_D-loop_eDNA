### MetabaR filtering of results from D-loop primers ####
#Author: Eldridge Wisely
#Date: 10-17-2023

### Load packages ###
library("metabaR")
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
# Load requested package for plotting
library("reshape2")
library("tidyverse")
library("ggplot2") 

#### load initial metabarlist from previous script ####
load('data/22-23_Climb2_amplified/05_MetabaR_filtered/Climb_22_23_Metabarlist.RData')


### following tutorial on https://metabarfactory.github.io/metabaR/articles/metabaRF-vignette.html ###

#initial plot of number of reads in each category of samples and controls

# Create an input table (named check1) for ggplot of 3 columns:
#  (i) control type
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus,
#  (iii) the corresponding values.

check1 <- melt(Climb_22_23$pcrs[,c("control_type", "nb_reads", "nb_motus")])

ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) +
  geom_boxplot() + theme_bw() +
  geom_jitter(alpha=0.2) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") +
  theme(axis.text.x = element_text(angle=45, h=1))

# Using the nb_reads and nb_motus defined previously in the Climb_22_23$pcrs table

ggplot(Climb_22_23$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) +
  geom_point() + theme_bw() +
  scale_y_log10() + scale_x_log10() +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")

#lack of correlation shows that the diversity is well covered by the sequencing depth.
#It looks like we need more sequencing depth to me.

##### identifying contaminants #####

# Identifying extraction contaminants
Climb_22_23 <- contaslayer(Climb_22_23,
                        control_types = "extraction",
                        output_col = "not_an_extraction_conta")

table(Climb_22_23$motus$not_an_extraction_conta)

## Identifying pcr contaminants
Climb_22_23 <- contaslayer(Climb_22_23,
                        control_types = "pcr",
                        output_col = "not_a_pcr_conta")

table(Climb_22_23$motus$not_a_pcr_conta)

## Identifying "sequencing" contaminants -I called the field negative a sequencing control since I didn't use a sequencing control in this set and metabar doesn't have a field negative option.
Climb_22_23 <- contaslayer(Climb_22_23,
                        control_types = "sequencing",
                        output_col = "not_a_field_neg_conta")

table(Climb_22_23$motus$not_a_field_neg_conta)

# Compute relative abundance of all extraction contaminants together
a <- data.frame(conta.relab = rowSums(Climb_22_23$reads[,!Climb_22_23$motus$not_a_field_neg_conta]) /
                  rowSums(Climb_22_23$reads))
# Add information on control types
a$control_type <- Climb_22_23$pcrs$control_type[match(rownames(a), rownames(Climb_22_23$pcrs))]

ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) +
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") +
  theme_bw() +
  scale_y_log10()

# flag pcrs with total extraction contaminant relative abundance > 10% of reads)
Climb_22_23$pcrs$low_contamination_level <-
  ifelse(a$conta.relab[match(rownames(Climb_22_23$pcrs), rownames(a))]>1e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(Climb_22_23$pcrs$low_contamination_level) / nrow(Climb_22_23$pcrs)

#Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa according to obitools3 assignment to all vertebrate database
Climb_22_23$motus$target_taxon <- grepl("TRUE", Climb_22_23$motus$ID_STATUS)

# Proportion of each of these over total number of MOTUs
table(Climb_22_23$motus$target_taxon) / nrow(Climb_22_23$motus)

# Intersection with extraction contaminant flags (not contaminant = T)
table(Climb_22_23$motus$target_taxon,
      Climb_22_23$motus$not_an_extraction_conta)

table(Climb_22_23$motus$target_taxon,
      Climb_22_23$motus$not_a_field_neg_conta)
# Plot the unweighted distribution of MOTUs similarity scores
a <-
  ggplot(Climb_22_23$motus, aes(x=Climb_22_23$motus$BEST_IDENTITY)) +
  geom_histogram(color="grey", fill="white", bins=200) +
  geom_vline(xintercept = 0.98, col="orange", lty=2) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x="% similarity against best match", y="# MOTUs")

# Same for the weighted distribution
b <-
  ggplot(Climb_22_23$motus,
         aes(x=Climb_22_23$motus$BEST_IDENTITY, y = ..count.., weight = Climb_22_23$motus$count)) +
  geom_histogram(color="grey", fill="white", bins=200) +
  geom_vline(xintercept = 0.98, col="orange", lty=2) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x="% similarity against best match", y="# Reads")

# Combine plots into one
library(cowplot)
ggdraw() +
  draw_plot(a, x=0, y=0, width = 0.5) +
  draw_plot(b, x=0.5, y=0, width = 0.5)


##### identifying degraded sequences #####

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
Climb_22_23$motus$not_degraded <-
  ifelse(Climb_22_23$motus$BEST_IDENTITY < 0.98, F, T)

# Proportion of each of these over total number of MOTUs
table(Climb_22_23$motus$not_degraded) / nrow(Climb_22_23$motus)

# Flag poorly-matched haplotypes (TRUE) vs. well-matched haplotypes (FALSE)
Climb_22_23$motus$hap_match <-
  ifelse(Climb_22_23$motus$annot_pctid <100, F, T)

Climb_22_23$motus$hap_length <-
  ifelse(Climb_22_23$motus$seq_len <393, F, T)

# Proportion of each of these over total number of MOTUs
table(Climb_22_23$motus$hap_match) / nrow(Climb_22_23$motus)



# Intersection with other flags
table(Climb_22_23$motus$target_taxon,
      Climb_22_23$motus$not_an_extraction_conta,
      Climb_22_23$motus$not_a_pcr_conta,
      Climb_22_23$motus$not_a_field_neg_conta,
      Climb_22_23$motus$not_degraded,
      Climb_22_23$motus$hap_match)

##### detecting PCR outliers #####

ggplot(Climb_22_23$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") +
  geom_vline(xintercept = 1e01, lty=2, color="orange") + # threshold
  scale_x_log10() +
  labs(x="# Reads (with all MOTUs and PCRs)",
       y="# PCRs") +
  theme_bw() +
  theme(panel.grid = element_blank())

###### subset the metabarlist based on criteria so far and look for PCR outliers again ######

Climb_22_23_no_contas <- subset_metabarlist(Climb_22_23, "motus", 
                                         indices = rowSums(Climb_22_23$motus[,c("not_an_extraction_conta", "not_a_pcr_conta", "not_a_field_neg_conta")]) == 3)

summary_metabarlist(Climb_22_23_no_contas)
#very important after a subsetting

# Compute the number of reads per pcr
Climb_22_23_no_contas$pcrs$nb_reads <- rowSums(Climb_22_23_no_contas$reads)

# Compute the number of motus per pcr
Climb_22_23_no_contas$pcrs$nb_motus <- rowSums(Climb_22_23_no_contas$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_22_23_no_contas <- subset_metabarlist(Climb_22_23_no_contas, table = "pcrs",
                                         indices = Climb_22_23_no_contas$pcrs$nb_reads>0)

summary_metabarlist(Climb_22_23)
summary_metabarlist(Climb_22_23_no_contas)

ggplot(Climb_22_23$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") +
  geom_vline(xintercept = 1e01, lty=2, color="orange") + # threshold
  scale_x_log10() +
  labs(x="# Reads (with all MOTUs and PCRs)",
       y="# PCRs") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
Climb_22_23_no_contas$pcrs$seqdepth_ok <- ifelse(Climb_22_23_no_contas$pcrs$nb_reads < 1e01, F, T)

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE) in full metabarlist too
Climb_22_23$pcrs$seqdepth_ok <- ifelse(Climb_22_23$pcrs$nb_reads < 1e01, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(Climb_22_23_no_contas$pcrs$seqdepth_ok[Climb_22_23_no_contas$pcrs$type=="sample"]) /
  nrow(Climb_22_23_no_contas$pcrs[Climb_22_23_no_contas$pcrs$type=="sample",])

### Lowering tag-jumps ####

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,5e-4,1e-3, 1e-2, 3e-2, 5e-2)

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(Climb_22_23,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- Climb_22_23$pcrs$control_type[match(tmp$sample, rownames(Climb_22_23$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) +
  geom_boxplot(color="grey40") +
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.03"), col="orange", lty=2) +
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=5) +
  theme_bw() +
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=40, h=1),
        legend.position = "none")
####Subset the no contas metabarlist to just annotated ASVs ####


Climb_22_23_hap_match <- subset_metabarlist(Climb_22_23_no_contas, "motus", 
                                         indices = Climb_22_23_no_contas$motus$seq_len >393)

summary_metabarlist(Climb_22_23_hap_match)
#very important after a subsetting

# Compute the number of reads per pcr
Climb_22_23_hap_match$pcrs$nb_reads <- rowSums(Climb_22_23_hap_match$reads)

# Compute the number of motus per pcr
Climb_22_23_hap_match$pcrs$nb_motus <- rowSums(Climb_22_23_hap_match$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_22_23_hap_match <- subset_metabarlist(Climb_22_23_hap_match, table = "pcrs",
                                         indices = Climb_22_23_hap_match$pcrs$nb_reads>0)

summary_metabarlist(Climb_22_23)
summary_metabarlist(Climb_22_23_no_contas)
summary_metabarlist(Climb_22_23_hap_match)

#looks like 0.03 is where the extraction and field negative controls drop off.


###### Pie charts of the noise in the dataset #####

# Create a table of MOTUs quality criteria
# noise is identified as FALSE in Climb_22_23, the "!" transforms it to TRUE
motus.qual <- !Climb_22_23$motus[,c("not_an_extraction_conta", "target_taxon", "not_degraded", "hap_match", "hap_length")]
colnames(motus.qual) <- c("extraction_conta", "untargeted_taxon", "degraded_seq", "no_hap_match","short_hap")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(Climb_22_23$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(Climb_22_23$motus$count[x])/sum(Climb_22_23$motus$count))

tmp.motus <-
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "not_artefactual", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

ggplot(tmp.motus, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") +
  coord_polar(theta="y") + theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.direction = "vertical") +
  ggtitle("MOTUs artefacts overview")

# Create a table of pcrs quality criteria


# noise is identified as FALSE in Climb_22_23, the "!" transforms it to TRUE
pcrs.qual <- !Climb_22_23$pcrs[,c("low_contamination_level", "seqdepth_ok")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[Climb_22_23$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[Climb_22_23$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[Climb_22_23$pcrs$type=="sample",])

tmp.pcrs <-
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[Climb_22_23$pcrs$type=="sample",x]==T,
           colnames(pcrs.qual)[x], NA)}), 1, function(x) {
             paste(sort(unique(x)), collapse = "|")
           })
tmp.pcrs <- as.data.frame(gsub("^$", "not_artefactual", tmp.pcrs))

colnames(tmp.pcrs) <- "artefact_type"

# Use tag-jump corrected metabarlist with the threshold identified above

##### Data cleaning and Aggregation #####

# Use tag-jump corrected metabarlist with the threshold identified above
tmp <- tests[["t_0.03"]]
#tmp <- tests[["t_0"]]
#The line above removes filtering based on putative tag-jumps since the diversity and abundance characteristics that filter is based on were just off-target reads.

# Subset on MOTUs: we keep motus that are defined as TRUE following the
# three criteria below (sum of three TRUE is equal to 3 with the rowSums function)
tmp <- subset_metabarlist(tmp, "motus",
                          indices = rowSums(tmp$motus[,c("not_an_extraction_conta", 
                                                         "target_taxon",#this is based on obitools assignment and we're going with usearch annotation
                                                         "not_a_pcr_conta",
                                                         "not_degraded", #this is based on obitools assignment and we're going with usearch annotation
                                                         "not_a_field_neg_conta"
                                                         #,"hap_match"
                          )]) == 5)
summary_metabarlist(tmp)

# Subset on pcrs and keep only samples, no controls
# Climb_22_23_clean <- subset_metabarlist(tmp, "pcrs",
#                                     indices = rowSums(tmp$pcrs[,c("low_contamination_level",
#                                                                   "seqdepth_ok")]) == 2 &
#                                       tmp$pcrs$type == "sample")

#subset without sequencing depth filter
Climb_22_23_clean <- subset_metabarlist(tmp, "pcrs",
                                     indices = tmp$pcrs$low_contamination_level ==TRUE &
                                       tmp$pcrs$type == "sample")
summary_metabarlist(Climb_22_23_clean)


if(sum(colSums(Climb_22_23_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(Climb_22_23_clean$reads)==0)>0){print("empty pcrs present")}

Climb_22_23_clean$motus$count = colSums(Climb_22_23_clean$reads)
Climb_22_23_clean$pcrs$nb_reads_postmetabaR = rowSums(Climb_22_23_clean$reads)
Climb_22_23_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(Climb_22_23_clean$reads>0, T, F))

#look at abundance and richness before and after metabar filtering

check <- melt(Climb_22_23_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR",
                                     "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))

#### Make a cleaned file that only requires usearch annotation (and not requiring obitools percent ID or assignment)####

# Use tag-jump corrected metabarlist with the threshold identified above
tmp_hap <- tests[["t_0.03"]]
#tmp_hap <- tests[["t_0"]]

#The line above removes filtering based on putative tag-jumps since the diversity and abundance characteristics that filter is based on were just off-target reads.

# Subset on MOTUs: we keep motus that are defined as TRUE following the
# three criteria below (sum of three TRUE is equal to 3 with the rowSums function)
tmp_hap <- subset_metabarlist(tmp_hap, "motus",
                              indices = rowSums(tmp_hap$motus[,c("not_an_extraction_conta", 
                                                                 #"target_taxon",#this is based on obitools assignment and we're going with usearch annotation
                                                                 "not_a_pcr_conta",
                                                                 #"not_degraded", #this is based on obitools assignment and we're going with usearch annotation
                                                                 "not_a_field_neg_conta",
                                                                 "hap_match",
                                                                 "hap_length")]) == 5)
summary_metabarlist(tmp_hap)

# Subset on pcrs and keep only samples, no controls
Climb_22_23_hap_matched_clean <- subset_metabarlist(tmp_hap, "pcrs",
                                                 indices = rowSums(tmp_hap$pcrs[,c("low_contamination_level",
                                                                                   "seqdepth_ok")]) == 2 &
                                                   tmp_hap$pcrs$type == "sample")
summary_metabarlist(Climb_22_23_hap_matched_clean)

#subset without sequencing depth filter
Climb_22_23_hap_matched_clean <- subset_metabarlist(tmp_hap, "pcrs",
                                                 indices = tmp$pcrs$low_contamination_level ==TRUE &
                                                   tmp$pcrs$type == "sample")
summary_metabarlist(Climb_22_23_hap_matched_clean)



if(sum(colSums(Climb_22_23_hap_matched_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(Climb_22_23_hap_matched_clean$reads)==0)>0){print("empty pcrs present")}

Climb_22_23_hap_matched_clean$motus$count = colSums(Climb_22_23_hap_matched_clean$reads)
Climb_22_23_hap_matched_clean$pcrs$nb_reads_postmetabaR = rowSums(Climb_22_23_hap_matched_clean$reads)
Climb_22_23_hap_matched_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(Climb_22_23_hap_matched_clean$reads>0, T, F))

#look at abundance and richness before and after metabar filtering

check <- melt(Climb_22_23_hap_matched_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR",
                                                 "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))


summary_metabarlist(Climb_22_23_hap_matched_clean)
summary_metabarlist(Climb_22_23_clean)



##### Data aggregation ##### 
# can comment this out because the three PCR replicates per sample were already pooled before sequencing, but keeping it doesn't hurt anything.
Climb_22_23_agg <- aggregate_motus(Climb_22_23_hap_matched_clean, groups = Climb_22_23_hap_matched_clean$motus$annot_hap)
summary_metabarlist(Climb_22_23_agg)
#very important after a subsetting

# Compute the number of reads per pcr
Climb_22_23_agg$pcrs$nb_reads <- rowSums(Climb_22_23_agg$reads)

# Compute the number of motus per pcr
Climb_22_23_agg$pcrs$nb_motus <- rowSums(Climb_22_23_agg$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_22_23_species <- subset_metabarlist(Climb_22_23_clean, table = "pcrs",
                                       indices = Climb_22_23_clean$pcrs$nb_reads>0)

summary_metabarlist(Climb_22_23_species)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_22_23_final <- subset_metabarlist(Climb_22_23_agg, table = "pcrs",
                                     indices = Climb_22_23_agg$pcrs$nb_reads>0)

summary_metabarlist(Climb_22_23_final)

#### Make a csv that substitutes motus$SCIENTIFIC_NAME for motu_id and puts it as colnames with pcr$rownames and read_numbers from Climb_22_23_final$reads ####
#Climb_22_23_final$pcrs$species <-
#  Climb_22_23_final$samples$`Species (scientific name)`[match(Climb_22_23_final$pcrs$sample_id,
#                                                                rownames(Climb_22_23_final$samples))]

df<-as.data.frame.array(Climb_22_23_final$reads)
rownames(df)
colnames(df)

library(data.table)
setnames(df, as.character(rownames(Climb_22_23_final$motus)), as.character(Climb_22_23_final$motus$annot_hap))


head(df)
library(tibble)
library(dplyr)
Climb_22_23_final_pcrs<-as.data.frame(Climb_22_23_final$pcrs)
Climb_22_23_final_pcrs<-tibble::rownames_to_column(.data = Climb_22_23_final_pcrs, var="SampleID")

Climb_22_23_final_samples<-as.data.frame(Climb_22_23_final$samples)
Climb_22_23_final_samples<-tibble::rownames_to_column(.data = Climb_22_23_final_samples, var="sample_id")


df1<-tibble::rownames_to_column(.data = df, var="SampleID")
df2<-left_join(df1,Climb_22_23_final_pcrs, by="SampleID")
df3<-left_join(df2,Climb_22_23_final_samples, by="sample_id")

df5<-column_to_rownames(.data = df3, var = "SampleID")

write.csv(df5, "data/blood_vs_water_haps/2022-2023_Climb_aggregated_haplotypes_with_sampling_data.csv")

