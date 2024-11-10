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
load('data/21_Climb2/05_MetabaR_filtered/Climb_21_Metabarlist.RData')


### following tutorial on https://metabarfactory.github.io/metabaR/articles/metabaRF-vignette.html ###

#initial plot of number of reads in each category of samples and controls

# Create an input table (named check1) for ggplot of 3 columns:
#  (i) control type
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus,
#  (iii) the corresponding values.

check1 <- melt(Climb_21$pcrs[,c("control_type", "nb_reads", "nb_motus")])

ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) +
  geom_boxplot() + theme_bw() +
  geom_jitter(alpha=0.2) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") +
  theme(axis.text.x = element_text(angle=45, h=1))

# Using the nb_reads and nb_motus defined previously in the Climb_21$pcrs table

ggplot(Climb_21$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) +
  geom_point() + theme_bw() +
  scale_y_log10() + scale_x_log10() +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")

#lack of correlation shows that the diversity is well covered by the sequencing depth.
#It looks like we need more sequencing depth to me.

##### identifying contaminants #####

# Identifying extraction contaminants
Climb_21 <- contaslayer(Climb_21,
                       control_types = "extraction",
                       output_col = "not_an_extraction_conta")

table(Climb_21$motus$not_an_extraction_conta)

## Identifying pcr contaminants
Climb_21 <- contaslayer(Climb_21,
                       control_types = "pcr",
                       output_col = "not_a_pcr_conta")

table(Climb_21$motus$not_a_pcr_conta)

## Identifying "sequencing" contaminants -I called the field negative a sequencing control since I didn't use a sequencing control in this set and metabar doesn't have a field negative option.
Climb_21 <- contaslayer(Climb_21,
                           control_types = "sequencing",
                           output_col = "not_a_field_neg_conta")

table(Climb_21$motus$not_a_field_neg_conta)

# Compute relative abundance of all extraction contaminants together
a <- data.frame(conta.relab = rowSums(Climb_21$reads[,!Climb_21$motus$not_an_extraction_conta]) /
                  rowSums(Climb_21$reads))
# Add information on control types
a$control_type <- Climb_21$pcrs$control_type[match(rownames(a), rownames(Climb_21$pcrs))]

ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) +
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") +
  theme_bw() +
  scale_y_log10()

# flag pcrs with total extraction contaminant relative abundance > 10% of reads)
Climb_21$pcrs$low_contamination_level <-
  ifelse(a$conta.relab[match(rownames(Climb_21$pcrs), rownames(a))]>1e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(Climb_21$pcrs$low_contamination_level) / nrow(Climb_21$pcrs)

#Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa according to obitools3 assignment to all vertebrate database
Climb_21$motus$target_taxon <- grepl("TRUE", Climb_21$motus$ID_STATUS)

# Proportion of each of these over total number of MOTUs
table(Climb_21$motus$target_taxon) / nrow(Climb_21$motus)

# Intersection with extraction contaminant flags (not contaminant = T)
table(Climb_21$motus$target_taxon,
      Climb_21$motus$not_an_extraction_conta)

table(Climb_21$motus$target_taxon,
      Climb_21$motus$not_a_field_neg_conta)
# Plot the unweighted distribution of MOTUs similarity scores
a <-
  ggplot(Climb_21$motus, aes(x=Climb_21$motus$BEST_IDENTITY)) +
  geom_histogram(color="grey", fill="white", bins=200) +
  geom_vline(xintercept = 0.98, col="orange", lty=2) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x="% similarity against best match", y="# MOTUs")

# Same for the weighted distribution
b <-
  ggplot(Climb_21$motus,
         aes(x=Climb_21$motus$BEST_IDENTITY, y = ..count.., weight = Climb_21$motus$count)) +
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
Climb_21$motus$not_degraded <-
  ifelse(Climb_21$motus$BEST_IDENTITY < 0.98, F, T)

# Proportion of each of these over total number of MOTUs
table(Climb_21$motus$not_degraded) / nrow(Climb_21$motus)

# Flag poorly-matched haplotypes (TRUE) vs. well-matched haplotypes (FALSE)
Climb_21$motus$hap_match <-
<<<<<<< Updated upstream
  ifelse(Climb_21$motus$annot_pctid < 95, F, T)
=======
  ifelse(Climb_21$motus$annot_pctid <100, F, T)

Climb_21$motus$hap_length <-
  ifelse(Climb_21$motus$seq_len <393, F, T)
>>>>>>> Stashed changes

# Proportion of each of these over total number of MOTUs
table(Climb_21$motus$hap_match) / nrow(Climb_21$motus)



# Intersection with other flags
table(Climb_21$motus$target_taxon,
      Climb_21$motus$not_an_extraction_conta,
      Climb_21$motus$not_a_pcr_conta,
      Climb_21$motus$not_a_field_neg_conta,
      Climb_21$motus$not_degraded,
      Climb_21$motus$hap_match)

##### detecting PCR outliers #####

ggplot(Climb_21$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") +
  geom_vline(xintercept = 1e02, lty=2, color="orange") + # threshold
  scale_x_log10() +
  labs(x="# Reads (with all MOTUs and PCRs)",
       y="# PCRs") +
  theme_bw() +
  theme(panel.grid = element_blank())

###### subset the metabarlist based on criteria so far and look for PCR outliers again ######

Climb_21_no_contas <- subset_metabarlist(Climb_21, "motus", 
                                        indices = rowSums(Climb_21$motus[,c("not_an_extraction_conta", "not_a_pcr_conta", "not_a_field_neg_conta")]) == 3)

summary_metabarlist(Climb_21_no_contas)
#very important after a subsetting

# Compute the number of reads per pcr
Climb_21_no_contas$pcrs$nb_reads <- rowSums(Climb_21_no_contas$reads)

# Compute the number of motus per pcr
Climb_21_no_contas$pcrs$nb_motus <- rowSums(Climb_21_no_contas$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_21_no_contas <- subset_metabarlist(Climb_21_no_contas, table = "pcrs",
                                        indices = Climb_21_no_contas$pcrs$nb_reads>0)

summary_metabarlist(Climb_21)
summary_metabarlist(Climb_21_no_contas)

ggplot(Climb_21$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") +
  geom_vline(xintercept = 1e02, lty=2, color="orange") + # threshold
  scale_x_log10() +
  labs(x="# Reads (with all MOTUs and PCRs)",
       y="# PCRs") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
Climb_21_no_contas$pcrs$seqdepth_ok <- ifelse(Climb_21_no_contas$pcrs$nb_reads < 1e02, F, T)

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE) in full metabarlist too
Climb_21$pcrs$seqdepth_ok <- ifelse(Climb_21$pcrs$nb_reads < 1e02, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(Climb_21_no_contas$pcrs$seqdepth_ok[Climb_21_no_contas$pcrs$type=="sample"]) /
  nrow(Climb_21_no_contas$pcrs[Climb_21_no_contas$pcrs$type=="sample",])

### Lowering tag-jumps ####

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,5e-4,1e-3, 1e-2, 3e-2, 5e-2)

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(Climb_21,x))
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
tmp$controls <- Climb_21$pcrs$control_type[match(tmp$sample, rownames(Climb_21$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) +
  geom_boxplot(color="grey40") +
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.05"), col="orange", lty=2) +
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


Climb_21_hap_match <- subset_metabarlist(Climb_21_no_contas, "motus", 
<<<<<<< Updated upstream
                                         indices = Climb_21_no_contas$motus$annot_pctid > 95)
=======
                                         indices = Climb_21_no_contas$motus$seq_len >393)
>>>>>>> Stashed changes

summary_metabarlist(Climb_21_hap_match)
#very important after a subsetting

# Compute the number of reads per pcr
Climb_21_hap_match$pcrs$nb_reads <- rowSums(Climb_21_hap_match$reads)

# Compute the number of motus per pcr
Climb_21_hap_match$pcrs$nb_motus <- rowSums(Climb_21_hap_match$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_21_hap_match <- subset_metabarlist(Climb_21_hap_match, table = "pcrs",
                                         indices = Climb_21_hap_match$pcrs$nb_reads>0)

summary_metabarlist(Climb_21)
summary_metabarlist(Climb_21_no_contas)
summary_metabarlist(Climb_21_hap_match)

#looks like 0.05 is where the extraction and field negative controls drop off.


###### Pie charts of the noise in the dataset #####

# Create a table of MOTUs quality criteria
# noise is identified as FALSE in Climb_21, the "!" transforms it to TRUE
motus.qual <- !Climb_21$motus[,c("not_an_extraction_conta", "target_taxon", "not_degraded", "hap_match", "hap_length")]
colnames(motus.qual) <- c("extraction_conta", "untargeted_taxon", "degraded_seq", "no_hap_match","short_hap")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(Climb_21$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(Climb_21$motus$count[x])/sum(Climb_21$motus$count))

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


# noise is identified as FALSE in Climb_21, the "!" transforms it to TRUE
pcrs.qual <- !Climb_21$pcrs[,c("low_contamination_level", "seqdepth_ok")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[Climb_21$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[Climb_21$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[Climb_21$pcrs$type=="sample",])

tmp.pcrs <-
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[Climb_21$pcrs$type=="sample",x]==T,
           colnames(pcrs.qual)[x], NA)}), 1, function(x) {
             paste(sort(unique(x)), collapse = "|")
           })
tmp.pcrs <- as.data.frame(gsub("^$", "not_artefactual", tmp.pcrs))

colnames(tmp.pcrs) <- "artefact_type"

ggplot(# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
# Use tag-jump corrected metabarlist with the threshold identified above
tmp <- tests[["t_0.05"]]tmp <- tests[["t_0.05"]]_hap.pcrs, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") +
  coord_polar(theta="y") + theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.direction = "vertical") +
  ggtitle("PCR artefacts overview")

##### Data cleaning and Aggregation #####

# Use tag-jump corrected metabarlist with the threshold identified above
<<<<<<< Updated upstream
tmp <- tests[["t_0.05"]]

=======
tmp <- tests[["t_0.03"]]
#tmp <- tests[["t_0"]]
>>>>>>> Stashed changes
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
Climb_21_clean <- subset_metabarlist(tmp, "pcrs",
                                    indices = rowSums(tmp$pcrs[,c("low_contamination_level",
                                                                  "seqdepth_ok")]) == 2 &
                                      tmp$pcrs$type == "sample")
summary_metabarlist(Climb_21_clean)

if(sum(colSums(Climb_21_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(Climb_21_clean$reads)==0)>0){print("empty pcrs present")}

Climb_21_clean$motus$count = colSums(Climb_21_clean$reads)
Climb_21_clean$pcrs$nb_reads_postmetabaR = rowSums(Climb_21_clean$reads)
Climb_21_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(Climb_21_clean$reads>0, T, F))

#look at abundance and richness before and after metabar filtering

check <- melt(Climb_21_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR",
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
<<<<<<< Updated upstream
tmp_hap <- tests[["t_0.05"]]
=======
tmp_hap <- tests[["t_0.03"]]
#tmp_hap <- tests[["t_0"]]
>>>>>>> Stashed changes

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
Climb_21_hap_matched_clean <- subset_metabarlist(tmp_hap, "pcrs",
                                     indices = rowSums(tmp_hap$pcrs[,c("low_contamination_level",
                                                                   "seqdepth_ok")]) == 2 &
                                       tmp_hap$pcrs$type == "sample")
summary_metabarlist(Climb_21_hap_matched_clean)

if(sum(colSums(Climb_21_hap_matched_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(Climb_21_hap_matched_clean$reads)==0)>0){print("empty pcrs present")}

Climb_21_hap_matched_clean$motus$count = colSums(Climb_21_hap_matched_clean$reads)
Climb_21_hap_matched_clean$pcrs$nb_reads_postmetabaR = rowSums(Climb_21_hap_matched_clean$reads)
Climb_21_hap_matched_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(Climb_21_hap_matched_clean$reads>0, T, F))

#look at abundance and richness before and after metabar filtering

check <- melt(Climb_21_hap_matched_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR",
                                     "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))


summary_metabarlist(Climb_21_hap_matched_clean)
summary_metabarlist(Climb_21_clean)



##### Data aggregation ##### 
# can comment this out because the three PCR replicates per sample were already pooled before sequencing, but keeping it doesn't hurt anything.
Climb_21_agg <- aggregate_pcrs(Climb_21_hap_matched_clean)
summary_metabarlist(Climb_21_agg)
#very important after a subsetting

# Compute the number of reads per pcr
Climb_21_agg$pcrs$nb_reads <- rowSums(Climb_21_agg$reads)

# Compute the number of motus per pcr
Climb_21_agg$pcrs$nb_motus <- rowSums(Climb_21_agg$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_21_species <- subset_metabarlist(Climb_21_clean, table = "pcrs",
                                    indices = Climb_21_clean$pcrs$nb_reads>0)

summary_metabarlist(Climb_21_species)

#subset the metabarlist to all PCRs with more than 0 reads
Climb_21_final <- subset_metabarlist(Climb_21_agg, table = "pcrs",
                                       indices = Climb_21_agg$pcrs$nb_reads>0)

summary_metabarlist(Climb_21_final)

