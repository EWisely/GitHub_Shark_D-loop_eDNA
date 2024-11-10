library(vegan)
library(tidyverse)

blood_and_water<-read.csv("data/blood_vs_water_haps/2021_Climb_aggregated_blood_and_water_haplotypes_with_sampling_data.csv")
blood_and_water<-blood_and_water%>% column_to_rownames(var="X")
blood_and_water_counts<-blood_and_water[1:2]
blood_and_water_groups<-blood_and_water[c(3:4, 6:7, 12:13, 21:23)]

#group by Site_Date and material, then calculate total proportions of each haplotype
prop_blood_and_water<-blood_and_water%>%
  mutate(Percent_Hap20 = Climb_Hap_20 /(Climb_Hap_20+Climb_Hap_22))%>%
  mutate(Percent_Hap22 =Climb_Hap_22 /(Climb_Hap_20+Climb_Hap_22))%>%
  group_by(Site_Date,material)%>%
  summarize(mean_prop_hap20= mean(Percent_Hap20),
            mean_prop_hap22 = mean(Percent_Hap22))%>%
  ungroup()

prop_blood_and_water

example_blood_v_water<-blood_and_water[c(1:4,21)]



# bar plot with Labels Inside bars
ggplot(data=blood_and_water, aes(x=sample_id, y=nb_motus, color =material)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()

# Stacked barplot with multiple groups
ggplot(data=blood_and_water, aes(x=Site, y=nb_motus, fill=material)) +
  geom_bar(stat="identity")

# Use position=position_dodge()
ggplot(data=blood_and_water, aes(x=Site, y=Percent_Hap20, fill=material)) +
  geom_bar(stat="identity", position=position_dodge())

#ggsave("data/21_Climb2/06_Phyloseq_viz/number_of_motus_in_blood_vs_water_by_site.jpg")

# Use position=position_dodge()
p<-ggplot(data=blood_and_water, aes(x=Site_Date, y=Percent_Hap20, fill=material)) +
  geom_bar(stat="identity", position=position_dodge())
p + theme(axis.text.x = element_text(angle = -90))

ggsave("data/21_Climb2/06_Phyloseq_viz/number_of_motus_in_blood_vs_water_by_site_date.jpg")

#convert data to percentages instead of raw read counts/shark counts

#following information at https://www.rdocumentation.org/packages/vegan/versions/1.11-0/topics/decostand
data(varespec)
sptrans <- decostand(varespec, "max")
apply(sptrans, 2, max)
sptrans <- wisconsin(varespec)

sptrans <-decostand(blood_and_water_counts, "max")
apply(sptrans, 2, max)
sptrans <- wisconsin(blood_and_water_counts)

sptrans <- decostand(blood_and_water_counts, "chi.square")
plot(procrustes(rda(sptrans), cca(blood_and_water_counts)))

#following tutorial at https://rpubs.com/CPEL/NMDS
blood_and_water_prop <- decostand(blood_and_water_counts, method = "total")

blood_and_water_prop_distmat <-vegdist(blood_and_water_prop, method = "bray")

blood_and_water_prop_distmat <- 
  as.matrix(blood_and_water_prop_distmat, labels = T)
write.csv(blood_and_water_prop_distmat, "data/blood_vs_water_haps/blood_and_water_prop_distmat.csv")

# Running NMDS in vegan (metaMDS)
blood_and_water_NMS <-
  metaMDS(blood_and_water_prop_distmat,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

# Shepards test/goodness of fit
goodness(blood_and_water_NMS) # Produces a results of test statistics for goodness of fit for each point

stressplot(blood_and_water_NMS) # Produces a Shepards diagram

# Plotting points in ordination space
plot(blood_and_water_NMS, "sites")   # Produces distance #here sites=sampleIDs
orditorp(blood_and_water_NMS, "sites")   # Gives points labels

colvec <- c("darkred", "blue")   # Identifies colors for group assignments
pchvec <- c(21, 22)   # Identifies character symbols for group assignments


plot(blood_and_water_NMS)
with(blood_and_water_groups,
     points(blood_and_water_NMS,
            display = "sites",
            col = "black",
            pch = pchvec[material],
            bg = colvec[material]))

ordihull(
  blood_and_water_NMS,
  blood_and_water_groups$material,
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  border = c("red", "aquamarine"),
  lty = c(1, 2),
  lwd = 2.5
)

# Calculating centroids 

# You can calculate centroids for your groups which can be viewed as the average position of observations in ordination space. 

# Calculating and plotting centroids of NMDS Result
scrs <-
  scores(blood_and_water_NMS, display = "sites", "species")
cent <-
  aggregate(scrs ~ material, data = blood_and_water_groups, FUN = "mean")
names(cent) [-1] <- colnames(scrs)
points(cent [,-1],
       pch = c( 8 , 8 ),
       col = c("darkred", "darkblue"),
       bg = c("black"),
       lwd = 3.0,
       cex = 2.0 # Plots centroids as points on ordination
)

orditorp(blood_and_water_NMS, "sites")   # Gives points labels




blood_and_water.spp_scrs <- 
  sppscores(blood_and_water_NMS) <- blood_and_water_prop

plot(blood_and_water_NMS, "sites")   # Produces distance based bi-plot
plot(blood_and_water_NMS, "species")   # Plots species scores
orditorp(blood_and_water_NMS, "sites")   # Gives points labels

blood_and_water.spp_cor <-
  cor(blood_and_water_prop,
      blood_and_water_NMS$points,
      use = "complete.obs",
      method = "pearson")
write.csv(blood_and_water.spp_cor, file = "data/blood_vs_water_haps/blood_and_water_spp_PearsonCor.csv")

####plot the proportion data ####

#add the groups data back in...
blood_and_water_prop1<-rownames_to_column(blood_and_water_prop, var="SampleID")
blood_and_water_groups1<-rownames_to_column(blood_and_water_groups, var = "SampleID")
meta_blood_and_water_prop<-left_join(blood_and_water_prop1, blood_and_water_groups1, by="SampleID")

write.csv(meta_blood_and_water_prop, file ="data/blood_vs_water_haps/metadata_blood_and_water_proportions_21.csv")

# Stacked barplot with multiple groups
ggplot(data=meta_blood_and_water_prop, aes(x=Site, y=, fill=material)) +
  geom_bar(stat="identity", position=position_dodge())
# Use position=position_dodge()
ggplot(data=blood_and_water, aes(x=Site, y=nb_motus, fill=material)) +
  geom_bar(stat="identity", position=position_dodge())

#following tutorial at https://larmarange.github.io/ggstats/articles/stat_prop.html
library(ggstats)
library(ggplot2)

d <- as.data.frame(Titanic)
p <- ggplot(d) +
  aes(x = Class, fill = Survived, weight = Freq, by = Class) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))
p
p + facet_grid(cols = vars(Sex))
###work on this more
df <- as.data.frame(blood_and_water)
p <- ggplot(df) +
  aes(x = Climb2, fill = material, weight = Date, by = Site) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5))
p
p + facet_grid(cols = vars(Sex))




# vegan tutorial at https://rdrr.io/cran/vegan/f/inst/doc/diversity-vegan.pdf
data("BCI")
rad<-radfit(BCI)

rad_Climb_21<- rad.lognormal(as.data.frame(blood_and_water_counts), family = poisson)
rad.lognormal(blood_and_water_counts, family = poisson, ...)
#look into Deseq2 (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)






#Mantel test looks at correlation between two different distance matrices.

#following tutorial from https://www.youtube.com/watch?v=_lFO_BA5nUs

set.seed(1)

#Basid NMDS
nmds_result<-metaMDS(blood_and_water)
plot(nmds_result)

#stress plot
stressplot(nmds_result)
nmds_result$stress

