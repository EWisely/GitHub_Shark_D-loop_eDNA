library(ggplot2)
library(ggshadow)
library(scatterpie)
library(geodata)  # spatial data download
library(rnaturalearth)  # more spatial data download
library(terra)
library(tidyterra)
library(car)  # 'pointLabel' function
library(dplyr)
library(basemaps)
library(ggrepel)

##### Prepare the data #####

#load data to make scatterpies with
blood_haps<- read.csv("metadata/percentage_of_haps_in_bloods.csv")
blood_haps<-blood_haps%>% arrange(total_bt_sharks)
blood_haps<-blood_haps[-6, ]
blood_haps$radius <- .0017* blood_haps$total_bt_sharks



#make scatterpies
ggplot() + geom_scatterpie(aes(x=Long, y=Lat, group=Site), data=blood_haps,
                           cols=c("count_hap1", "count_hap2")) + coord_equal()



#plot blood haps on low resolution world map to test
world <- map_data('world')
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill="forestgreen", color="black") +
  coord_quickmap()

#add radius information related to number of sharks sampled and play with the colors and opacity
blood_haps$radius <- .0017* blood_haps$total_bt_sharks

p+ geom_scatterpie(aes(x=Long, y=Lat, group=Site, r=radius), data=blood_haps,
                   cols=c("percent_hap1", "percent_hap2"), alpha=.7) + 
  #coord_equal()+
  scale_fill_manual(values = c("orange","violet"),
                    labels=c("haplotype 1","haplotype 2"))+
  coord_sf(xlim = c(-89, -90), ylim = c(-.6, -1))


#Download Ecuador vector map

ecuador_adm <- geodata::gadm(country = "Ecuador", level = 2, path = "maps/")
#takes a really long time to plot ecuador_adm
#plot(ecuador_adm)
#view columns and values
ecuador_adm
#view values for NAME_1
ecuador_adm$NAME_1 #these look like State-level groups
#view values for NAME_2
ecuador_adm$NAME_2 #These look like county-level groups

#Subset the Ecuador vector map to just the Galápagos and just San Cristóbal.

San_cristobal_adm <- subset(ecuador_adm, ecuador_adm$NAME_2 =="San Cristóbal")
plot (San_cristobal_adm)

Galapagos_adm<-subset(ecuador_adm, ecuador_adm$NAME_1 =="Galápagos")
plot(Galapagos_adm)

#Load bathymetry raster for the whole Galapagos region downloaded from GEBCO: https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2019/grid_terms_of_use.html

bathy <- terra::rast("../../Physalia_GIS_in_R_2023/course_materials/practicals/data/bathymetry_galapagos/GEBCO_30_Oct_2023_672896985ee1/gebco_2023_n2.3839_s-2.9182_w-93.4343_e-85.3989.tif")
plot(bathy, col = hcl.colors(100, "blues"))
#add Galapagos vector on top
plot(Galapagos_adm, add=TRUE, col="forestgreen")

#Load bathymetry raster for the San Cristóbal region downloaded from GEBCO: https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2019/grid_terms_of_use.html

sancristobal_bathy<- terra::rast("../Maps/GEBCO_30_Oct_2023_519f2751c841/gebco_2023_n-0.5651_s-1.0455_w-89.721_e-89.1533.tif")
plot(sancristobal_bathy, col = hcl.colors(100, "Blues"))

plot(San_cristobal_adm, add=TRUE, col= "forestgreen")

##### Plot the Island with the Scatterpies (still can't add bathymetry raster) #####

#read in depth raster
sancristobal_bathy<- terra::rast("../Maps/GEBCO_30_Oct_2023_519f2751c841/gebco_2023_n-0.5651_s-1.0455_w-89.721_e-89.1533.tif")
#change names of bathymetry raster values
names(sancristobal_bathy)[1]<-"Depth"

#read in Ecuador vector
ecuador_adm <- geodata::gadm(country = "Ecuador", level = 2, path = "maps/")

#crop to San Cristóbal and surrounding minor islands
San_cristobal_adm <- subset(ecuador_adm, ecuador_adm$NAME_2 =="San Cristóbal")

#find extent of bathymetry
bathy_ext<-terra::ext(sancristobal_bathy)

#crop administrative vector to bathymetry extent
San_cristobal_island<-terra::crop(x=San_cristobal_adm,
                                  y=bathy_ext)

#set up basemap with bathymetry
base_map<-ggplot(data=sancristobal_bathy)+
  geom_raster(mapping=aes(x=x,y=y,fill=Depth))
base_map

#add the island vector
island<-base_map+
  geom_spatvector(data=San_cristobal_island,color="tan",fill="forestgreen")
island

#make the pies
blood_haps<- read.csv("metadata/percentage_of_haps_in_bloods.csv")
blood_haps<-blood_haps%>% arrange(total_bt_sharks)
blood_haps$radius <- .0017* blood_haps$total_bt_sharks
blood_haps_SCzoom<-blood_haps[-c(2, 3), ]

#add the scatterpies

#The code below creates the following error:Error: Discrete value supplied to continuous scale
#island+ geom_scatterpie(data=blood_haps,
#                       mapping = aes(x=Long, y=Lat, group=Site, r=radius),
#                       cols=c("percent_hap1", "percent_hap2"))

#So, try without bathymetry raster
island_vect<-ggplot()+
  geom_spatvector(data=San_cristobal_island,fill="forestgreen", color="tan")
island_vect

#add pies

pie_map<-island_vect+ geom_scatterpie(aes(x=Long, y=Lat, group=Site, r=radius), data=blood_haps_SCzoom,
                             cols=c("percent_hap1", "percent_hap2"), alpha=.8, color="NA") + 
  geom_label_repel(aes(x=blood_haps_SCzoom$Long, y=blood_haps_SCzoom$Lat, 
                       label = blood_haps_SCzoom$Site), 
                   size = 3, nudge_x = -.04, nudge_y = 0.06,point.padding=2.8,
                   segment.colour = "tan")+
  scale_fill_manual(values = c("red","pink"),
                    labels=c("haplotype 1","haplotype 2"))+
  geom_shadowpoint()+
  xlab("Longitude")+
  ylab("Latitude")+
  #ggtitle("Short D-loop Haplotypes of C.limbatus from Catch-and-Release Sampling")+
  theme(panel.background=element_rect(fill = "darkblue"),panel.grid = element_blank())

pie_map

ggsave("viz/SanCristóbal2021_Climb2_bloodhaps.jpg")

#Next steps:
#add a Zoomed-in box for Rosa Blanca 1 and 2 separately in the bottom right of the figure.


#install.packages("ggmagnify", repos = c("https://hughjonesd.r-universe.dev", 
#                                       "https://cloud.r-project.org"))
library(ggmagnify)

#pie_map + geom_magnify(aes(from = blood_haps$Site == c("Rosa Blanca 1", "Rosa Blanca 2")), 
#                    to = c(89.35, 0.96, 0.96, 0.89))

#this isn't working yet, so try cowplot
library(cowplot)

##### Make a zoomed-in map of Rosa Blanca 1 and 2 ####

RosaBlanca_blood_haps<-blood_haps[-c(1,4,5,6), ]
RosaBlanca_blood_haps$radius <- .0006* RosaBlanca_blood_haps$total_bt_sharks


island_vect+
  coord_sf(xlim = c(-89.38, -89.3), ylim = c(-.85, -.8))

RosaBlanca_pie_map<-island_vect+ 
  coord_sf(xlim = c(-89.38, -89.32), ylim = c(-.85, -.8))+
  geom_scatterpie(aes(x=Long, y=Lat, group=Site, r=radius), data=RosaBlanca_blood_haps,
                  cols=c("percent_hap1", "percent_hap2"), alpha=.8, color="NA") + 
  geom_label_repel(aes(x=RosaBlanca_blood_haps$Long, y=RosaBlanca_blood_haps$Lat, 
                       label = RosaBlanca_blood_haps$Site), 
                   size = 8,point.padding=2.8,nudge_x = .01, nudge_y = -0.015,
                   segment.colour = "tan")+
  scale_fill_manual(values = c("red","pink"),
                    labels=c("haplotype 1","haplotype 2"))+
  geom_shadowpoint()+
  xlab("Longitude")+
  ylab("Latitude")+
  #ggtitle("Short D-loop Haplotypes of C.limbatus from Catch-and-Release Sampling")+
  theme(panel.background=element_rect(fill = "darkblue"),panel.grid = element_blank())

RosaBlanca_pie_map
ggsave2("viz/RosaBlanca_pie_map.png")
RosaBlanca_insert<-readPNG("viz/RosaBlanca_pie_map.png")

ggdraw(pie_map) + 
  draw_image(RosaBlanca_insert, x = 1, y = 0.4, hjust = 1, vjust = 1, width = 0.3, height = 0.3)

#### DONE!  Just need to change the size of the output to the zoomed level and ideally add the bathymetry map instead of just the darkblue background.

##### End of working section #####

##### Map the eDNA results ######









#trying to add the bathymetry as a png 
#not working because I would need to adjust the lat and long extent of the scatterpies to approximate where it should put them on the map that is now just a picture.  

#following Jeff's suggestions about specifying aes for each layer explicitly

#set up basemap with bathymetry
bathy_map<-ggplot(data=sancristobal_bathy)+
  geom_raster(mapping=aes(x=x,y=y,fill=Depth))

island_map<-bathy_map+
  geom_spatvector(data=San_cristobal_island,fill="forestgreen", color="tan")

island_map

pie_map2<-island_map +
  geom_scatterpie(aes(x=blood_haps_SCzoom$Long, y=blood_haps_SCzoom$Lat, group=Site, r=radius), data=blood_haps_SCzoom,
                  cols=c("percent_hap1", "percent_hap2"), alpha=1, color="NA") + 
  #geom_label_repel(aes(x=blood_haps_SCzoom$Long, y=blood_haps_SCzoom$Lat, 
  #                     label = blood_haps_SCzoom$Site), 
  #                 size = 3, nudge_x = -.04, nudge_y = 0.06,point.padding=2.8,
  #                 segment.colour = "tan")+
  scale_fill_manual(values = c("red","pink"),
                    labels=c("haplotype 1","haplotype 2"))+
  #geom_shadowpoint()+
  xlab("Longitude")+
  ylab("Latitude")+
  #ggtitle("Short D-loop Haplotypes of C.limbatus from Catch-and-Release Sampling")+
  theme(panel.background=element_rect(fill = "darkblue"),panel.grid = element_blank())

pie_map2
#still not working




#export the basemap with desired extend as a png
library(png)
library(grid)
base_map
ggsave("viz/bathypic.png")
#import the bathymetry as a png
basepic<-readPNG("viz/bathypic.png")

#add the png as a background image

island_vect+ geom_scatterpie(aes(x=Long, y=Lat, group=Site, r=radius), data=blood_haps,
                             cols=c("percent_hap1", "percent_hap2"), alpha=1) + 
  scale_fill_manual(values = c("orange","violet"),
                    labels=c("haplotype 1","haplotype 2"))+
  coord_sf(xlim = c(-89, -90), ylim = c(-.6, -1))+
  annotation_custom(rasterGrob(image = basepic,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")))

#add the basepic first, then overlay the scatterpies on top of it... island not working now
ggplot()+
  annotation_custom(rasterGrob(image = basepic,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")))+
  geom_scatterpie(aes(x=Long, y=Lat, group=Site, r=radius), data=blood_haps,
                  cols=c("percent_hap1", "percent_hap2"), alpha=1) + 
  scale_fill_manual(values = c("orange","violet"),
                    labels=c("haplotype 1","haplotype 2"))+
  coord_sf(xlim = c(-89, -90), ylim = c(-.6, -1))

  

island_vect+ geom_scatterpie(aes(x=Long, y=Lat, group=Site, r=radius), data=blood_haps,
                             cols=c("percent_hap1", "percent_hap2"), alpha=1) + 
  scale_fill_manual(values = c("orange","violet"),
                    labels=c("haplotype 1","haplotype 2"))+
  coord_sf(xlim = c(-89, -90), ylim = c(-.6, -1))

