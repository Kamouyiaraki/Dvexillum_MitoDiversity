### Chapter 4 Analysis
setwd("../Desktop/Ch4")

library(ggplot2)
library(ggrepel)
library(measurements)
library(rgdal)
library(pegas)
library(haplotypes)
library(ggmap)
library(phylotools)
library(tidyr)
library(purrr)
library(dplyr)
library(ggthemes)


### table of location, location name, frequencies of haplotype1, frequencies of haplotype2
#### for COI
#### for NAD1
#### for combined

### df of coordinates

##### Estimate centroid of coordinates
#df<-read.csv("../maria scottish Dvex_cleaned.csv", header = T)
df<-read.csv("./Dvex_coords_cleaned.csv", header = T)
samples<-read.csv("./samples.csv", header = T)
samp<-samples$Sample.ID
sampcheck<-df$sampleID

##Not run
#df$coords<-paste(df$lat_Northing, "N ", "-", df$long_West, "E", sep="")
#df$coords<- gsub("NAN -NAE", "NA", df$coords)
#write.csv(df, "../Dvex_coords_cleaned.csv", row.names = F)

###to double check that this contains all the correct samples
wsamp<-which(!is.element(samp,sampcheck))
#samples$Sample.ID[wsamp]


## Convert deg with decimal minutes to decimal degrees
df$lat2<- measurements::conv_unit(df$lat_Northing, from = "deg_dec_min", to ="dec_deg")
df$long2<- measurements::conv_unit(df$long_West, from = "deg_dec_min", to ="dec_deg")
df$lat2<-as.numeric(df$lat2)
df$long2<-as.numeric(df$long2)
df$long2<-df$long2*(-1)

lgdf<-df[df$loc == "Largs ",]
lgdf<-lgdf[complete.cases(lgdf),]

lcdf<-df[df$loc== "Creran",]
lcdf<-lcdf[complete.cases(lcdf),]

fofdf<-df[df$loc == "Clyde area",]
fofdf<-fofdf[complete.cases(fofdf),]

pvdf<-df[df$loc == "Portavadie",]
pvdf<-pvdf[complete.cases(pvdf),]

####both functions taken from https://livefreeordichotomize.com/2018/06/27/bringing-the-family-together-finding-the-center-of-geographic-points-in-r/

geographic_average <- function(lon, lat, weight = NULL) {
  if (is.null(weight)) {
    weight <- rep(1, length(lon))
  }
  lon <- weighted.mean(lon, w = weight)
  lat <- weighted.mean(lat, w = weight)
  data.frame(lon = lon, lat = lat)
}

geographic_midpoint <- function(lon, lat, weight = NULL) {
  if (is.null(weight)) {
    weight <- rep(1, length(lon))
  }
  # degrees to radians
  lat <- lat * pi / 180
  lon <- lon * pi / 180
  # cartesian coordinates
  x <- cos(lat) * cos(lon)
  y <- cos(lat) * sin(lon)
  z <- sin(lat)
  # weighted mean
  x <- weighted.mean(x, w = weight)
  y <- weighted.mean(y, w = weight)
  z <- weighted.mean(z, w = weight)
  # convert to lat and lon
  lon <- atan2(y, x) * 180 / pi
  hyp <- sqrt(x * x + y * y)
  lat <- atan2(z, hyp) * 180 / pi
  
  data.frame(lon = lon, lat = lat)
}

dfmid<-data.frame(site = c("Fairlie", "Largs", "Loch Creran", "Portavadie"), lat = rep(NA, 4), lon = rep(NA,4))
dfmid[1,2:3]<-geographic_midpoint(fofdf$lat2, fofdf$long2, weight=NULL)
dfmid[2,2:3]<-geographic_midpoint(lgdf$lat2, lgdf$long2, weight=NULL)
dfmid[3,2:3]<-geographic_midpoint(lcdf$lat2, lcdf$long2, weight=NULL)
dfmid[4,2:3]<-geographic_midpoint(pvdf$lat2, pvdf$long2, weight=NULL)

### Add in haplotype frequency data
df<- dfmid


fas<- "NAD1_aln.fas"
x<-read.fas(fas)
h<-haplotypes::haplotype(x,indels = "missing")
seqnames<-names(x)
seqnames[startsWith(seqnames,"Lc")]<- "Loch Creran"
seqnames[startsWith(seqnames,"Lg")]<- "Largs"
seqnames[startsWith(seqnames,"PvT")]<- "Portavedie"
seqnames[startsWith(seqnames,"FOF")] <- "Fairlie"
g<-grouping(h,factors=seqnames) 

df$hapA<-as.numeric(g$hapmat[1,])
df$hapE<-as.numeric(g$hapmat[2,])



fas<- "coi_haplotypes_aln.fas"
x<-read.fas(fas)
h<-haplotypes::haplotype(x,indels = "missing")
seqnames<-names(x)
seqnames[startsWith(seqnames,"Lc")]<- "Loch Creran"
seqnames[startsWith(seqnames,"Lg")]<- "Largs"
seqnames[startsWith(seqnames,"PvT")]<- "Portavedie"
seqnames[startsWith(seqnames,"FOF")] <- "Fairlie"
g<-grouping(h,factors=seqnames) 

df$hap3<-as.numeric(g$hapmat[1,])
df$hap2<-as.numeric(g$hapmat[2,])




df$N<-(df$hapA + df$hapE)
df$radius<-df$N
df$radius<-(log10(df$N)/10)

df$sitetype <- c("Oyster Farm", "Marina", "Oyster Farm", "Marina")


### Note: Fairlie and Largs overlap 
### manually shift up Largs and shift down Fairlie (lat)
df$lat[df$site == "Fairlie"]<- 55.65
df$lat[df$site == "Largs"]<- 55.8


head(df)
#         site      lat      long hapA hapE hap3 hap2  N    radius    sitetype
#1     Fairlie 55.65000 -4.870946    4    0    4    0  4 0.0602060 Oyster Farm
#2       Largs 55.80000 -4.857967    6    2    6    2  8 0.0903090      Marina
#3 Loch Creran 56.51321 -5.387106   36    0   36    0 36 0.1556303 Oyster Farm
#4  Portavadie 55.87270 -5.303017   10    3   10    3 13 0.1113943      Marina


#### Map of locations and Sample size (final)

#if(!requireNamespace("devtools")) install.packages("devtools")
#devtools::install_github("dkahle/ggmap", ref = "tidyup") #from the tidyup branch on Github

library(ggmap)
register_google(key= "xxxxx")

centre<- c(lon = -5, lat = 56)
map<-ggmap(get_googlemap(center = centre, matype = "terrain", zoom = 8, style = 'element:labels|visibility:off')) #terrain without labels


legend_size <- c(5, 8, 12, 16)

map +
  geom_label_repel(data = df, aes(x = lon, y = lat,  label = site),family = "serif", size = 5, inherit.aes = FALSE, size = 3.5, fontface = "bold", nudge_x = c(1, 1, 1, 1), nudge_y = c(0.25,0.25, 0.25, 0.25))+
  geom_point(data = df, aes(size = N, color = sitetype))+
  scale_color_manual(values=c("#a6611a","#018571"))+
  scale_size(breaks = sort(df$N), range = c(5, 16), name="N")+
  guides(
    size=guide_legend(override.aes = list(size = legend_size)))+
  labs(y= "Latitude", x = "Longitude")+
  theme_tufte()+
  theme(text = element_text(size=20), axis.line= element_line(size = 0.5, color = "black"))



#### haplotype map (final)

#pie_col<-c("#5ab4ac", "#d8b365")
pie_col<-c("#264653", "#e76f51")
scales::show_col(pie_col)

pie<- df
pie<-pie[c("lat", "lon", "hapA", "hapE", "radius")]
colnames(pie)<-c("lat", "lon", "A", "E", "radius")

#### code from: https://stackoverflow.com/questions/51398344/r-pie-charts-distorted-when-adding-to-projected-map-using-ggplot

pie.list <- pie %>% 
  tidyr::gather(type, value, -lon, -lat, -radius) %>% #turns into long format 
  tidyr::nest(type, value) %>%##turns it into a tibble of lat, lon, radius and data (data = nested list of type (haplotype) and value (N))
  
  # make a pie chart from each row, & convert to grob
  dplyr::mutate(pie.grob = purrr::map(data,
                                      function(d) ggplotGrob(ggplot(d, 
                                                                    aes(x = 1, y = value, fill = type)) +
                                                               geom_col(color = NA,
                                                                        show.legend = FALSE) +
                                                               coord_polar(theta = "y") +
                                                               scale_fill_manual(values=pie_col)+
                                                               theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = lon - radius, xmax = lon + radius,
                                          ymin = lat - radius, ymax = lat + radius)))


m<-map+
  geom_label_repel(data = df, aes(x = lon, y = lat,  label = site),inherit.aes = FALSE, family = "serif", size = 5, fontface = "bold", nudge_x = c(1, 1, 1, 1), nudge_y = c(0.25,0.25, 0.25, 0.25))+
  coord_quickmap(xlim = c(-6.7, -3.25), ylim = c(55, 56.9))+
  theme_tufte()+
  theme(text = element_text(size=20), axis.line= element_line(size = 0.5, color = "black"))


m+
  # Optional. this hides some tiles of the corresponding color scale BEHIND the
  # pie charts, in order to create a legend for them
  geom_tile(data = pie %>% tidyr::gather(type, value, -lon, -lat, -radius),
            aes(x = lon, y = lat, fill = type)) +
  pie.list$subgrob+
  scale_fill_manual(values=pie_col)+
  labs(fill = "Haplotype ID", y= "Latitude", x = "Longitude")


#### Haplotype Network
source("hap.custom.R")

phylotools::rm.sequence.fasta("NAD1_aln_to_ch3.fas", "NAD1_aln_to_ch3_v2.fas", to.rm = c("ch3_A19_mito", "ch3_H2_mito")) #removed haplotypes A and E

#nad1_hap<-pegas::haplotype(seq, strict = T, labels = c("B", "C", "D", "F", "G", "A","E"))
nad1_hap<- hap.custom("NAD1_aln_to_ch3_v2.fas", labels = c("B", "C", "D", "F", "G", "A","E"))


labels(nad1_hap)
nad1_net<-pegas::haploNet(nad1_hap)

fasta<-read.dna("NAD1_aln_to_ch3_v2.fas" , format="fasta")
rownames(fasta)<-gsub("Lc","LochCreran_Lc", rownames(fasta))
rownames(fasta)<-gsub("Lg","Largs_Lg", rownames(fasta))
rownames(fasta)<-gsub("PvT","Portavadie_Pvt", rownames(fasta))
rownames(fasta)<-gsub("FOF","Fairlie_FOF", rownames(fasta))

ind.hap<-with(utils::stack(setNames(attr(nad1_hap, "index"), rownames(nad1_hap))), table(hap=ind, pop=rownames(fasta)[values]))
#rownames(ind.hap)<- c("B", "A", "C", "D", "E", "F", "G")
mydata<-as.data.frame(ind.hap)
loc<-mydata[mydata$Freq ==1,]
locations<-strsplit(as.character(loc$pop), "_")
locations<-sapply(locations, "[[",1)
#locations[locations=="clade"]<-"ref"
newhap_nad1<-table(loc$hap, locations)

colval<-c("ch3" = "white", "Fairlie" = "#dfc27d", "Largs"= "#018571", "LochCreran" = "#a6611a", "Portavadie"="#80cdc1")
scales::show_col(colval)
colval_legend<-c("Fairlie" = "#dfc27d", "LochCreran" = "#a6611a", "Largs"= "#018571", "Portavadie"="#80cdc1")
legend_names<-c("Fairlie" , "Loch Creran", "Largs", "Portavadie")

#plot(nad1_net, scale.ratio = 0.5, cex = 0.8, pie=newhap_nad1, bg=colval, show.mutation=2)  #option for legend = T
par(family = "serif")
plot(nad1_net, size=log(attr(nad1_net, "freq")), scale.ratio = 1, cex = 1, pie=newhap_nad1, bg=colval, show.mutation=2)  #option for legend = T
replot()
legend("right", legend_names, col=colval_legend, pch=19, cex = 1.3) ##not all to remove chapter 3 in the key

coi_hap<-hap.custom(fastafile = "coi_haplotypes_aln_allhaps.fas", labels = c(3,2,1,4,5,6))

labels(coi_hap)
coi_net<-pegas::haploNet(coi_hap)

fasta<-read.dna("coi_haplotypes_aln_allhaps.fas" , format="fasta")
rownames(fasta)<-gsub("Lc","LochCreran_Lc", rownames(fasta))
rownames(fasta)<-gsub("Lg","Largs_Lg", rownames(fasta))
rownames(fasta)<-gsub("PvT","Portavadie_Pvt", rownames(fasta))
rownames(fasta)<-gsub("FOF","Fairlie_FOF", rownames(fasta))

ind.hap<-with(utils::stack(setNames(attr(coi_hap, "index"), rownames(coi_hap))), table(hap=ind, pop=rownames(fasta)[values]))
#rownames(ind.hap)<- c("B", "A", "C", "D", "E", "F", "G")
mydata<-as.data.frame(ind.hap)
loc<-mydata[mydata$Freq ==1,]
locations<-strsplit(as.character(loc$pop), "_")
locations<-sapply(locations, "[[",1)
#locations[locations=="clade"]<-"ref"
newhap_coi<-table(loc$hap, locations)

colval<-c("Fairlie" = "#dfc27d", "haplotype" = "white", "Largs"= "#018571", "LochCreran" = "#a6611a", "Portavadie"="#80cdc1")
scales::show_col(colval)
colval_legend<-c("Fairlie" = "#dfc27d", "LochCreran" = "#a6611a", "Largs"= "#018571", "Portavadie"="#80cdc1")
legend_names<-c("Fairlie" , "Loch Creran", "Largs", "Portavadie")

#plot(coi_net, scale.ratio = 0.5, cex = 0.8, pie=newhap_coi, bg=colval, show.mutation=2)  #option for legend = T
plot(coi_net, size=log(attr(coi_net, "freq")), scale.ratio = 1, cex = 1, pie=newhap_coi, bg=colval, show.mutation=2)  #option for legend = T
replot()
legend("right", legend_names, col=colval_legend, pch=19, cex = 1) ##not all to remove chapter 3 in the key

#### Haplotype diversity

source("hap.custom.R")


x<-read.fas("NAD1_aln.fas")
seqnames<-names(x)
seqnames[startsWith(seqnames,"Lc")]<- "Loch Creran"
seqnames[startsWith(seqnames,"Lg")]<- "Largs"
seqnames[startsWith(seqnames,"PvT")]<- "Portavedie"
seqnames[startsWith(seqnames,"FOF")] <- "Fairlie"


nad1cust<-hap.custom(fastafile = "NAD1_aln.fas", indel="missing", hapdiv = TRUE, pops = seqnames, variance = TRUE)

#Pop   Hap.Div   Variance
#1     Fairlie 0.0000000 0.00000000
#2       Largs 0.4285714 0.02762277
#3 Loch Creran 0.0000000 0.00000000
#4  Portavedie 0.3846154 0.01730439

x2<-read.fas("coi_haplotypes_aln.fas")
seqnames<-names(x)
seqnames[startsWith(seqnames,"Lc")]<- "Loch Creran"
seqnames[startsWith(seqnames,"Lg")]<- "Largs"
seqnames[startsWith(seqnames,"PvT")]<- "Portavedie"
seqnames[startsWith(seqnames,"FOF")] <- "Fairlie"


coicust<-hap.custom(fastafile = "coi_haplotypes_aln.fas", indel="missing", hapdiv = TRUE, pops = seqnames, variance = TRUE)
#Pop   Hap.Div   Variance
#1     Fairlie 0.0000000 0.00000000
#2       Largs 0.4285714 0.02762277
#3 Loch Creran 0.0000000 0.00000000
#4  Portavedie 0.3846154 0.01730439



#### Punnett square

fasta<-read.dna("coi_haplotypes_aln_allhaps.fas" , format="fasta")
rownames(fasta)<-gsub("Lc","LochCreran_Lc", rownames(fasta))
rownames(fasta)<-gsub("Lg","Largs_Lg", rownames(fasta))
rownames(fasta)<-gsub("PvT","Portavadie_Pvt", rownames(fasta))
rownames(fasta)<-gsub("FOF","Fairlie_FOF", rownames(fasta))

coi_hap<-hap.custom(fastafile = "coi_haplotypes_aln_allhaps.fas", labels = c(3,2,1,4,5,6))

labels(coi_hap)
coi_net<-pegas::haploNet(coi_hap)

ind.hap<-with(utils::stack(setNames(attr(coi_hap, "index"), rownames(coi_hap))), table(hap=ind, pop=rownames(fasta)[values]))
#rownames(ind.hap)<- c("B", "A", "C", "D", "E", "F", "G")
mydata<-as.data.frame(ind.hap)
loc<-mydata[mydata$Freq ==1,]


sq<-data.frame(Ncoi3 =  c(sum(loc[loc$hap ==3, "Freq"]), 0), Ncoi2 = c(0, sum(loc[loc$hap ==2, "Freq"]) ))
sq$nad<- c("hapA", "hapE")

sq_long<-reshape2::melt(sq, id.vars=c("nad"))
sq_long$variable<- gsub("Ncoi2", "2", sq_long$variable)
sq_long$variable<- gsub("Ncoi3", "3", sq_long$variable)
sq_long$nad<- gsub("hapA", "A", sq_long$nad)
sq_long$nad<- gsub("hapE", "E", sq_long$nad)

library(ggplot2)
library(ggthemes)
  
ggplot(sq_long, aes(nad, variable)) + 
  geom_tile(colour='black', fill='white')+
  geom_text(aes(label = value),size = 8, family="serif") +
  ylab("Partial-COI Haplotype ID")+
  xlab("NAD1 Haplotype ID") +
  theme_tufte()+
  theme(text = element_text(size=20))


