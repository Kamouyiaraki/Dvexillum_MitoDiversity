#### Alingment formats - no deletion 2020-06-11

Mito<-seqinr::read.fasta("Mito_aln.fas")
Win30 <-seqinr::read.fasta("Rooted_Win_30.fas")
Amp1 <-seqinr::read.fasta("Rooted_Dv1F2.fas")
Amp2 <-seqinr::read.fasta("Rooted_Dv1F3.fas")
barcode<- seqinr::read.fasta("Rooted_TUN_COI.fas")
Reg1<-seqinr::read.fasta("Rooted_Region1.fas")


### create PHYLIP format (PhyML, JModelTest, RAXML and PartitionFinder requirement)
ape::write.dna(Mito, file="./PHYLIP/Rooted_Mito.phy", format="interleaved",  nbcol=1, colsep="", colw=1000000)
ape::write.dna(Win30, file="./PHYLIP/Rooted_Win30.phy", format="interleaved",  nbcol=1,colsep="", colw=1000000)
ape::write.dna(Amp1, file="./PHYLIP/Rooted_Amp1.phy", format="interleaved",  nbcol=1,colsep="", colw=1000000)
ape::write.dna(Amp2, file="./PHYLIP/Rooted_Amp2.phy", format="interleaved",  nbcol=1,colsep="", colw=1000000)
ape::write.dna(Reg1, file="./PHYLIP/Rooted_Reg1.phy", format="interleaved",  nbcol=1,colsep="", colw=1000000)
ape::write.dna(barcode, file="./PHYLIP/Rooted_TUN_COI.phy", format="interleaved",  nbcol=1, colsep="", colw=1000000)



### Create NEXUS format (BEAuTi and BEAST requirement)
ape::write.nexus.data(Mito, file="./NEXUS/Rooted_Mito.nexus", format="dna")
ape::write.nexus.data(Win30, file="./NEXUS/Rooted_Win30.nexus", format="dna")
ape::write.nexus.data(Amp1, file="./NEXUS/Rooted_Amp1.nexus", format="dna")
ape::write.nexus.data(Amp2, file="./NEXUS/Rooted_Amp2.nexus", format="dna")
ape::write.nexus.data(Reg1, file="./NEXUS/Rooted_Reg1.nexus", format="dna")
ape::write.nexus.data(barcode, file="./NEXUS/Rooted_TUN_COI.nexus", format="dna")



#### Library

library(seqinr)
library(ape)
library(adegenet)
library(pegas)
library(Biostrings)
library(plotrix)
library(reshape2)
library(geosphere)
library(scales)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(phytools)
library(ggrepel)
library(stringr)
library(stringi)
library(abind)
library(treeio)
library(ggplot2)


### Colour palette & legend order
altcol<-RColorBrewer::brewer.pal(10, "Paired")
altcol[9]<- "grey"

names(altcol)<- c("England", "Wales", "Scotland", "N.Ireland", "Alaska", "Japan", "Germany", "France", "Canada", "New.Zealand")
scales::show_col(altcol)

colval<-altcol[order(factor(names(altcol)))]
#colval<-c(colval, "black")
#names(colval)[11]<- "Clade B"


##colourblind palette
colval<-c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")
names(colval)<- c("England", "Wales", "Scotland", "N.Ireland",  "Germany", "France", "Japan","Alaska", "New.Zealand", "Canada")
scales::show_col(colval)
colval<-colval[order(factor(names(colval)))]

legend_ord<-c("Alaska", "Canada","England", "Scotland","Wales", "N.Ireland", "Germany", "France",  "New.Zealand","Japan" )
#legend_ord<-c("Alaska", "Canada", "Japan", "New.Zealand", "N.Ireland", "England", "Scotland","Wales", "Germany", "France", "Clade B")

alllabels<-c("Scotland_58A","Scotland_FOF1b", "Scotland_66A", "France_FBM05", "France_FBM06" , "France_FBM09", "England_D1", "England_G1", "England_G2" , "England_KB", "England_BK01", "Alaska_A02",  "Alaska_A09", "Alaska_A19", "Canada_CAN49", "Canada_CAN50", "Canada_CAN52", "New Zealand_NZ034","New Zealand_NZ108","New Zealand_NZ224", "Japan_J296","Japan_J390", "Japan_J431", "Japan_J424", "Germany_GS01", "Germany_GS02", "Germany_GS03", "N. Ireland_Ire2", "Wales_H2")
#alllabels<-c("Scotland_58A","Scotland_FOF1b", "Scotland_66A", "France_FBM05", "France_FBM06" , "France_FBM09", "England_D1", "England_G1", "England_G2" , "England_KB", "England_BK01", "Alaska_A02",  "Alaska_A09", "Alaska_A19", "Canada_CAN49", "Canada_CAN50", "Canada_CAN52", "New Zealand_NZ034","New Zealand_NZ108","New Zealand_NZ224", "Japan_J296","Japan_J390", "Japan_J431", "Japan_J424", "Germany_GS01", "Germany_GS02", "Germany_GS03", "N. Ireland_Ire2", "Wales_H2", "Japan_J387")


### Import phylogenetic tree
coi_tree<-ape::read.tree("Rooted_TUN_COI.contree")     ##tree specific

coi_tree <- ape::drop.tip(coi_tree, "J387_mito")
#coi_tree<-root(coi_tree, outgroup = "J387_mito")

### Tree data
coi_tree$tip.label<-gsub("58A_mito","Scotland_58A (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("FOF1b_mito","Scotland_FOF1b (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("66A_mito","Scotland_66A (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("FBM05_mito","France_FBM05 (2)", coi_tree$tip.label)  ####
coi_tree$tip.label<-gsub("FBM06_mito","France_FBM06 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("FBM09_mito","France_FBM09 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("D1_mito","England_D1 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("G1_mito","England_G1 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("G2_mito","England_G2 (5)", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("KB_mito","England_KB (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("BK01_mito","England_BK01 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("A02_mito","Alaska_A02 (6)", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("A09_mito","Alaska_A09 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("A19_mito","Alaska_A19 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("CAN49_mito","Canada_CAN49 (5)", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("CAN50_mito","Canada_CAN50 (5)", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("CAN52_mito","Canada_CAN52 (5)", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("NZ034_mito","New.Zealand_NZ034 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("NZ108_mito","New.Zealand_NZ108 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("NZ224_mito","New.Zealand_NZ224 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("J296_mito","Japan_J296 (3)", coi_tree$tip.label)
#coi_tree$tip.label<-gsub("J387_mito","CladeB_Japan_J387", coi_tree$tip.label)
coi_tree$tip.label<-gsub("J390_mito","Japan_J390 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("J424_mito","Japan_J424 (2)", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("J431_mito","Japan_J431 (5)", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("GS01_mito","Germany_GS01 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("GS02_mito","Germany_GS02 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("GS03_mito","Germany_GS03 (3)", coi_tree$tip.label)
coi_tree$tip.label<-gsub("Ire2_mito","N.Ireland_Ire2 (1)", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("H2_mito","Wales_H2 (2)", coi_tree$tip.label) ####

##get tree data to manipulate bootstrap values (remove <50%) and add posterior probabilities
q1 <- ggtree(coi_tree)                                   #### tree specfic
a <- q1$data    # extract bootstraps in node metadata
#d<-q1$data

# Specify group based on location
a$group<-a$label
a$group[a$isTip == FALSE] <- NA
a$group<-strsplit(as.character(a$group), "_")
a$group<-sapply(a$group, "[[",1)

# Add colours to use for each
a$color<-a$group
a$color<-gsub("Alaska","#FB9A99",a$color)
a$color<-gsub("Canada", "grey" , a$color)
a$color<-gsub("Japan", "#E31A1C", a$color)
a$color<-gsub("New.Zealand", "#6A3D9A", a$color)
a$color<-gsub("N.Ireland", "#33A02C", a$color)
a$color<-gsub("England", "#A6CEE3", a$color)
a$color<-gsub("Wales", "#1F78B4", a$color)
a$color<-gsub("Scotland", "#B2DF8A", a$color)
a$color<-gsub("Germany", "#FDBF6F", a$color)
a$color<-gsub("France", "#FF7F00", a$color)
#a$color<-gsub("CladeB", "black", a$color)


# Remove low bootstraps and tip labels
a$label[a$isTip == TRUE] <- NA  #keep only values for non-tip nodes
a$label[a$label == ""]<- NA   #change root node to a value (by default empty)
#a$label[a$label == "Root"]<- NA  
a$label<-as.numeric(a$label)   #change to numeric
a$label[a$label < 65] <- NA  #remove all bootstraps of less than 50%

# Add posterior probabilities                              #### tree specific
a$label[a$label ==91]<-"91/1.000"
a$label[a$label ==66]<- "66/0.734"
a$label[a$label ==62]<-NA
a$label[a$label ==75]<-"75/1.000"
a$label[a$label ==61]<-NA
a$label[a$label ==77]<-NA
a$label[a$label ==81]<-NA
a$label[a$label ==84]<-NA
a$label[a$label ==87]<-NA
a$label[a$label ==74]<-NA
a$label[a$label ==65]<-NA
a$label[a$label ==72]<-NA
a$label[a$label ==59]<-NA
a$label[a$label ==64]<-NA
a$label[a$node == 50]<- "55/0.724"
a$label[a$node == 42]<-NA

# Format node labels
#a$label<-ifelse(!grepl("/", a$label), paste0(a$label, "/-"), paste(a$label))
#a$label<-gsub("0/-", "NA", a$label)
#a$label<-gsub("0", NA, a$label)
#a<-a %>% replace_with_na(replace = list(label = "NA"))

# Create a df of metadata for each location for colours
meta_data1<-as.data.frame(a[,c("group", "label")])
meta_data1$group<-as.character(meta_data1$group)
meta_data1<-meta_data1[1:length(alllabels),]
rownames(meta_data1)<-coi_tree$tip.label                #### tree specfic
meta_data1$label<-NA

# Plot tree with node support
p1 <- ggtree(coi_tree) +                                   #### tree specfic
  geom_treescale()+
  geom_rootedge(rootedge = 0.001)+
  geom_tiplab(size = 4, align=FALSE, family = "serif") +
  geom_text2(size = 4, aes(label=a$label), color="black", family = "serif", nudge_x = -0.0008, nudge_y = 0.2)



# Add colours for each location
hm1<-gheatmap(p1,meta_data1, offset = 0.002, width=0.1, font.size=4, colnames=FALSE, hjust = 0)+ scale_fill_manual(breaks = legend_ord, values=colval) +
  theme(legend.position = "left", text = element_text(family = "serif", size = 15), legend.title = element_blank(), legend.key = element_blank()) # no keys
plot(hm1)


###############################################################################################################################################

reg_tree<-ape::read.tree("../rooted_reg2.contree")

reg_tree <- ape::drop.tip(reg_tree, "J387_mito")
#reg_tree<-root(reg_tree, outgroup = "J387_mito")

### Tree data
reg_tree$tip.label<-gsub("58A_mito","Scotland_58A", reg_tree$tip.label)
reg_tree$tip.label<-gsub("FOF1b_mito","Scotland_FOF1b", reg_tree$tip.label)
reg_tree$tip.label<-gsub("66A_mito","Scotland_66A", reg_tree$tip.label)
reg_tree$tip.label<-gsub("FBM05_mito","France_FBM05", reg_tree$tip.label)
reg_tree$tip.label<-gsub("FBM06_mito","France_FBM06", reg_tree$tip.label)
reg_tree$tip.label<-gsub("FBM09_mito","France_FBM09", reg_tree$tip.label)
reg_tree$tip.label<-gsub("D1_mito","England_D1", reg_tree$tip.label)
reg_tree$tip.label<-gsub("G1_mito","England_G1", reg_tree$tip.label)
reg_tree$tip.label<-gsub("G2_mito","England_G2", reg_tree$tip.label)
reg_tree$tip.label<-gsub("KB_mito","England_KB", reg_tree$tip.label)
reg_tree$tip.label<-gsub("BK01_mito","England_BK01", reg_tree$tip.label)
reg_tree$tip.label<-gsub("A02_mito","Alaska_A02", reg_tree$tip.label)
reg_tree$tip.label<-gsub("A09_mito","Alaska_A09", reg_tree$tip.label)
reg_tree$tip.label<-gsub("A19_mito","Alaska_A19", reg_tree$tip.label)
reg_tree$tip.label<-gsub("CAN49_mito","Canada_CAN49", reg_tree$tip.label)
reg_tree$tip.label<-gsub("CAN50_mito","Canada_CAN50", reg_tree$tip.label)
reg_tree$tip.label<-gsub("CAN52_mito","Canada_CAN52", reg_tree$tip.label)
reg_tree$tip.label<-gsub("NZ034_mito","New.Zealand_NZ034", reg_tree$tip.label)
reg_tree$tip.label<-gsub("NZ108_mito","New.Zealand_NZ108", reg_tree$tip.label)
reg_tree$tip.label<-gsub("NZ224_mito","New.Zealand_NZ224", reg_tree$tip.label)
reg_tree$tip.label<-gsub("J296_mito","Japan_J296", reg_tree$tip.label)
#reg_tree$tip.label<-gsub("J387_mito","CladeB_Japan_J387", reg_tree$tip.label)
reg_tree$tip.label<-gsub("J390_mito","Japan_J390", reg_tree$tip.label)
reg_tree$tip.label<-gsub("J424_mito","Japan_J424", reg_tree$tip.label)
reg_tree$tip.label<-gsub("J431_mito","Japan_J431", reg_tree$tip.label)
reg_tree$tip.label<-gsub("GS01_mito","Germany_GS01", reg_tree$tip.label)
reg_tree$tip.label<-gsub("GS02_mito","Germany_GS02", reg_tree$tip.label)
reg_tree$tip.label<-gsub("GS03_mito","Germany_GS03", reg_tree$tip.label)
reg_tree$tip.label<-gsub("Ire2_mito","N.Ireland_Ire2", reg_tree$tip.label)
reg_tree$tip.label<-gsub("H2_mito","Wales_H2", reg_tree$tip.label)

##get tree data to manipulate bootstrap values (remove <50%) and add posterior probabilities
q2 <- ggtree(reg_tree)                                   #### tree specfic
b <- q2$data    # extract bootstraps in node metadata
#b$branch.length[b$branch.length < 0.0000020226] <- 0.0000000000
#d<-q2$data


# Specify group based on location
b$group<-b$label
b$group[b$isTip == FALSE] <- NA
b$group<-strsplit(as.character(b$group), "_")
b$group<-sapply(b$group, "[[",1)

# Add colours to use for each
b$color<-b$group
b$color<-gsub("Alaska","#FB9A99",b$color)
b$color<-gsub("Canada", "grey" , b$color)
b$color<-gsub("Japan", "#E31A1C", b$color)
b$color<-gsub("New.Zealand", "#6A3D9A", b$color)
b$color<-gsub("N.Ireland", "#33A02C", b$color)
b$color<-gsub("England", "#A6CEE3", b$color)
b$color<-gsub("Wales", "#1F78B4", b$color)
b$color<-gsub("Scotland", "#B2DF8A", b$color)
b$color<-gsub("Germany", "#FDBF6F", b$color)
b$color<-gsub("France", "#FF7F00", b$color)
#b$color<-gsub("CladeB", "black", b$color)


# Remove low bootstraps and tip labels
b$label[b$isTip == TRUE] <- NA  #keep only values for non-tip nodes
b$label[b$label == ""]<- NA   #change root node to a value (by default empty)
#b$label[b$label == "Root"]<- NA  
b$label<-as.numeric(b$label)   #change to numeric
b$label[b$label < 65] <- NA  #remove all bootstraps of less than 50%

# Add posterior probabilities                              #### tree specific
b$label[b$label == 75]<- "75/0.657"  ##FBM05-H2
b$label[b$label == 86]<- "86/0.841"  ##J431 & G2
b$label[b$label == 80]<- "80/0.791"  ##canadas
b$label[b$label == 68]<- "68/-"  ##sub-clade 2 from sub-clade 1
b$label[b$label == 76]<- "76/-" ##sub-clade 1
b$label[b$label == 66]<- "66/0.728"  ##Alaskas
b$label[b$label == 72]<-NA
b$label[b$label == 80]<-NA
b$label[b$label == 74]<-NA
b$label[b$label == 65]<-NA
b$label[b$label == 77]<-NA
b$label[b$label == 79]<-NA
b$label[b$label == 70]<-NA
b$label[b$label == 72]<-NA
b$label[b$node == 39]<-NA
b$label[b$node == 56]<- "55/-"

# Create a df of metadata for each location for colours
meta_data2<-as.data.frame(a[,c("group", "label")])
meta_data2$group<-as.character(meta_data2$group)
meta_data2<-meta_data2[1:length(alllabels),]
rownames(meta_data2)<-reg_tree$tip.label                #### tree specfic
meta_data2$label<-NA

# Plot tree with node support
p2<-ggtree(reg_tree) +                                   #### tree specfic
  geom_treescale()+
  geom_rootedge(rootedge = 0.001)+
  geom_tiplab(size =3, align=FALSE) +
  geom_text2(size = 3,aes(label=b$label), color="black", nudge_x = -0.0003, nudge_y = 0.2) + 
  geom_cladelabel(node=50, label="Sub-clade 2", align=TRUE) +
  geom_cladelabel(node=55, label="Sub-clade 3", align=TRUE, offset = 0.002) + 
  geom_strip('Germany_GS01', 'England_D1', label = "Sub-clade 1") 


p2<-ggtree(reg_tree) +                                   #### tree specfic
  geom_treescale(width = 0.001)+
  geom_rootedge(rootedge = 0.001)+
  geom_tiplab(size = 4, align=FALSE, family = "serif") +
  geom_text2(size = 4, aes(label=a$label), color="black", family = "serif", nudge_x = -0.0008, nudge_y = 0.2)

p<-ggtree(reg_tree) +
  geom_text2(size = 3,aes(label=d$node), color="black")

# Add colours for each location
hm2<-gheatmap(p2,meta_data2, offset = 0.002, width=0.1, font.size=3, colnames=FALSE, hjust = 0)+ scale_fill_manual(breaks = legend_ord, values=colval) +
  theme(legend.position = "left", text = element_text(family = "serif", size = 15), legend.title = element_blank(), legend.key = element_blank()) # no keys
plot(hm2)


###############################################################################################################################################

mito_tree<-ape::read.tree("Rooted_Mito_partitioned.contree")

mito_tree <- ape::drop.tip(mito_tree, "J387_mito")
#mito_tree<-root(mito_tree, outgroup = "J387_mito")

### Tree data
mito_tree$tip.label<-gsub("58A_mito","Scotland_58A", mito_tree$tip.label)
mito_tree$tip.label<-gsub("FOF1b_mito","Scotland_FOF1b", mito_tree$tip.label)
mito_tree$tip.label<-gsub("66A_mito","Scotland_66A", mito_tree$tip.label)
mito_tree$tip.label<-gsub("FBM05_mito","France_FBM05", mito_tree$tip.label)
mito_tree$tip.label<-gsub("FBM06_mito","France_FBM06", mito_tree$tip.label)
mito_tree$tip.label<-gsub("FBM09_mito","France_FBM09", mito_tree$tip.label)
mito_tree$tip.label<-gsub("D1_mito","England_D1", mito_tree$tip.label)
mito_tree$tip.label<-gsub("G1_mito","England_G1", mito_tree$tip.label)
mito_tree$tip.label<-gsub("G2_mito","England_G2", mito_tree$tip.label)
mito_tree$tip.label<-gsub("KB_mito","England_KB", mito_tree$tip.label)
mito_tree$tip.label<-gsub("BK01_mito","England_BK01", mito_tree$tip.label)
mito_tree$tip.label<-gsub("A02_mito","Alaska_A02", mito_tree$tip.label)
mito_tree$tip.label<-gsub("A09_mito","Alaska_A09", mito_tree$tip.label)
mito_tree$tip.label<-gsub("A19_mito","Alaska_A19", mito_tree$tip.label)
mito_tree$tip.label<-gsub("CAN49_mito","Canada_CAN49", mito_tree$tip.label)
mito_tree$tip.label<-gsub("CAN50_mito","Canada_CAN50", mito_tree$tip.label)
mito_tree$tip.label<-gsub("CAN52_mito","Canada_CAN52", mito_tree$tip.label)
mito_tree$tip.label<-gsub("NZ034_mito","New.Zealand_NZ034", mito_tree$tip.label)
mito_tree$tip.label<-gsub("NZ108_mito","New.Zealand_NZ108", mito_tree$tip.label)
mito_tree$tip.label<-gsub("NZ224_mito","New.Zealand_NZ224", mito_tree$tip.label)
mito_tree$tip.label<-gsub("J296_mito","Japan_J296", mito_tree$tip.label)
#mito_tree$tip.label<-gsub("J387_mito","CladeB_Japan_J387", mito_tree$tip.label)
mito_tree$tip.label<-gsub("J390_mito","Japan_J390", mito_tree$tip.label)
mito_tree$tip.label<-gsub("J424_mito","Japan_J424", mito_tree$tip.label)
mito_tree$tip.label<-gsub("J431_mito","Japan_J431", mito_tree$tip.label)
mito_tree$tip.label<-gsub("GS01_mito","Germany_GS01", mito_tree$tip.label)
mito_tree$tip.label<-gsub("GS02_mito","Germany_GS02", mito_tree$tip.label)
mito_tree$tip.label<-gsub("GS03_mito","Germany_GS03", mito_tree$tip.label)
mito_tree$tip.label<-gsub("Ire2_mito","N.Ireland_Ire2", mito_tree$tip.label)
mito_tree$tip.label<-gsub("H2_mito","Wales_H2", mito_tree$tip.label)

##get tree data to manipulate bootstrap values (remove <50%) and add posterior probabilities
q3 <- ggtree(mito_tree)                                   #### tree specfic
c <- q3$data    # extract bootstraps in node metadata
d<-q3$data

# Specify group based on location
c$group<-c$label
c$group[c$isTip == FALSE] <- NA
c$group<-strsplit(as.character(c$group), "_")
c$group<-sapply(c$group, "[[",1)

# Add colours to use for each
c$color<-c$group
c$color<-gsub("Alaska","#FB9A99",c$color)
c$color<-gsub("Canada", "grey" , c$color)
c$color<-gsub("Japan", "#E31A1C", c$color)
c$color<-gsub("New.Zealand", "#6A3D9A", c$color)
c$color<-gsub("N.Ireland", "#33A02C", c$color)
c$color<-gsub("England", "#A6CEE3", c$color)
c$color<-gsub("Wales", "#1F78B4", c$color)
c$color<-gsub("Scotland", "#B2DF8A", c$color)
c$color<-gsub("Germany", "#FDBF6F", c$color)
c$color<-gsub("France", "#FF7F00", c$color)
#c$color<-gsub("CladeB", "black", c$color)


# Remove low bootstraps and tip labels
c$label[c$isTip == TRUE] <- NA  #keep only values for non-tip nodes
c$label[c$label == ""]<- NA   #change root node to a value (by default empty)
#c$label[c$label == "Root"]<- NA  
c$label<-as.numeric(c$label)   #change to numeric
c$label[c$label < 65] <- NA  #remove all bootstraps of less than 50%

# Add posterior probabilities                              #### tree specific
c$label[c$label ==99]<-"99/1.000" #sub 3, sub 2 & alaskas
c$label[c$label ==62]<-"62/0.638" #H2-Ire2-FBM05
c$label[c$label ==100]<- "100/1.000" #J431
c$label[c$label ==92]<- "92/1.000" #canadas
c$label[c$label ==72]<-NA
c$label[c$label ==78]<-"78/0.521"
c$label[c$label ==79]<-NA
c$label[c$label ==97]<- "97/1.000" #sub1
c$label[c$label ==86]<-"86/0.983" #alaska-nz
c$label[c$label ==91]<-"91/1.000"  ##d1-kb
c$label[c$label ==55]<-NA
c$label[c$node ==57]<- "-/0.614"
c$label[c$node ==56]<- "62/0.638"


# Create a df of metadata for each location for colours
meta_data3<-as.data.frame(a[,c("group", "label")])
meta_data3$group<-as.character(meta_data3$group)
meta_data3<-meta_data3[1:length(alllabels),]
rownames(meta_data3)<-mito_tree$tip.label                #### tree specfic
meta_data3$label<-NA

# Plot tree with node support
p3<-ggtree(mito_tree) +                                   #### tree specfic
  geom_treescale()+
  geom_rootedge(rootedge = 0.001)+
  geom_tiplab(size =3, align=FALSE) +
  geom_text2(size = 3,aes(label=c$label), color="black", nudge_x = -0.0002, nudge_y = 0.2) + 
  geom_cladelabel(node=49, label="Sub-clade 2", align=TRUE, offset = 0.002) +
  geom_cladelabel(node=55, label="Sub-clade 3", align=TRUE, offset = 0.002) + 
  geom_cladelabel(node=32, label="Sub-clade 1", align=TRUE, offset = 0.002)  
#+ geom_strip('New.Zealand_NZ108', 'Japan_J296', label = "Sub-clade 1") 



p3<-ggtree(mito_tree) +                                   #### tree specfic
  geom_treescale(width = 0.001)+
  geom_rootedge(rootedge = 0.001)+
  geom_tiplab(size = 4, align=FALSE, family = "serif") +
  geom_text2(size = 4, aes(label=a$label), color="black", family = "serif", nudge_x = -0.0008, nudge_y = 0.2)


p<-ggtree(mito_tree) +
  geom_text2(size = 3,aes(label=d$node), color="black")


# Add colours for each location
hm3<-gheatmap(p3,meta_data3, offset = 0.002, width=0.1, font.size=3, colnames=FALSE, hjust = 0)+ scale_fill_manual(breaks = legend_ord, values=colval) +
  theme(legend.position = "left", text = element_text(family = "serif", size = 15),legend.title = element_blank(), legend.key = element_blank()) # no keys
plot(hm3)




###############################################################################################################################################


mito_tree<-ape::read.tree("Rooted_Mito_partitioned.contree")

mito_tree <- ape::drop.tip(mito_tree, "J387_mito")
#mito_tree<-root(mito_tree, outgroup = "J387_mito")

### Tree data
mito_tree$tip.label<-gsub("58A_mito","Scotland_58A", mito_tree$tip.label)
mito_tree$tip.label<-gsub("FOF1b_mito","Scotland_FOF1b", mito_tree$tip.label)
mito_tree$tip.label<-gsub("66A_mito","Scotland_66A", mito_tree$tip.label)
mito_tree$tip.label<-gsub("FBM05_mito","France_FBM05", mito_tree$tip.label)
mito_tree$tip.label<-gsub("FBM06_mito","France_FBM06", mito_tree$tip.label)
mito_tree$tip.label<-gsub("FBM09_mito","France_FBM09", mito_tree$tip.label)
mito_tree$tip.label<-gsub("D1_mito","England_D1", mito_tree$tip.label)
mito_tree$tip.label<-gsub("G1_mito","England_G1", mito_tree$tip.label)
mito_tree$tip.label<-gsub("G2_mito","England_G2", mito_tree$tip.label)
mito_tree$tip.label<-gsub("KB_mito","England_KB", mito_tree$tip.label)
mito_tree$tip.label<-gsub("BK01_mito","England_BK01", mito_tree$tip.label)
mito_tree$tip.label<-gsub("A02_mito","Alaska_A02", mito_tree$tip.label)
mito_tree$tip.label<-gsub("A09_mito","Alaska_A09", mito_tree$tip.label)
mito_tree$tip.label<-gsub("A19_mito","Alaska_A19", mito_tree$tip.label)
mito_tree$tip.label<-gsub("CAN49_mito","Canada_CAN49", mito_tree$tip.label)
mito_tree$tip.label<-gsub("CAN50_mito","Canada_CAN50", mito_tree$tip.label)
mito_tree$tip.label<-gsub("CAN52_mito","Canada_CAN52", mito_tree$tip.label)
mito_tree$tip.label<-gsub("NZ034_mito","New.Zealand_NZ034", mito_tree$tip.label)
mito_tree$tip.label<-gsub("NZ108_mito","New.Zealand_NZ108", mito_tree$tip.label)
mito_tree$tip.label<-gsub("NZ224_mito","New.Zealand_NZ224", mito_tree$tip.label)
mito_tree$tip.label<-gsub("J296_mito","Japan_J296", mito_tree$tip.label)
#mito_tree$tip.label<-gsub("J387_mito","CladeB_Japan_J387", mito_tree$tip.label)
mito_tree$tip.label<-gsub("J390_mito","Japan_J390", mito_tree$tip.label)
mito_tree$tip.label<-gsub("J424_mito","Japan_J424", mito_tree$tip.label)
mito_tree$tip.label<-gsub("J431_mito","Japan_J431", mito_tree$tip.label)
mito_tree$tip.label<-gsub("GS01_mito","Germany_GS01", mito_tree$tip.label)
mito_tree$tip.label<-gsub("GS02_mito","Germany_GS02", mito_tree$tip.label)
mito_tree$tip.label<-gsub("GS03_mito","Germany_GS03", mito_tree$tip.label)
mito_tree$tip.label<-gsub("Ire2_mito","N.Ireland_Ire2", mito_tree$tip.label)
mito_tree$tip.label<-gsub("H2_mito","Wales_H2", mito_tree$tip.label)

##get tree data to manipulate bootstrap values (remove <50%) and add posterior probabilities
q3 <- ggtree(mito_tree)                                   #### tree specfic
c <- q3$data    # extract bootstraps in node metadata
d<-q3$data

# Specify group based on location
c$group<-c$label
c$group[c$isTip == FALSE] <- NA

# Add colours to use for each
c$group<-gsub("New.Zealand_NZ108","A",c$group)
c$group<-gsub("Alaska_A19", "A" , c$group)


c$group<-gsub("England_KB","B",c$group)
c$group<-gsub("England_D1", "B" , c$group)

c$group<-gsub("Japan_J296","C",c$group)

c$group<-gsub("Canada_CAN49", "D" , c$group)
c$group<-gsub("Canada_CAN50", "D" , c$group)
c$group<-gsub("Canada_CAN52", "D" , c$group)

c$group<-gsub("Japan_J431", "E" , c$group)
c$group<-gsub("England_G2", "E" , c$group)

c$group<-gsub("Alaska_A02", "F" , c$group)
c$group<-gsub("Alaska_A09", "F" , c$group)


c$group<-gsub("Wales_H2", "G" , c$group)
c$group<-gsub("N.Ireland_Ire2", "G" , c$group)
c$group<-gsub("France_FBM05", "G" , c$group)
c$group<-gsub("Japan_J424", "G" , c$group)


c$group<-gsub("Germany_GS03", "H" , c$group)
c$group<-gsub("England_BK01", "H" , c$group)
c$group<-gsub("Germany_GS02", "H" , c$group)
c$group<-gsub("Germany_GS01", "H" , c$group)
c$group<-gsub("New.Zealand_NZ034", "H" , c$group)
c$group<-gsub("New.Zealand_NZ224", "H" , c$group)
c$group<-gsub("France_FBM06", "H" , c$group)
c$group<-gsub("France_FBM09", "H" , c$group)
c$group<-gsub("Scotland_58A", "H" , c$group)
c$group<-gsub("Scotland_FOF1b", "H" , c$group)
c$group<-gsub("Scotland_66A", "H" , c$group)
c$group<-gsub("England_G1", "H" , c$group)
c$group<-gsub("Japan_J390", "H" , c$group)

subset1<-subset(c, !is.na(c$group))
dd<-data.frame(seq = subset1$label , cat = subset1$group)
dd$col <- dd$cat
dd$col[dd$col == "A"] <- compcol[1]
dd$col[dd$col == "B"] <- compcol[2]
dd$col[dd$col == "C"] <- compcol[3]
dd$col[dd$col == "D"] <- compcol[4]
dd$col[dd$col == "E"] <- compcol[5]
dd$col[dd$col == "F"] <- compcol[6]
dd$col[dd$col == "G"] <- compcol[7]
dd$col[dd$col == "H"] <- compcol[8]

# Remove low bootstraps and tip labels
c$label[c$isTip == TRUE] <- NA  #keep only values for non-tip nodes
c$label[c$label == ""]<- NA   #change root node to a value (by default empty)
#c$label[c$label == "Root"]<- NA  
c$label<-as.numeric(c$label)   #change to numeric
c$label[c$label < 65] <- NA  #remove all bootstraps of less than 50%

# Add posterior probabilities                              #### tree specific
c$label[c$label ==99]<-"*"
c$label[c$label ==62]<- NA
c$label[c$label ==100]<- "*" #J431
c$label[c$label ==92]<- "*" #canadas
c$label[c$label ==72]<- NA
c$label[c$label ==78]<-"*"
c$label[c$label ==79]<- NA
c$label[c$label ==97]<- "*" #sub1
c$label[c$label ==86]<-"*" #alaska-nz
c$label[c$label ==91]<-"*"  ##d1-kb
c$label[c$label ==55]<- NA
c$label[c$node ==57]<- NA
c$label[c$node ==56]<- NA



pm<-ggtree(mito_tree) +                                   #### tree specfic
  geom_rootedge(rootedge = 0.001)

pm <- pm  %<+% dd +
  geom_tiplab(aes(fill= col), geom="label", family = "serif", align=FALSE) +
  geom_text2(size = 5 ,aes(label=c$label), color="black", nudge_x = -0.00002)+
  xlim(-0.0001,0.005)+
  scale_fill_manual(values = compcol)+
  theme(legend.position="none")


#scales::show_col(RColorBrewer::brewer.pal(8, "Accent"))
#RColorBrewer::brewer.pal(8, "Set2")
#geom_label(aes(label= c$label , alpha=0.8, fill= c$group ))

##colourblind friendly palette: 
compcol<- c( '#0077BB', '#33BBEE', '#009988','#EE7733','#CC3311', '#EE3377',  '#BBBBBB', "white")
scales::show_col(compcol)



###############################################################################################################################################

coi_tree<-ape::read.tree("Rooted_TUN_COI.contree")     ##tree specific

coi_tree <- ape::drop.tip(coi_tree, "J387_mito")
#coi_tree<-root(coi_tree, outgroup = "J387_mito")

### Tree data
coi_tree$tip.label<-gsub("58A_mito","Scotland_58A", coi_tree$tip.label)
coi_tree$tip.label<-gsub("FOF1b_mito","Scotland_FOF1b", coi_tree$tip.label)
coi_tree$tip.label<-gsub("66A_mito","Scotland_66A", coi_tree$tip.label)
coi_tree$tip.label<-gsub("FBM05_mito","France_FBM05", coi_tree$tip.label)  ####
coi_tree$tip.label<-gsub("FBM06_mito","France_FBM06", coi_tree$tip.label)
coi_tree$tip.label<-gsub("FBM09_mito","France_FBM09", coi_tree$tip.label)
coi_tree$tip.label<-gsub("D1_mito","England_D1", coi_tree$tip.label)
coi_tree$tip.label<-gsub("G1_mito","England_G1", coi_tree$tip.label)
coi_tree$tip.label<-gsub("G2_mito","England_G2", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("KB_mito","England_KB", coi_tree$tip.label)
coi_tree$tip.label<-gsub("BK01_mito","England_BK01", coi_tree$tip.label)
coi_tree$tip.label<-gsub("A02_mito","Alaska_A02", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("A09_mito","Alaska_A09", coi_tree$tip.label)
coi_tree$tip.label<-gsub("A19_mito","Alaska_A19", coi_tree$tip.label)
coi_tree$tip.label<-gsub("CAN49_mito","Canada_CAN49", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("CAN50_mito","Canada_CAN50", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("CAN52_mito","Canada_CAN52", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("NZ034_mito","New.Zealand_NZ034", coi_tree$tip.label)
coi_tree$tip.label<-gsub("NZ108_mito","New.Zealand_NZ108", coi_tree$tip.label)
coi_tree$tip.label<-gsub("NZ224_mito","New.Zealand_NZ224", coi_tree$tip.label)
coi_tree$tip.label<-gsub("J296_mito","Japan_J296", coi_tree$tip.label)
#coi_tree$tip.label<-gsub("J387_mito","CladeB_Japan_J387", coi_tree$tip.label)
coi_tree$tip.label<-gsub("J390_mito","Japan_J390", coi_tree$tip.label)
coi_tree$tip.label<-gsub("J424_mito","Japan_J424", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("J431_mito","Japan_J431", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("GS01_mito","Germany_GS01", coi_tree$tip.label)
coi_tree$tip.label<-gsub("GS02_mito","Germany_GS02", coi_tree$tip.label)
coi_tree$tip.label<-gsub("GS03_mito","Germany_GS03", coi_tree$tip.label)
coi_tree$tip.label<-gsub("Ire2_mito","N.Ireland_Ire2", coi_tree$tip.label) ####
coi_tree$tip.label<-gsub("H2_mito","Wales_H2", coi_tree$tip.label) ####


###############################################################################################################################################

pc<-ggtree(coi_tree) +                                   #### tree specfic
  geom_rootedge(rootedge = 0.001)

pc<-pc %<+% dd +
  geom_tiplab(aes(fill= col), family = "serif", geom="label", align=FALSE) +
  scale_fill_manual(values = compcol)+
  theme(legend.position="none")


pr<-ggtree(reg_tree) +                                   #### reg tree specfic
  geom_rootedge(rootedge = 0.001)

pr<- pr  %<+% dd +
  geom_tiplab(aes(fill= col), family = "serif", geom="label", align=FALSE) +
  scale_fill_manual(values = compcol)+
  theme(legend.position="none")


cowplot::plot_grid(pm, pr, pc, ncol=3, labels=c("Mitochondrial Genome", "Candidate Region", "Partial-COI"), rel_widths = c(1,2,2) )
ggarrange(pm, pr, pc,  nrow = 1, widths = c(1, 2, 2), labels = c("A","B","C"), family = "serif")


##### Haplotype Networks
################# FUNCTION HAP.custom()

hap.custom<- function(fastafile, indel="5th", labels=NULL){
  #indels = argument in haplotypes::distance() on how to treat indels
  #labels = what to use to label each haplotype ID (pegas::haplotype), default = roman numerals
  
  ### Read FASTA file and calculate distance matrix ("N")  
  fas<-haplotypes::read.fas(fastafile)
  
  d<-distance(fas,indels=indel)
  if(length(d)==0) stop("at least two DNA sequences are required") #total number of comparisons
  
  
  ### Run modified haplotype function based on combined haplotypes and pegas packages
  
  ## Part 1=Haplotypes::haplotypes
  
  x <- as.matrix(d) #turn distance (class dist()) into matrix
  diag(x)<-0   #fill diagonal
  nseq<-nrow(x)  #total number of sequences
  whap<-x[1,]==0  #T/F of which in the top row have zero distance  
  haploind<-list(which(whap))   #column ID of which are the which are the same haplotype at 58A (as a list)
  haplovec<-which(whap)  ##column ID of which are the which are the same haplotype at 58A (as a vector)
  
  if(length(x)>1){   #total number of cells in matrix
    
    for(i in 2:nseq)  #repeat pairwise comparison for the rest of the rows (i.e. 66A-NZ224)
    { 
      whap<-x[i,]==0
      whap[haplovec]<-FALSE
      haploind<-c(haploind,list(which(whap)))  #combines list of haploind (column ID of ones the same at 58A), with ID's that are identical to the rest of the samples (list of 30 containing ID matches)
      haplovec<-unique(c(haplovec,which(whap)))  #c(haplovec, which(whap)) = combines vector for first row (58A) with all other sample comparisons - length=30 - and makes sure there are no duplicates
    }	
  }
  
  empthaplo<-sapply(haploind,length) #vector of frequencies (i.e. number of ID's in each item of the list)
  haploind<-unique(haploind[empthaplo>0])  #makes sure that you are only keeping unique comparisons
  haplolistint<-lapply(haploind, as.integer)
  haplolist<-lapply(haploind,names)   #list originally contains names + row number; this replaces that with just names
  
  uniqhapindex<-sapply(haploind,"[",1)  #unlists into the first value of every haplotype (Na's for remaining samples)
  
  hapdistmat<-as.matrix(x[uniqhapindex,uniqhapindex])  #turns into matrix of pairwise comparisons 
  
  
  freq<-sapply(haplolist,length) #vector of frequencies with number of ID's in each item of the list > 0
  hapnum<-length(freq)
  
  
  ## Part 2= change for pegas format:  class(obj) = "haplotype" "DNAbin"  
  
  dnabinFAS<-read.dna(fastafile, format="fasta")
  
  obj <- dnabinFAS[uniqhapindex, ] #extract 6 unique sequences
  
  if (is.null(labels)) labels <- as.character(as.roman(seq_along(haploind)))  #creates labels for each haplotype identified
  rownames(obj) <- labels  #label sequences by haplotype ID
  class(obj) <- c("haplotype", "DNAbin")
  #attr(obj, "index") <- lapply(i, function(x) which(h == x))
  attr(obj, "index") <- haplolistint
  
  nms.x <- deparse(substitute(fastafile))
  attr(obj, "from") <- nms.x
  obj
}

#fastafile = files[1]
#indel = "5th"
#labels = NULL

change.sampleID<-function(fastafiles, torm){
  seq<-read.dna(fastafiles, format="fasta")
  
  rownames(seq)<-gsub("58A_mito","Scotland_58A", rownames(seq))
  rownames(seq)<-gsub("FOF1b_mito","Scotland_FOF1b", rownames(seq))
  rownames(seq)<-gsub("66A_mito","Scotland_66A", rownames(seq))
  rownames(seq)<-gsub("FBM05_mito","France_FBM05", rownames(seq))
  rownames(seq)<-gsub("FBM06_mito","France_FBM06", rownames(seq))
  rownames(seq)<-gsub("FBM09_mito","France_FBM09", rownames(seq))
  rownames(seq)<-gsub("D1_mito","England_D1", rownames(seq))
  rownames(seq)<-gsub("G1_mito","England_G1", rownames(seq))
  rownames(seq)<-gsub("G2_mito","England_G2", rownames(seq))
  rownames(seq)<-gsub("KB_mito","England_KB", rownames(seq))
  rownames(seq)<-gsub("BK01_mito","England_BK01", rownames(seq))
  rownames(seq)<-gsub("A02_mito","Alaska_A02", rownames(seq))
  rownames(seq)<-gsub("A09_mito","Alaska_A09", rownames(seq))
  rownames(seq)<-gsub("A19_mito","Alaska_A19", rownames(seq))
  rownames(seq)<-gsub("CAN49_mito","Canada_CAN49", rownames(seq))
  rownames(seq)<-gsub("CAN50_mito","Canada_CAN50", rownames(seq))
  rownames(seq)<-gsub("CAN52_mito","Canada_CAN52", rownames(seq))
  rownames(seq)<-gsub("NZ034_mito","New.Zealand_NZ034", rownames(seq))
  rownames(seq)<-gsub("NZ108_mito","New.Zealand_NZ108", rownames(seq))
  rownames(seq)<-gsub("NZ224_mito","New.Zealand_NZ224", rownames(seq))
  rownames(seq)<-gsub("J296_mito","Japan_J296", rownames(seq))
  rownames(seq)<-gsub("J387_mito","CladeB_Japan_J387", rownames(seq))
  rownames(seq)<-gsub("J390_mito","Japan_J390", rownames(seq))
  rownames(seq)<-gsub("J424_mito","Japan_J424", rownames(seq))
  rownames(seq)<-gsub("J431_mito","Japan_J431", rownames(seq))
  rownames(seq)<-gsub("GS01_mito","Germany_GS01", rownames(seq))
  rownames(seq)<-gsub("GS02_mito","Germany_GS02", rownames(seq))
  rownames(seq)<-gsub("GS03_mito","Germany_GS03", rownames(seq))
  rownames(seq)<-gsub("Ire2_mito","N.Ireland_Ire2", rownames(seq))
  rownames(seq)<-gsub("H2_mito","Wales_H2", rownames(seq))
  
  seq
}
setwd("../")

files<-c("NoRef.Rooted_TUN_COI.fas","NoRef.Rooted_Region1.fas", "NoRef.Rooted_Mito.fas")

### COI
coi_hap<-hap.custom(fastafile = files[1], indel = "5th", labels = NULL)
coi_net<-pegas::haploNet(coi_hap)
attr(coi_net, "labels") <- c("3", "6", "5", "2", "1")
coi_seq<-change.sampleID(files[1])
ind.hap<-with(utils::stack(setNames(attr(coi_hap, "index"), rownames(coi_hap))), table(hap=ind, pop=rownames(coi_seq)[values]))
mydata<-as.data.frame(ind.hap)
loc<-mydata[mydata$Freq ==1,]
locations<-strsplit(as.character(loc$pop), "_")
locations<-sapply(locations, "[[",1)
#locations[locations=="clade"]<-"ref"
newhap_coi<-table(loc$hap, locations)
plot(coi_net, size=attr(coi_net, "freq"), scale.ratio = 2, cex = 1.5, pie=newhap_coi, bg=colval, show.mutation=2)  #option for legend = T
replot()
legend("topright", colnames(newhap_coi), col=colval, pch=19, cex = 1.5)
#rect(-14, -10, 10, 10, lty= "dashed") #left bottom right top
#rect(15, -5, 19, 2, lty= "dashed")


### Region
reg_hap<-hap.custom(fastafile = files[2], indel = "5th", labels = NULL)
reg_net<-pegas::haploNet(reg_hap)
reg_seq<-change.sampleID(files[2])
ind.hap<-with(utils::stack(setNames(attr(reg_hap, "index"), rownames(reg_hap))), table(hap=ind, pop=rownames(reg_seq)[values]))
mydata<-as.data.frame(ind.hap)
loc<-mydata[mydata$Freq ==1,]
locations<-strsplit(as.character(loc$pop), "_")
locations<-sapply(locations, "[[",1)
#locations[locations=="clade"]<-"ref"
newhap_reg<-table(loc$hap, locations)
plot(reg_net, size=attr(reg_net, "freq"), pie=newhap_reg, bg=colval, show.mutation=2, scale.ratio =3, cex = 1)  #option for legend = T
legend("topleft", colnames(newhap_reg), col=colval, pch=19, cex = 1.5)



### Mitogenome
mito_hap<-hap.custom(fastafile = files[3], indel = "5th", labels = NULL)
mito_net<-pegas::haploNet(mito_hap)
mito_seq<-change.sampleID(files[3])
ind.hap<-with(utils::stack(setNames(attr(mito_hap, "index"), rownames(mito_hap))), table(hap=ind, pop=rownames(mito_seq)[values]))
mydata<-as.data.frame(ind.hap)
loc<-mydata[mydata$Freq ==1,]
locations<-strsplit(as.character(loc$pop), "_")
locations<-sapply(locations, "[[",1)
#locations[locations=="clade"]<-"ref"
newhap_mito<-table(loc$hap, locations)
plot(mito_net, size=attr(mito_net, "freq"), scale.ratio = 0.3, cex = 0.8, pie=newhap_mito, bg=colval, show.mutation=2)  #option for legend = T
replot()

#$x
#[1]  -0.1473569   2.2984349  -0.2854542  -8.5532775  -7.2552227  -0.3012630  -3.1626759 -22.1970013
#[9]   4.1936200   8.0877842  -0.2594931 -20.0652501 -16.0291660   3.5038853 -18.0894808   1.3108098

#$y
#[1] -0.25224583 12.59368160 12.66160703  0.87527010 -0.03336821 20.64655637 -7.42920330 17.34282351
#[9] -6.54960295  0.77142572 17.56825390 17.34282351 17.42081440  6.53478871 17.36882047 17.57679619
legend("bottomleft", colnames(newhap_mito), col=colval, pch=19, cex = 1.5)



##### Mantel Tests 
### Calculate Geographical distance using WGS84 ellipsoid
setwd("IQTree_2020-06-11")



#coordinates<-coordinates[complete.cases(coordinates),]

samplesJP<-c("Japan_J296", "Japan_J390", "Japan_J424", "Japan_J431")

geodist.formatted<-function(coordinates, pairwiserm = vector(), samplesrm = vector()){
  
  coordinates<-read.csv(coordinates, header=TRUE)
  coordinates$Sequence.ID<-gsub("NZealand","New.Zealand", coordinates$Sequence.ID)
  coordinates$Sequence.ID<-gsub("Japan_J387","CladeB_Japan_J387", coordinates$Sequence.ID)
  coordinates$Sequence.ID<-gsub("NIreland_Ire2","N.Ireland_Ire2", coordinates$Sequence.ID)
  coordinates$Sequence.ID<-gsub("England_KA","England_KB", coordinates$Sequence.ID)
  rownames(coordinates)<-coordinates$Sequence.ID
  coordinates<-coordinates[order(coordinates$Sequence.ID),]  #order alphabetically
  
  latlon_only<-as.matrix(coordinates[,c("Longitude","Latitude")])    ##geosphere::distGeo first column = longitude, second column = latitude
  geodistmat<-geosphere::distm(latlon_only) #uses WGS84 ellipsoid to calculate distances by default, output in meters
  rownames(geodistmat)<-rownames(latlon_only)
  colnames(geodistmat)<-rownames(latlon_only)
  
  if(length(pairwiserm) > 0){
    dx <- which(rownames(geodistmat) %in% pairwiserm)
    geodistmat<-geodistmat[-dx,]
  }
  
  if(length(samplesrm) > 0){
    rx <- which(rownames(geodistmat) %in% samplesrm)
    geodistmat<-geodistmat[-rx,]
    
    rx2<- which(colnames(geodistmat) %in% samplesrm)
    geodistmat<-geodistmat[,-rx2]
  }
  
  geodistmat
}

geodistmat<-geodist.formatted("../Geo_coordinates.csv")
geodistmatA<-geodist.formatted("../Geo_coordinates.csv", samplesrm = c("CladeB_Japan_J387") )
geodistmatJP<-geodist.formatted("../Geo_coordinates.csv", pairwiserm = samplesJP, samplesrm = c("CladeB_Japan_J387"))

### Cartesian geographic coordinates (XY instead of pairwise matrix)

#coordinates$x = 6371* cos(coordinates$Latitude) * cos(coordinates$Longitude)
#coordinates$y = 6371* cos(coordinates$Latitude) * sin(coordinates$Longitude)
#coordinates$z = 6371* sin(coordinates$Latitude)

#coordinatesxy<-coordinates[,5:6]


### Calculate genetic distance using trees:

gendist.formatted<-function(contree, pairwiserm = vector(), samplesrm = vector()){
  
  tree<-ape::read.tree(contree)
  
  tree$tip.label<-gsub("58A_mito","Scotland_58A", tree$tip.label)
  tree$tip.label<-gsub("FOF1b_mito","Scotland_FOF1b", tree$tip.label)
  tree$tip.label<-gsub("66A_mito","Scotland_66A", tree$tip.label)
  tree$tip.label<-gsub("FBM05_mito","France_FBM05", tree$tip.label)
  tree$tip.label<-gsub("FBM06_mito","France_FBM06", tree$tip.label)
  tree$tip.label<-gsub("FBM09_mito","France_FBM09", tree$tip.label)
  tree$tip.label<-gsub("D1_mito","England_D1", tree$tip.label)
  tree$tip.label<-gsub("G1_mito","England_G1", tree$tip.label)
  tree$tip.label<-gsub("G2_mito","England_G2", tree$tip.label)
  tree$tip.label<-gsub("KB_mito","England_KB", tree$tip.label)
  tree$tip.label<-gsub("BK01_mito","England_BK01", tree$tip.label)
  tree$tip.label<-gsub("A02_mito","Alaska_A02", tree$tip.label)
  tree$tip.label<-gsub("A09_mito","Alaska_A09", tree$tip.label)
  tree$tip.label<-gsub("A19_mito","Alaska_A19", tree$tip.label)
  tree$tip.label<-gsub("CAN49_mito","Canada_CAN49", tree$tip.label)
  tree$tip.label<-gsub("CAN50_mito","Canada_CAN50", tree$tip.label)
  tree$tip.label<-gsub("CAN52_mito","Canada_CAN52", tree$tip.label)
  tree$tip.label<-gsub("NZ034_mito","New.Zealand_NZ034", tree$tip.label)
  tree$tip.label<-gsub("NZ108_mito","New.Zealand_NZ108", tree$tip.label)
  tree$tip.label<-gsub("NZ224_mito","New.Zealand_NZ224", tree$tip.label)
  tree$tip.label<-gsub("J296_mito","Japan_J296", tree$tip.label)
  tree$tip.label<-gsub("J387_mito","CladeB_Japan_J387", tree$tip.label)
  tree$tip.label<-gsub("J390_mito","Japan_J390", tree$tip.label)
  tree$tip.label<-gsub("J424_mito","Japan_J424", tree$tip.label)
  tree$tip.label<-gsub("J431_mito","Japan_J431", tree$tip.label)
  tree$tip.label<-gsub("GS01_mito","Germany_GS01", tree$tip.label)
  tree$tip.label<-gsub("GS02_mito","Germany_GS02", tree$tip.label)
  tree$tip.label<-gsub("GS03_mito","Germany_GS03", tree$tip.label)
  tree$tip.label<-gsub("Ire2_mito","N.Ireland_Ire2", tree$tip.label)
  tree$tip.label<-gsub("H2_mito","Wales_H2", tree$tip.label)
  
  gendist<-ape::cophenetic.phylo(tree)       #calculates genetic distance
  
  
  if(length(pairwiserm) > 0){                    #options to remove Japan samples
    ex <- which(rownames(gendist) %in% pairwiserm)
    gendist<-gendist[-ex,]
  }
  
  if(length(samplesrm) > 0){
    rx <- which(rownames(gendist) %in% samplesrm)
    gendist<-gendist[-rx,]
    
    rx2<- which(colnames(gendist) %in% samplesrm)
    gendist<-gendist[,-rx2]
  }
  
  
  
  distmat<-as.data.frame(gendist)
  distmat<-distmat[order(rownames(gendist)),]
  distmat<-distmat[,order(colnames(gendist))]
  
  as.matrix(distmat)
}

## All samples
coi_gendistmat<- gendist.formatted("Rooted_TUN_COI.contree") 
win_gendistmat<- gendist.formatted("Rooted_Win30.contree")
amp1_gendistmat<- gendist.formatted("Rooted_Amp1.contree")
amp2_gendistmat<- gendist.formatted("Rooted_Amp2.contree")
reg_gendistmat<- gendist.formatted("Rooted_Reg1.contree")
mito_gendistmat<- gendist.formatted("Rooted_Mito_partitioned.contree")

## Without Clade B
coi_gendistmatA<- gendist.formatted("Rooted_TUN_COI.contree", samplesrm= c("CladeB_Japan_J387")) 
win_gendistmatA<- gendist.formatted("Rooted_Win30.contree", samplesrm= c("CladeB_Japan_J387"))
amp1_gendistmatA<- gendist.formatted("Rooted_Amp1.contree", samplesrm= c("CladeB_Japan_J387"))
amp2_gendistmatA<- gendist.formatted("Rooted_Amp2.contree", samplesrm= c("CladeB_Japan_J387"))
reg_gendistmatA<- gendist.formatted("Rooted_Reg1.contree", samplesrm= c("CladeB_Japan_J387"))
mito_gendistmatA<- gendist.formatted("Rooted_Mito_partitioned.contree", samplesrm= c("CladeB_Japan_J387"))


## Distances to Japan without Clade B
coi_gendistmatJP<- gendist.formatted("Rooted_TUN_COI.contree",pairwiserm = samplesJP, samplesrm = c("CladeB_Japan_J387")) 
win_gendistmatJP<- gendist.formatted("Rooted_Win30.contree", pairwiserm = samplesJP, samplesrm = c("CladeB_Japan_J387"))
amp1_gendistmatJP<- gendist.formatted("Rooted_Amp1.contree", pairwiserm = samplesJP, samplesrm = c("CladeB_Japan_J387"))
amp2_gendistmatJP<- gendist.formatted("Rooted_Amp2.contree", pairwiserm = samplesJP, samplesrm = c("CladeB_Japan_J387"))
reg_gendistmatJP<- gendist.formatted("Rooted_Reg1.contree", pairwiserm = samplesJP, samplesrm = c("CladeB_Japan_J387"))
mito_gendistmatJP<- gendist.formatted("Rooted_Mito_partitioned.contree", pairwiserm = samplesJP, samplesrm = c("CladeB_Japan_J387"))

## Make sure matrix names are the same: (columns and rows) These should all be true
colnames(coi_gendistmat) == colnames(geodistmat)
rownames(coi_gendistmat) == rownames(geodistmat)

colnames(win_gendistmat) == colnames(geodistmat)
rownames(win_gendistmat) == rownames(geodistmat)

colnames(amp1_gendistmat) == colnames(geodistmat)
rownames(amp1_gendistmat) == rownames(geodistmat)

colnames(amp2_gendistmat) == colnames(geodistmat)
rownames(amp2_gendistmat) == rownames(geodistmat)

colnames(reg_gendistmat) == colnames(geodistmat)
rownames(reg_gendistmat) == rownames(geodistmat)

colnames(mito_gendistmat) == colnames(geodistmat)
rownames(mito_gendistmat) == rownames(geodistmat)



###Mantel test (unstratified using Pearson correlation test) into a dataframe
mantel_df <- data.frame(Dataset=as.character(c("TUN_COI", "Window", "Amplicon1", "Amplicon2", "Region", "CompMito", "TUN_COI_A", "Window_A", "Amplicon1_A", "Amplicon2_A", "Region_A", "CompMito_A", "TUN_COI_JP", "Window_JP", "Amplicon1_JP", "Amplicon2_JP", "Region_JP", "CompMito_JP")),
                        Mantel.R=as.numeric(rep(NA, 18)),
                        P.value=as.numeric(rep(NA, 18)),
                        stringsAsFactors=FALSE)

vegan_mantel_coi<-vegan::mantel(as.dist(coi_gendistmat), as.dist(geodistmat), permutations  = 10000)
vegan_mantel_win<-vegan::mantel(as.dist(win_gendistmat), as.dist(geodistmat), permutations  = 10000)
vegan_mantel_amp1<-vegan::mantel(as.dist(amp1_gendistmat), as.dist(geodistmat), permutations  = 10000)
vegan_mantel_amp2<-vegan::mantel(as.dist(amp2_gendistmat), as.dist(geodistmat), permutations  = 10000)
vegan_mantel_reg<-vegan::mantel(as.dist(reg_gendistmat), as.dist(geodistmat), permutations  = 10000)
vegan_mantel_mit<-vegan::mantel(as.dist(mito_gendistmat), as.dist(geodistmat), permutations  = 10000)

mantel_df$Mantel.R[1]<-vegan_mantel_coi$statistic
mantel_df$P.value[1]<-vegan_mantel_coi$signif

mantel_df$Mantel.R[2]<-vegan_mantel_win$statistic
mantel_df$P.value[2]<-vegan_mantel_win$signif

mantel_df$Mantel.R[3]<-vegan_mantel_amp1$statistic
mantel_df$P.value[3]<-vegan_mantel_amp1$signif

mantel_df$Mantel.R[4]<-vegan_mantel_amp2$statistic
mantel_df$P.value[4]<-vegan_mantel_amp2$signif

mantel_df$Mantel.R[5]<-vegan_mantel_reg$statistic
mantel_df$P.value[5]<-vegan_mantel_reg$signif

mantel_df$Mantel.R[6]<-vegan_mantel_mit$statistic
mantel_df$P.value[6]<-vegan_mantel_mit$signif


vegan_mantel_coi<-vegan::mantel(as.dist(coi_gendistmatA), as.dist(geodistmatA), permutations  = 10000)
vegan_mantel_win<-vegan::mantel(as.dist(win_gendistmatA), as.dist(geodistmatA), permutations  = 10000)
vegan_mantel_amp1<-vegan::mantel(as.dist(amp1_gendistmatA), as.dist(geodistmatA), permutations  = 10000)
vegan_mantel_amp2<-vegan::mantel(as.dist(amp2_gendistmatA), as.dist(geodistmatA), permutations  = 10000)
vegan_mantel_reg<-vegan::mantel(as.dist(reg_gendistmatA), as.dist(geodistmatA), permutations  = 10000)
vegan_mantel_mit<-vegan::mantel(as.dist(mito_gendistmatA), as.dist(geodistmatA), permutations  = 10000)

mantel_df$Mantel.R[7]<-vegan_mantel_coi$statistic
mantel_df$P.value[7]<-vegan_mantel_coi$signif

mantel_df$Mantel.R[8]<-vegan_mantel_win$statistic
mantel_df$P.value[8]<-vegan_mantel_win$signif

mantel_df$Mantel.R[9]<-vegan_mantel_amp1$statistic
mantel_df$P.value[9]<-vegan_mantel_amp1$signif

mantel_df$Mantel.R[10]<-vegan_mantel_amp2$statistic
mantel_df$P.value[10]<-vegan_mantel_amp2$signif

mantel_df$Mantel.R[11]<-vegan_mantel_reg$statistic
mantel_df$P.value[11]<-vegan_mantel_reg$signif

mantel_df$Mantel.R[12]<-vegan_mantel_mit$statistic
mantel_df$P.value[12]<-vegan_mantel_mit$signif


vegan_mantel_coi<-vegan::mantel(as.dist(coi_gendistmatJP), as.dist(geodistmatJP), permutations  = 10000)
vegan_mantel_win<-vegan::mantel(as.dist(win_gendistmatJP), as.dist(geodistmatJP), permutations  = 10000)
vegan_mantel_amp1<-vegan::mantel(as.dist(amp1_gendistmatJP), as.dist(geodistmatJP), permutations  = 10000)
vegan_mantel_amp2<-vegan::mantel(as.dist(amp2_gendistmatJP), as.dist(geodistmatJP), permutations  = 10000)
vegan_mantel_reg<-vegan::mantel(as.dist(reg_gendistmatJP), as.dist(geodistmatJP), permutations  = 10000)
vegan_mantel_mit<-vegan::mantel(as.dist(mito_gendistmatJP), as.dist(geodistmatJP), permutations  = 10000)

mantel_df$Mantel.R[13]<-vegan_mantel_coi$statistic
mantel_df$P.value[13]<-vegan_mantel_coi$signif

mantel_df$Mantel.R[14]<-vegan_mantel_win$statistic
mantel_df$P.value[14]<-vegan_mantel_win$signif

mantel_df$Mantel.R[15]<-vegan_mantel_amp1$statistic
mantel_df$P.value[15]<-vegan_mantel_amp1$signif

mantel_df$Mantel.R[16]<-vegan_mantel_amp2$statistic
mantel_df$P.value[16]<-vegan_mantel_amp2$signif

mantel_df$Mantel.R[17]<-vegan_mantel_reg$statistic
mantel_df$P.value[17]<-vegan_mantel_reg$signif

mantel_df$Mantel.R[18]<-vegan_mantel_mit$statistic
mantel_df$P.value[18]<-vegan_mantel_mit$signif

#write.csv(mantel_df, file="MantelResults_2020-07-07.csv", row.names=FALSE)


mantel.formatted<-function(distmat, value = "Dist"){
  distmat[upper.tri(distmat)]<-NA
  distmelted<-reshape2::melt(distmat)
  colnames(distmelted)<-c("A", "B", value)
  distmelted$A<-as.character(distmelted$A)
  distmelted$B<-as.character(distmelted$B)
  distmelted<-na.omit(distmelted)
  distmelted
}

geodistmelted <- mantel.formatted(geodistmat)
coi_gendistmelted <- mantel.formatted(coi_gendistmat, value = "COI")
reg_gendistmelted <- mantel.formatted(reg_gendistmat, value = "Region")
mito_gendistmelted <- mantel.formatted(mito_gendistmat, value = "Mito")


### double check that labels all match
unique(coi_gendistmelted$A) == unique(geodistmelted$A)
#unique(win_gendistmelted$A) == unique(geodistmelted$A)
#unique(amp1_gendistmelted$A) == unique(geodistmelted$A)
#unique(amp2_gendistmelted$A) == unique(geodistmelted$A)
unique(reg_gendistmelted$A) == unique(geodistmelted$A)
unique(mito_gendistmelted$A) == unique(geodistmelted$A)



geodistmeltedA <- mantel.formatted(geodistmatA)
coi_gendistmeltedA <- mantel.formatted(coi_gendistmatA, value = "COI")
reg_gendistmeltedA <- mantel.formatted(reg_gendistmatA, value = "Region")
mito_gendistmeltedA <- mantel.formatted(mito_gendistmatA, value = "Mito")


geodistmeltedJP <- mantel.formatted(geodistmatJP)
coi_gendistmeltedJP <- mantel.formatted(coi_gendistmatJP, value = "COI")
reg_gendistmeltedJP <- mantel.formatted(reg_gendistmatJP, value = "Region")
mito_gendistmeltedJP <- mantel.formatted(mito_gendistmatJP, value = "Mito")



merge.plot<-function(mylist, mantelresult){
  
  ### Merge data frames
  
  mergeddf<-plyr::join_all(mylist, by = c("A", "B"))   #465 obs.
  
  ### change m to km
  mergeddf$Dist<-mergeddf$Dist*0.001
  
  dx <- which(!(colnames(mergeddf) %in% c("A","B")))  #only columns of distances
  mergeddf2<-mergeddf[, dx]
  
  melted_all<-reshape2::melt(mergeddf2, id= "Dist")
  
  
  p<-ggplot(melted_all, aes(x=Dist, y=value)) + 
    geom_point() +
    geom_tile(colour = "gray20") +
    xlab("Geographic Distance (km)") +
    ylab(expression(paste(italic("D")," (substitutions/site)"))) +
    geom_text(mapping = aes(x = Inf, y = Inf, label = paste("r[M] == ", mantelresult)),  parse=T, hjust = 1, vjust = 1, check_overlap = T) +
    theme_tufte()+
    theme(text = element_text(size=20), axis.line= element_line(size = 0.5, color = "black"))
  
  p
}

mantel_df$Mantel.R<-round(mantel_df$Mantel.R, digits = 3)
mantel_df$P.value<-round(mantel_df$P.value, digits=3)
mantel_df$label<- paste("'", mantel_df$Mantel.R, ", p = ", mantel_df$P.value, "'", sep="")


#mylist<-list(geodistmelted, coi_gendistmelted, win_gendistmelted, amp1_gendistmelted, amp2_gendistmelted, reg_gendistmelted, mito_gendistmelted)
mylist<-list(geodistmeltedA, coi_gendistmeltedA)
mylist2<- list (geodistmeltedA, reg_gendistmeltedA)
mylist3<- list(geodistmeltedA, mito_gendistmeltedA)


mylistJP<-list(geodistmeltedJP, coi_gendistmeltedJP)
mylistJP2<- list (geodistmeltedJP, reg_gendistmeltedJP)
mylistJP3<- list(geodistmeltedJP, mito_gendistmeltedJP)


mylist_all<-list(geodistmelted, coi_gendistmelted)
pCOI<-merge.plot(mylist = mylist_all, mantelresult = mantel_df$label[ mantel_df$Dataset == "TUN_COI"])


pCOIA<-merge.plot(mylist = mylist, mantelresult = mantel_df$label[mantel_df$Dataset == "TUN_COI_A"])
pRegA<-merge.plot(mylist = mylist2, mantelresult = mantel_df$label[ mantel_df$Dataset == "Region_A"])
pMitoA<-merge.plot(mylist = mylist3, mantelresult = mantel_df$label[ mantel_df$Dataset == "CompMito_A"])

pCOIJP<-merge.plot(mylist = mylistJP , mantelresult = mantel_df$label[ mantel_df$Dataset == "TUN_COI_JP"])
pRegJP<-merge.plot(mylist = mylistJP2, mantelresult = mantel_df$label[ mantel_df$Dataset == "Region_JP"])
pMitoJP<-merge.plot(mylist = mylistJP3, mantelresult = mantel_df$label[ mantel_df$Dataset == "CompMito_JP"])

margin = theme(plot.margin = unit(c(0.2,0.3,0.2,0.3), "cm"))  #top, right, bottom, left

pl<- list(pCOIA, pCOIJP)
pl2<- list(pRegA, pRegJP)
pl3<- list(pMitoA, pMitoJP)


library(ggpubr)
ggarrange(pCOIA,NULL, pCOIJP, nrow = 1, widths = c(1, 0.05, 1), labels = c("A","","B"))
ggarrange(pRegA,NULL, pRegJP, nrow = 1, widths = c(1, 0.05, 1), labels = c("A","","B"))
ggarrange(pMitoA,NULL, pMitoJP, nrow = 1, widths = c(1, 0.05, 1), labels = c("A","","B"), family = "serif")




## svg saved at 1500 x 900
##IBD_COI_2021-09-10

