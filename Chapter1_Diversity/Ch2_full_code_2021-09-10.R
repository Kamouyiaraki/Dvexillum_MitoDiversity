library(ggplot2)
library(seqinr)
library(vhica)
library(plot3D)
library(ggthemes)
library(reshape2)
library(dplyr)
library(scales)
library(seqmagick) 
library(ggmsa) 
library(patchwork)
library(gggenes)
library(alignfigR)
library(sjPlot)
library(gggenes)
library(circlize)
library(ape)
library(pegas)
library(seqRFLP)
library(plotrix)
library(adegenet)
library(Biostrings)
library(tidyr)
library(cowplot)


setwd("C:/Users/Maraki/Desktop/Shared_Drive_2020/D.vexillum_Ch2_mitodiversity")

# as.integer(ape::as.DNAbin("-")) #4
# as.integer(ape::as.DNAbin("a")) #136
# as.integer(ape::as.DNAbin("g")) #72
# as.integer(ape::as.DNAbin("c")) #40
# as.integer(ape::as.DNAbin("t")) #24
# as.integer(ape::as.DNAbin("n")) #240


FastaAli <- alignfigR::read_alignment("Mito_aln_cladeA.fas")
names(FastaAli)<-c("Scotland_58A", 
                   "Scotland_66A",
                   "Alaska_A02",
                   "Alaska_A09",
                   "Alaska_A19",
                   "England_BK01",
                   "Canada_CAN49",
                   "Canada_CAN50",
                   "Canada_CAN52",    
                   "England_D1",   
                   "France_FBM05",
                   "France_FBM06",
                   "France_FBM09",
                   "Scotland_FOF1b",
                   "England_G1",
                   "England_G2", 
                   "Germany_GS01",
                   "Germany_GS02",
                   "Germany_GS03",   
                   "Wales_H2",      
                   "N.Ireland_Ire2",   
                   "Japan_J269",
                   "Japan_J390",
                   "Japan_J424",
                   "Japan_J431",
                   "England_KB",   
                   "New Zealand_NZ034",
                   "New Zealand_NZ108",
                   "New Zealand_NZ224"
)

orderx<-c("Alaska_A02", "Alaska_A09", "Alaska_A19", "Canada_CAN49", "Canada_CAN50", "Canada_CAN52", "England_D1", "England_G1","England_G2","England_BK01","England_KB", "Wales_H2" ,"N.Ireland_Ire2", "Scotland_58A","Scotland_66A","Scotland_FOF1b","France_FBM05","France_FBM06","France_FBM09","Germany_GS01","Germany_GS02","Germany_GS03","Japan_J269","Japan_J390","Japan_J424","Japan_J431", "New Zealand_NZ034","New Zealand_NZ108","New Zealand_NZ224" )
FastaAli<-FastaAli[orderx]

ntcols <-
  c(
    "A" = "#2F2F2F",
    "C" = "#2F2F2F",
    "G" = "#2F2F2F" ,
    "T" = "#2F2F2F",
    "t" = "#2F2F2F",
    "-" = "#2F2F2F",
    "K" = "#2F2F2F",
    "n" = NULL,
    "N" = NULL,
    "R" = "#2F2F2F",
    "S" = "#2F2F2F",
    "W" = "#2F2F2F",
    "Y" = "#2F2F2F"
  )

# Extract desired alignment subset
plot_frame <- alignfigR::extract_subalign(FastaAli)

paln<- ggplot() + 
  geom_rect(plot_frame, mapping=aes(xmin=x1-1, xmax=x2-1, ymin = y1-1, ymax=y2-1.1, fill = seq), linetype=0) +
  geom_path(size= 0.8)+
  scale_fill_manual(values=ntcols, name = "") +
  scale_y_discrete(limits = rev(names(FastaAli)))+
  scale_x_continuous(breaks = c(1, seq(1000, 13438, by = 1000)))+
  geom_hline(yintercept=seq(0.5,29,1), color="#2F2F2F", size=0.5)+
  geom_segment(aes(y = -0.5, x=10004, yend=-0.5, xend= 10589), colour = "red", size = 2.5)+
  theme(panel.grid.major= element_blank(), text = element_text(family = "serif"), axis.line=element_blank(), axis.ticks.y = element_blank(), panel.background = element_blank(), axis.title.y =  element_blank(), axis.title.x = element_blank())+
  guides(fill=FALSE)

paln


geom_text(aes(label = raw), size = 8, family="serif") +
  
  ########

seqfs <- adegenet::fasta2DNAbin("./Mito_aln_cladeA.fas")
coverageest<-data.frame("SampleID" = rownames(seqfs), "Length" = rep(13438, 29), "Missing" = rep(NA, 29))

for(i in 1:nrow(seqfs)){
  foo <- function(x) sum(x == 240)   #same as x[x==240], i.e. number of N's in alignment
  g <- apply(seqfs[i,], 2, foo)   #apply the function foo to all the columns (2) in alignment (x); g <- vector of number of N's in each column (length = number of columns) 
  coverageest$Missing[i] <- length(which(g==1)) 
}

coverageest$Perc.Coverage <- ((coverageest$Length - coverageest$Missing) / coverageest$Length) * 100
sjPlot::tab_df(coverageest, file="sj_Coverage.doc", alternate.rows = T, describe = FALSE)

SNPseq<-adegenet::fasta2DNAbin("Mito_aln_cladeA.fas",snpOnly=TRUE)
dimseq<-(dim(SNPseq)) 
paste("Total SNPs", dimseq[2], sep="=")

ape::checkAlignment(seqfs, what=1)


##Gene plot

genecoords<-read.csv("./gene_coords_plot2.csv") #50 obs

genecoords<-genecoords[!is.na(genecoords$GeneID),]  #remove all non-coding regions (36)

genecoords$GeneID2<-as.character(genecoords$GeneID)  #create new simplified geneID classification
genecoords[genecoords$GeneType == "NPCG" & genecoords$GeneID!= "rrnS" & genecoords$GeneID!= "rrnL", "GeneID2"]<-as.character(seq(1,22,1))
genecoords$GeneID2<-as.factor(genecoords$GeneID2)


genecoords$GeneType2<-as.character(genecoords$GeneType)
#genecoords[genecoords$GeneType == "PCG", "GeneType2"]<-c("COs", "NADs", "NADs", "NADs", "NADs", "CYTB", "ATP6", "COs", "NADs", "NADs", "COs", "NADs") #group PCGs by domain
genecoords[genecoords$GeneType == "NPCG" & genecoords$GeneID!= "rrnS" & genecoords$GeneID!= "rrnL", "GeneType2"]<-as.character(rep("tRNA", 22)) #group tRNAs
genecoords[genecoords$GeneType2 == "NPCG", "GeneType2" ]<-c("rRNA", "rRNA") #rename last two NPCGs to rRNAs
genecoords$GeneType2<- as.factor(genecoords$GeneType2)

genecoords$Length <- genecoords$Stop - genecoords$Start + 1

#gencolpal<- c("#3cb44b", "#01665e", "#4363d8", "#469990", "#8c510a","#dfc27d")
gencolpal<- c("#469990", "#8c510a","#dfc27d")
show_col(gencolpal)


geneplot<-
  ggplot(genecoords, aes(xmin = Start, xmax = Stop,y = Name, fill = GeneType2, label = GeneID2)) +
  geom_gene_arrow(arrowhead_height = unit(6,"mm"), arrowhead_width = unit(5, "mm"), arrow_body_height = unit(5, "mm")) +
  geom_gene_label(align =  "left", padding.y = grid::unit(0.01, "lines"), min.size = 4) +
  scale_fill_manual(values = gencolpal) +
  labs(y = "")+
  theme(axis.text.y = element_blank() , text = element_text(family = "serif"), axis.ticks.y = element_blank(), panel.background = element_blank(), axis.title.y =  element_blank(), axis.line = element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0,0,0,0), "mm") )+
  guides(fill=FALSE)


geneplot

##NAD4L, NAD6


### Circlize plot
library(circlize)
library(ape)
library(seqinr)
library(pegas)
library(seqRFLP)
library(plotrix)
library(adegenet)
library(Biostrings)
library(dplyr)

genecoordsDF <- read.csv("gene_coords_plot2.csv", stringsAsFactors = F)
genecoordsDF <- na.exclude(genecoordsDF)

### Change gene colours to be colourblind-friendly
genecoordsDF$itemRGB[genecoordsDF$itemRGB == "#f58231"] <- "#dfc27d" #tRNAs
genecoordsDF$itemRGB[genecoordsDF$itemRGB == "#e6194B"] <- "#8c510a" #rRNAs
genecoordsDF[genecoordsDF$GeneType == "PCG", "itemRGB"] <- "#469990" #PCGs

coloursgen<-genecoordsDF$itemRGB

circos.clear()
op <- par(family = "serif")
circos.par(cell.padding = c(0, 0.02, 0, 0.02))
#circos.par(track.margin=c(0,0))
circos.genomicInitialize(genecoordsDF,
                         plotType = "axis", major.by = 1000,
                         tickLabelsStartFromZero = FALSE, axis.labels.cex = 1)


circos.info()


#### Gene ID track 
circos.track(
  ylim = c(0, 1),
  bg.col = NA,
  bg.border = NA,
  track.height = 0.1,
  track.margin = c(0, 0),
  panel.fun = function(region, value) {
    circos.genomicRect(
      region = genecoordsDF[, 2:3],
      value = genecoordsDF,
      col = genecoordsDF$itemRGB,
      border = "black"
    )
  }
)

##### Add labels
genecoordsGENES <- genecoordsDF[genecoordsDF$itemRGB != "#dfc27d", ]

circos.genomicText(
  genecoordsGENES[, 2:3],
  value = genecoordsGENES,
  labels.column = "GeneID",
  y = 0.5,
  facing = "inside",
  cex = 1,
  track.index = 2
)

genecoordsTRNA <- genecoordsDF[genecoordsDF$itemRGB == "#dfc27d", ]
genecoordsTRNA$GeneID <- gsub("tRNA_", "", genecoordsTRNA$GeneID)


circos.genomicText(
  genecoordsTRNA[, 2:3],
  value = genecoordsTRNA,
  labels.column = "GeneID",
  y = 0.5,
  facing = "clockwise",
  cex = 0.6,
  track.index = 2
)


##### Adding lines for amplicons (32) - track 1
xy<-1
amplicons <- read.csv("amplicon_coords_whole.csv", header = TRUE, stringsAsFactors = F)

ampID1<-seq(1,27,1)
ampID2<-seq(1,27,1)
ampID<-c(ampID1,ampID2)
ampID<-sort(ampID)
ampID<-c(27,27,ampID)
amplicons$ampID<-ampID

amplicons$yaxis_mod[amplicons$yaxis_mod==1]<-0
amplicons$yaxis_mod[amplicons$yaxis_mod==2]<-0.3
amplicons$yaxis_mod[amplicons$yaxis_mod==3]<-0.6
amplicons$yaxis_mod[amplicons$yaxis_mod==4]<-0.9

amplicons$yaxis_mod[amplicons$GeneID == "Dv13F1-Dv0R4"] <-0
# circos.track(
#   ylim = c(0, 1),
#   bg.col = NA,
#   bg.border = NA,
#   track.height = 0.1,
#   track.margin = c(0, 0),
#   panel.fun = function(x, y) {
#     for (xy in 1:length(amplicons$xaxis_mod)) {
#       if ((amplicons$GeneID[xy] == amplicons$GeneID[xy + 1]) == TRUE) {
#         circos.lines(
#           x = c(amplicons$xaxis_mod[xy], amplicons$xaxis_mod[xy + 1]),
#           y = c(amplicons$yaxis_mod[xy], amplicons$yaxis_mod[xy + 1]),
#           type = "l",
#           col = "#543005",
#           lwd = 4,
#           track.index = 3,
#           baseline = "bottom"
#         )
#         
#       } else{
#         next
#       }
#     }
#   }
# )

amplicons$GeneID <- gsub( "Dvmt8854", "Dvmt8127", amplicons$GeneID)


xy<-1
for (xy in 1:length(amplicons$xaxis_mod)) {
  if ((amplicons$GeneID[xy] == amplicons$GeneID[xy + 1]) == TRUE) {
    amplicons$Start[xy]<-amplicons$xaxis_mod[xy]
    amplicons$Stop[xy]<-amplicons$xaxis_mod[xy + 1]
    amplicons$Start[xy+1]<-amplicons$xaxis_mod[xy]
    amplicons$Stop[xy+1]<-amplicons$xaxis_mod[xy + 1]
  } else{
    next
  }
}


ampliconsmod<-amplicons %>% distinct(amplicons$GeneID, .keep_all = TRUE)
ampliconsmod$ampID<-rownames(ampliconsmod)
ampliconsmod[ampliconsmod$ampID == 28, "ampID"] <- 1

amplicons1<-ampliconsmod[ampliconsmod$yaxis_mod==0,]
amplicons2<-ampliconsmod[ampliconsmod$yaxis_mod==0.3,]
amplicons3<-ampliconsmod[ampliconsmod$yaxis_mod==0.6,]
amplicons4<-ampliconsmod[ampliconsmod$yaxis_mod==0.9,]

circos.track(
  ylim = c(0, 1),
  bg.col = NA,
  bg.border = NA,
  track.height = 0.03,
  track.margin = c(0, 0),
  panel.fun = function(region, value) {
    circos.genomicRect(
      region = amplicons1[,5:6],
      value = amplicons1,
      col = "#543005",
      border = "black"
    )
  }
)

circos.info()

circos.genomicText(
  amplicons1[, 5:6],
  value = amplicons1,
  labels.column = "ampID",
  y = 0.5,
  facing = "inside",
  cex = 0.8,
  track.index = 3,
  col="white"
)

circos.track(
  ylim = c(0, 1),
  bg.col = NA,
  bg.border = NA,
  track.height = 0.03,
  track.margin = c(0, 0),
  panel.fun = function(region, value) {
    circos.genomicRect(
      region = amplicons2[,5:6],
      value = amplicons2,
      col = "#543005",
      border = "black"
    )
  }
)

circos.genomicText(
  amplicons2[, 5:6],
  value = amplicons2,
  labels.column = "ampID",
  y = 0.5,
  facing = "inside",
  cex = 0.8,
  track.index = 4,
  col="white"
)

circos.track(
  ylim = c(0, 1),
  bg.col = NA,
  bg.border = NA,
  track.height = 0.03,
  track.margin = c(0, 0),
  panel.fun = function(region, value) {
    circos.genomicRect(
      region = amplicons3[,5:6],
      value = amplicons3,
      col = "#543005",
      border = "black"
    )
  }
)

circos.genomicText(
  amplicons3[, 5:6],
  value = amplicons3,
  labels.column = "ampID",
  y = 0.5,
  facing = "inside",
  cex = 0.8,
  track.index = 5,
  col="white"
)

circos.track(
  ylim = c(0, 1),
  bg.col = NA,
  bg.border = NA,
  track.height = 0.03,
  track.margin = c(0, 0),
  panel.fun = function(region, value) {
    circos.genomicRect(
      region = amplicons4[,5:6],
      value = amplicons4,
      col = "#543005",
      border = "black"
    )
  }
)

circos.genomicText(
  amplicons4[, 5:6],
  value = amplicons4,
  labels.column = "ampID",
  y = 0.5,
  facing = "inside",
  cex = 0.8,
  track.index = 6,
  col="white"
)



## Gene Diversity
setwd("./PCGenes")
genes<-("../gene_coords_plot2.csv")
GenePart <- function (alignment, genecoords, genetype, gene)
{
  fastaFile <- Biostrings::readDNAStringSet(alignment)
  metadata <- read.csv(genecoords, header = T)
  CGpartitions <- metadata[metadata$GeneType == genetype, ]
  Fastasubseq <-
    subseq(fastaFile, start = CGpartitions$Start[gene], end = CGpartitions$Stop[gene])
  outfilen <-
    file.path(getwd(), file.name = (paste(CGpartitions$GeneID[gene], "fas", sep =
                                            ".")))
  seqRFLP::dataframe2fas(Fastasubseq, file = outfilen)
}

## to separate out into individual fastas
for (k in 1:nrow(genecoords[genecoords$GeneType == "PCG",])) {
  GenePart(
    alignment = "../Mito_aln_cladeA.fas",
    genecoords = genes,
    genetype = "PCG",
    gene = k
  )
}


genefiles<-list.files(pattern=".fas")

for(i in 1:length(genefiles)){
  codon_del<-adegenet::fasta2DNAbin(genefiles[i])
  codon_del<-del.codongaps(codon_del, threshold = 0.11)
  ape::write.dna(codon_del, file= paste(genefiles[i]), format = "fasta", colw=9999)
}


for(i in 1:length(genefiles)){
  coldel_diversityStats(genefiles, threshold = 0.11)
}

gendiv<-read.csv("Coldel_DivStats_0.11.csv")
gendiv$GeneID<-genecoords$GeneID[1:11]
sjPlot::tab_df(gendiv, file="sj_PCGDiv.doc", alternate.rows = T, describe = FALSE)



### tRNA Diversity

setwd("../NPCGenes")
genes<-("../gene_coords_plot2.csv")
genecoords<-read.csv("../gene_coords_plot2.csv")
GenePart <- function (alignment, genecoords, genetype, gene)
{
  fastaFile <- Biostrings::readDNAStringSet(alignment)
  metadata <- read.csv(genecoords, header = T)
  CGpartitions <- metadata[metadata$GeneType == genetype, ]
  Fastasubseq <-
    subseq(fastaFile, start = CGpartitions$Start[gene], end = CGpartitions$Stop[gene])
  outfilen <-
    file.path(getwd(), file.name = (paste(CGpartitions$GeneID[gene], "fas", sep =
                                            ".")))
  seqRFLP::dataframe2fas(Fastasubseq, file = outfilen)
}

## to separate out into individual fastas
for (k in 1:nrow(genecoords[genecoords$GeneType == "NPCG",])) {
  GenePart(
    alignment = "../Mito_aln_cladeA.fas",
    genecoords = genes,
    genetype = "NPCG",
    gene = k
  )
}

genefiles<-list.files(pattern=".fas")

for(i in 1:length(genefiles)){
  coldel_diversityStats(genefiles, threshold = 0.11)
}

NPCGs<-genecoords[genecoords$GeneType == "NPCG",]
gendiv<-read.csv("Coldel_DivStats_0.11.csv")
gendiv$GeneID<-NPCGs$GeneID
sjPlot::tab_df(gendiv, file="sj_NPCGDiv.doc", alternate.rows = T, describe = FALSE)


###### Sliding Window
# Sliding Window Final Analysis

source("./DivStats.R")

## Contig boundaries

starts<-c(231, 2516, 3714, 4154, 5488, 6905, 9162, 10030, 12710)
stops<-c(2436, 3619, 4034, 5210, 6418, 8642, 9987, 11717, 13537)
lengths<-stops-starts+1
contigL_coord<-data.frame(Start = starts, Stop = stops, Length = lengths)
contigL500<-contigL_coord[contigL_coord$Length > 499,]


## Sliding Windows
for (k in 1:length(stops)) {
  SlidPart(
    alignment = "Mito_aln.fas",
    coords = contigL_coord,
    ID = k
  )
}


## Run diversity stats 
fastas<-list.files(pattern = "Win_")
coldel_diversityStats(fastas, threshold = 0.11) 


## Missing data outside of boundaries

missLongStart<-c(1,stops[1:8],13011) #add in NA for window data that is removed at the end (last window is 13010-13510 with a single point at 13010)
missLongStop<-c(230,starts[2:9],13438) #add in NA for window data that is removed at the end

m1<-seq(missLongStart[1],missLongStop[1],50)
m2<-seq(missLongStart[2],missLongStop[2],50)
m3<-seq(missLongStart[3],missLongStop[3],50)
m4<-seq(missLongStart[4],missLongStop[4],50)
m5<-seq(missLongStart[5],missLongStop[5],50)
m6<-seq(missLongStart[6],missLongStop[6],50)
m7<-seq(missLongStart[7],missLongStop[7],50)
m8<-seq(missLongStart[8],missLongStop[8],50)
m9<-seq(missLongStart[9],missLongStop[9],50)
m10<-seq(missLongStart[10],missLongStop[10],50)

mall<-c(m1, m2, m3, m4, m5, m6, m7, m8, m9,m10, 13438)

missingcontigsLONG<-data.frame(fastafiles=rep(NA,length(mall)), Len=rep(NA,length(mall)), N=rep(NA,length(mall)), Nuc.Div=rep(NA,length(mall)), GC=rep(NA,length(mall)), PIS=rep(NA,length(mall)), Num.Hap=rep(NA,length(mall)), SNVs=rep(NA,length(mall)), SNPs=rep(NA,length(mall)), Singletons=rep(NA,length(mall)), mtDNA=mall)


## Plots
contigL500_start1<-seq(contigL500$Start[1],contigL500$Stop[1],50)
contigL500_stop1<-seq((contigL500$Start[1]+499),contigL500$Stop[1],50)
contigL500_start1<-contigL500_start1[1:length(contigL500_stop1)]

contigL500_start2<-seq(contigL500$Start[2],contigL500$Stop[2],50)
contigL500_stop2<-seq((contigL500$Start[2]+499),contigL500$Stop[2],50)
contigL500_start2<-contigL500_start2[1:length(contigL500_stop2)]

contigL500_start3<-seq(contigL500$Start[3],contigL500$Stop[3],50)
contigL500_stop3<-seq((contigL500$Start[3]+499),contigL500$Stop[3],50)
contigL500_start3<-contigL500_start3[1:length(contigL500_stop3)]

contigL500_start4<-seq(contigL500$Start[4],contigL500$Stop[4],50)
contigL500_stop4<-seq((contigL500$Start[4]+499),contigL500$Stop[4],50)
contigL500_start4<-contigL500_start4[1:length(contigL500_stop4)]

contigL500_start5<-seq(contigL500$Start[5],contigL500$Stop[5],50)
contigL500_stop5<-seq((contigL500$Start[5]+499),contigL500$Stop[5],50)
contigL500_start5<-contigL500_start5[1:length(contigL500_stop5)]

contigL500_start6<-seq(contigL500$Start[6],contigL500$Stop[6],50)
contigL500_stop6<-seq((contigL500$Start[6]+499),contigL500$Stop[6],50)
contigL500_start6<-contigL500_start6[1:length(contigL500_stop6)]

contigL500_start7<-seq(contigL500$Start[7],contigL500$Stop[7],50)
contigL500_stop7<-seq((contigL500$Start[7]+499),contigL500$Stop[7],50)
contigL500_start7<-contigL500_start7[1:length(contigL500_stop7)]

contigL500_start8<-seq(contigL500$Start[8],contigL500$Stop[8],50)
contigL500_stop8<-seq((contigL500$Start[8]+499),contigL500$Stop[8],50)
contigL500_start8<-contigL500_start8[1:length(contigL500_stop8)]

contigL500_Starts<-c(contigL500_start1, contigL500_start2, contigL500_start3, contigL500_start4, contigL500_start5, contigL500_start6, contigL500_start7, contigL500_start8)  
contigL500_Stops<-c(contigL500_stop1, contigL500_stop2, contigL500_stop3, contigL500_stop4, contigL500_stop5, contigL500_stop6, contigL500_stop7, contigL500_stop8)

contigL500_windowcoords<-data.frame(Start = contigL500_Starts, Stop = contigL500_Stops)

########

### strsplit wins in each 
strsplitWinID<- function(x){
  y<-base::strsplit(as.character(x), "_")
  x<-sapply(y, "[[",2)
  x<-gsub(".fas", "", x)
  x<-as.numeric(x)  
}

setwd("../Sliding_Window_500x50")
cl500<-read.csv("Coldel_DivStats_0.11.csv")
#cl500$fastafiles<-strsplitWinID(cl500$fastafiles)

cl500<-cl500[order(cl500$fastafiles),] 

cl500$mtDNA<-contigL500_windowcoords$Start


### rbind with missing windows created

cl500<-rbind(cl500,missingcontigsLONG)

cl500$GC<-cl500$GC*100




cl500$Start<-cl500$mtDNA
cl500[!is.na(cl500$fastafiles), "Start"] <- contigL500_windowcoords$Start
cl500$Stop<-cl500$Start
cl500[!is.na(cl500$fastafiles), "Stop"] <- contigL500_windowcoords$Stop
MissingStops <- cl500[is.na(cl500$fastafiles), "Stop"] +499
cl500[is.na(cl500$fastafiles), "Stop"]<- MissingStops

cl500$mtDNA<-ceiling((cl500$Start+cl500$Stop)/2)
cl500[is.na(cl500$fastafiles), "mtDNA"]<- MissingStops-500
cl500$binwidth<-500


## make sure they are all as.numeric()

cols.num<-c("PIS", "Num.Hap", "SNVs", "SNPs", "Singletons")

cl500[cols.num]<-sapply(cl500[cols.num], as.numeric)

str(cl500)

min(cl500$Nuc.Div, na.rm = T) #0
max(cl500$Nuc.Div, na.rm = T) #0.007

min(cl500$Num.Hap, na.rm = T) #1
max(cl500$Num.Hap, na.rm = T) #6

min(cl500$SNVs, na.rm = T) #0
max(cl500$SNVs, na.rm = T) #14


### Add partial-COI details. 
setwd("../TUN_COI")
amplicon_coords<-read.csv("../Amplicon_coords_whole.csv")
COI_coords<-amplicon_coords[amplicon_coords$GeneID == "TUN",]

COI_coords<-data.frame("GeneID" = "COI", "Start" = 10004, "Stop" = 10589) #10030 is where the alignment starts from my own sequences

fastaFile <- Biostrings::readDNAStringSet("../Mito_aln_cladeA.fas")

Fastasubseq <-subseq(fastaFile, start = COI_coords$Start[1], end = COI_coords$Stop[1])
outfilen <- file.path(getwd(), file.name = "TUN_COI.fas")
seqRFLP::dataframe2fas(Fastasubseq, file = outfilen)

TUNfas<- list.files(pattern=".fas")
coldel_diversityStats(TUNfas, threshold = 0.11)

TUNdiv<-read.csv("Coldel_DivStats_0.11.csv")
TUNdiv$GC<-TUNdiv$GC*100

# mySequences <- Biostrings::readAAStringSet(TUNfas)
# TUNAln<-msa(mySequences)
# print(TUNAln, showConsensus = FALSE)
# msaPrettyPrint(TUNAln, output="pdf",
#                showNames="none", showLogo="none", file = "Aln.pdf", askForOverwrite=FALSE)
# library(msa)
# 
# tools::texi2pdf("TUNAln.tex")

### Add concatenated Mito details

Mitofas<-"../Mito_aln_cladeA.fas"

setwd("../Whole_Mito")
coldel_diversityStats(Mitofas, threshold = 0.11)

Mitodiv<-read.csv("Coldel_DivStats_0.11.csv")
Mitodiv$GC<-Mitodiv$GC*100

sjPlot::tab_df(Mitodiv, file="sj_Mitodiv.doc", alternate.rows = T, describe = FALSE)

### Plots
setwd("../Sliding_Window_500x50")

pnuc<-ggplot(data=cl500, aes(x=Start, y=Nuc.Div))+
  geom_line()+
  labs(y= bquote(""~ pi ~""), x="")+
  scale_x_continuous(breaks = c(1, seq(1000, 13528, by = 1000)))+
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.002))+
  geom_line(aes(y= TUNdiv$Nuc.Div), color = "red", linetype = "dashed")+
  geom_text(aes(color = "red", x= 10, y = TUNdiv$Nuc.Div, label = paste(signif(TUNdiv$Nuc.Div, digits=2))), vjust = -0.5, show.legend = FALSE, family = "serif")+
  theme_tufte()+
  theme(text = element_text(size=14, family = "serif"), axis.line= element_line(size = 0.5, color = "black"))



pSNV<-ggplot(data=cl500, aes(x=Start, y=SNPs))+
  geom_line()+
  labs(y= "VS", x="")+
  scale_x_continuous(breaks = c(1, seq(1000, 13528, by = 1000)))+
  scale_y_continuous(breaks = seq(0, 10, by = 2))+
  geom_line(aes(y= TUNdiv$SNVs), color = "red", linetype = "dashed")+
  geom_text(aes(color = "red", x= 10, y = TUNdiv$SNVs, label = paste(signif(TUNdiv$SNVs, digits=2))), vjust = -0.5, show.legend = FALSE, , family = "serif")+
  theme_tufte()+
  theme(text = element_text(size=14, family = "serif"), axis.line= element_line(size = 0.5, color = "black"))



pNH<-ggplot(data=cl500, aes(x=Start, y=Num.Hap))+
  geom_line(color = "black")+
  labs(y=expression(paste(H[N]) ), x="mtDNA")+
  scale_x_continuous(breaks = c(1, seq(1000, 13528, by = 1000)))+
  scale_y_continuous(breaks = seq(1, 6, by = 2))+
  geom_line(aes(y= TUNdiv$Num.Hap), color = "red", linetype = "dashed")+
  geom_text(aes( x= 10, y = TUNdiv$Num.Hap, label = paste(signif(TUNdiv$Num.Hap, digits=2))), color = "red", vjust = -0.5, show.legend = FALSE, family = "serif")+
  theme_tufte()+
  theme(text = element_text(size=14, family = "serif"), axis.line= element_line(size = 0.5, color = "black"))


plotlist<-list(geneplot, paln, pNH , pSNV, pnuc)
cowplot::plot_grid(plotlist=plotlist, labels = c("", "", "A", "B", "C"), label_size=12, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.5,3,2,2,2,2))

##without GC content
plotlist<-list(geneplot, paln, pNH , pSNV,pnuc)


cowplot::plot_grid(plotlist=plotlist, labels = c("", "", "A", "B", "C"), label_size=12, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.5,4,2,2,2))




### Candidate amplicons

### Candidate amplicons
setwd("../")
candidate_coord<- data.frame(candidate = c("NAD1", "COI"), Start = c(1469, 9314), Stop = c(2444, 11280))
rownames(candidate_coord)
str(candidate_coord)
source("DivStats.R")


k<-1
for (k in 1:nrow(candidate_coord)) {
  SlidPart(
    alignment = "Mito_aln_cladeA.fas",
    genecoords = candidate_coord,
    gene = k,
    prefix = "Candidate"
  )
}


## Run diversity stats 
fastas<-list.files(pattern = "Candidate")
coldel_diversityStats(fastas, threshold = 0.11) 

