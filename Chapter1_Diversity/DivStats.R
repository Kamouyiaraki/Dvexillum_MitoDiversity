## Not run (https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/#:~:text=Install%20%7C%20install%20%26%20load%20packages&text=check()%20function%20basically%20goes,Otherwise%2C%20load%20it)

# packages = c("adegenet", "ape",
#              "Biostrings", "dplyr", "haplotypes", "pegas", "phylotools", 
#              "plotrix", "plyr", "seqinr", "seqRFLP", "sidier")
# 
# ## Now load or install&load all
# lapply(
#   packages,
#   FUN = function(x) {
#     if (!require(x, character.only = TRUE)) {
#       install.packages(x, dependencies = TRUE)
#       library(x, character.only = TRUE)
#     }
#   }
# )


##Library
library(adegenet)
library(ape)
library(Biostrings)
library(dplyr)
library(haplotypes)
library(pegas)
library(phylotools)
library(plotrix)
library(plyr)
library(seqinr)
library(seqRFLP)
library(sidier)
#library(phyloch)  ##not available for R v4.0 as of March 2020


# as.integer(ape::as.DNAbin("-")) #4
# as.integer(ape::as.DNAbin("a")) #136
# as.integer(ape::as.DNAbin("g")) #72
# as.integer(ape::as.DNAbin("c")) #40
# as.integer(ape::as.DNAbin("t")) #24
# as.integer(ape::as.DNAbin("n")) #240

##Sliding window - individual fastas 

SlidPart <- function (alignment, genecoords, gene, prefix = "Win_")
{
  fastaFile <- Biostrings::readDNAStringSet(alignment)
  CGpartitions <- genecoords
  Fastasubseq <-subseq(fastaFile, start = CGpartitions$Start[gene], end = CGpartitions$Stop[gene])
  outfilen <- file.path(getwd(), file.name = (paste(prefix, rownames(CGpartitions)[gene], ".fas", sep="")))
  seqRFLP::dataframe2fas(Fastasubseq, file = outfilen)
} 


## String split function for changing fasta file ID to as.numeric() Window number: 
strsplitWinID<- function(x){
  y<-base::strsplit(as.character(x), "_")
  x<-sapply(y, "[[",2)
  x<-gsub(".fas", "", x)
  x<-as.numeric(x)  
}

## Delete column missing data in DNAbin
### Edited version of ape::del.colgapsonly (changed to remove columns with 'N' (x==240) instead of '-' (x==4), and changed `which(g / nrow(x) >= threshold)` to `which(g / nrow(x) > threshold)`)
del.colgapsonly <- function(x, threshold = 1, freq.only = FALSE)
{
  if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)    #object-oriented method based on class of argument
  if (!is.matrix(x)) stop("DNA sequences not in a matrix")   #make sure DNAbin is in matrix format
  foo <- function(x) sum(x == 240)   #same as x[x==240], i.e. number of N's in alignment
  g <- apply(x, 2, foo)   #apply the function foo to all the columns (2) in alignment (x); g <- vector of number of N's in each column (length = number of columns) 
  if (freq.only) return(g)
  i <- which(g / nrow(x) > threshold) # g/nrow(x) = the proportion of samples (rows) in alignment (x) that are n's (g); `which... > threshold` specifies columns that are above the threshold (those where g/nrow(x)>threshold == TRUE)
  if (length(i)) x <- x[, -i] #if there are any columns above the threshold (i.e. length(i) is a valid function), remove those columns (i) from alignment (x)
  x #updated alignment (x)
}

### Edited version of ape::del.rowgapsonly(changed to remove rows with 'N' (x == 240) instead of '-' (x == 4), removed threshold argument and changed `which(g / nrow(x) >= threshold)` to `which(g / nrow(x) == 1)`)
del.rowgapsfull <- function(x, freq.only=FALSE)
{
  if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)
  if (!is.matrix(x)) stop("DNA sequences not in a matrix")
  foo <- function(x) sum(x == 240)
  g <- apply(x, 1, foo) #apply function to all rows (1) in alignment (x)
  if (freq.only) return(g)
  i <- which(g / ncol(x) == 1) #NOTE. this is a threshold that can be changed to accomodate nucleotide diversity estimates
  if (length(i)) x <- x[-i, ]
  x
}

### Edited version of ape::del.rowgapsonly(changed to remove rows with 'N' (x == 240) instead of '-' (x == 4), removed threshold argument and added in conditional loops for identifying columns to delete in codon patter. This assumes that the start of DNAbin is codon 1)
del.rowgapsonly <- function(x, freq.only=FALSE, threshold=1)
{
  if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)
  if (!is.matrix(x)) stop("DNA sequences not in a matrix")
  foo <- function(x) sum(x == 240)
  g <- apply(x, 1, foo) #apply function to all rows (1) in alignment (x)
  if (freq.only) return(g)
  i <- which(g  >= threshold) 
  if (length(i)) x <- x[-i, ]
  x
}

### Edited version of ape::del.rowgapsonly(changed to remove rows with 'N' (x == 240) instead of '-' (x == 4), removed threshold argument and added in conditional loops for identifying columns to delete in codon patter. This assumes that the start of DNAbin is codon 1)
del.codongaps <- function(x, threshold) {
  
  if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)
  if (!is.matrix(x))stop("DNA sequences not in a matrix") 
  
  foofunct <- function(y) sum(y ==  240)
  m <- apply(x, 2, foofunct)
  removed <- which(m / nrow(x) > threshold)
  coltoremove <- as.vector(removed)   #24 with 3 missing
  
  if (length(removed)){
  for (i in 1:length(removed)) {
    if (removed[i] %% 3 == 0) {
      thirdx <- removed[removed %% 3 == 0]
      funthird <- function(x) c(x - 1, x - 2)
      coltoremove <- append(coltoremove, funthird(thirdx))
      
    } else if ((removed[i] + 1) %% 3 == 0) {
      secondx <- removed[(removed + 1) %% 3 == 0]
      funsecond <-  function(x) c(x - 1, x + 1)
      coltoremove <- append(coltoremove, funsecond(secondx))
      
    } else {
      firstx <- removed[i]
      funfirst <- function(x) c(x + 1, x + 2)
      coltoremove <- append(coltoremove, funfirst(firstx))
    }
  }
  coltoremove <- sort(coltoremove)
  coltoremove <- unique(coltoremove)
  
  if (length(coltoremove)) x <- x[, -coltoremove]
  x
  } else{
  print("There's nothing to remove here")
  }
  
}

##modified from ape::del.colgapsonly()


coldel_diversityStats <- function(fastafiles, threshold, nprefix="Coldel_DivStats_") {
  WindowStats <- data.frame(fastafiles)
  
  pb = txtProgressBar(min = 0,
                      max = length(fastafiles),
                      initial = 0)
  
  
  for (i in 1:length(fastafiles)) {
    seqfs <- adegenet::fasta2DNAbin(fastafiles[i])
    
    Align1<-del.colgapsonly(seqfs, threshold = threshold, freq.only = FALSE)
    
    if ((ncol(Align1) > 0)== TRUE){
      removedls<-list()
      WindowID<-fastafiles[i]
      foofunct<-function(y) sum(y ==  240) 
      m <-apply(Align1, 1, foofunct)
      removed<-which(m/ncol(Align1) == 1) #NOTE. this is a threshold that can be changed to accomodate nucleotide diversity estimates
      if (length(removed)) removedls[[i]]<-list("WindowID" = WindowID, "Removed_Samples" = removed ) 
      
      Align<-del.rowgapsfull(Align1, freq.only = FALSE) 
      
      setTxtProgressBar(pb, i)
      message('File ', i, ' of ', length(fastafiles))
      
      WindowStats$Len[i]<-ncol(Align)
      
      WindowStats$Nuc.Div[i] <- pegas::nuc.div(Align, pairwise.deletion=TRUE)
      
      WindowStats$GC[i]<- ape::GC.content(Align)
      
      #WindowStats$PIS[i] <- phyloch::pis(Align)   #number of PIS
      
      hap<-haplotypes::as.dna(Align)
      WindowStats$Num.Hap[i]<-length(haplotypes::haplotype(hap)) #includes haplotype based on indel ("-")
      #hap <- pegas::haplotype(Align)
      #WindowStats$Num.Hap[i] <- nrow(hap)
      
      WindowStats$SNVs[i] <- if(length(adegenet::DNAbin2genind(Align, exp.char = c("a", "g", "c", "t", "-")))==0) as.numeric(0) else length(adegenet::locNames(adegenet::DNAbin2genind(Align, exp.char = c("a", "g", "c", "t", "-"))))
      
      WindowStats$SNPs[i] <- if(length(adegenet::DNAbin2genind(Align, exp.char = c("a", "g", "c", "t")))==0) as.numeric(0) else length(adegenet::locNames(adegenet::DNAbin2genind(Align, exp.char = c("a", "g", "c", "t"))))
      
      #WindowStats$Singletons <- WindowStats$SNPs - WindowStats$PIS
      
      
    } else{
      WindowStats$Len[i]<-ncol(Align1)  ##keep the length of window for now
      WindowStats$Nuc.Div[i] <- NA
      WindowStats$GC[i]<-NA
     # WindowStats$PIS[i] <- NA  
      WindowStats$Num.Hap[i] <- NA
      WindowStats$SNVs[i] <- NA
      WindowStats$SNPs[i] <- NA
     # WindowStats$Singletons <- NA
    }  
  }
  
  outfilen<-paste(nprefix, threshold, ".csv", sep="")  
  write.csv(WindowStats, outfilen, row.names = FALSE, append = TRUE)
  print(WindowStats)
  
  outfilelist<-file("RemovedSamples.txt")
  sink(outfilelist, append=TRUE)
  print(removedls)
  sink()
  close(outfilelist)
}


## Window by Window
diversityStats <- function(fastafiles, threshold) {
  WindowStats <- data.frame(fastafiles)
  
  pb = txtProgressBar(min = 0,
                      max = length(fastafiles),
                      initial = 0)
  
  
  for (i in 1:length(fastafiles)) {
    Align <- adegenet::fasta2DNAbin(fastafiles[i])
    
    setTxtProgressBar(pb, i)
    message('File ', i, ' of ', length(fastafiles))
    
    if ((length(as.integer(Align[Align == 240]))) / (length(Align)) <=
        threshold) {
      WindowStats$Nuc.Div[i] <- pegas::nuc.div(Align, pairwise.deletion=TRUE)
      
      WindowStats$GC[i]<- ape::GC.content(Align)
      
      #WindowStats$PIS[i] <- phyloch::pis(Align)   #number of PIS
      
      hap<-haplotypes::as.dna(Align)
      WindowStats$Num.Hap[i]<-length(haplotypes::haplotype(hap)) #includes haplotype based on indel ("-")
      #hap <- pegas::haplotype(Align)
      #WindowStats$Num.Hap[i] <- nrow(hap)
      
      WindowStats$SNVs[i] <- if(length(adegenet::DNAbin2genind(Align, exp.char = c("a", "g", "c", "t", "-")))==0) 0 else length(adegenet::locNames(adegenet::DNAbin2genind(Align, exp.char = c("a", "g", "c", "t", "-"))))
      
      WindowStats$SNPs[i] <- if(length(adegenet::DNAbin2genind(Align, exp.char = c("a", "g", "c", "t")))==0) 0 else length(adegenet::locNames(adegenet::DNAbin2genind(Align, exp.char = c("a", "g", "c", "t"))))
      
      #WindowStats$Singletons[i] <- WindowStats$SNPs[i] - WindowStats$PIS[i]
      
      
    } else{
      WindowStats$Nuc.Div[i] <- NA
      WindowStats$GC[i]<-NA
      #WindowStats$PIS[i] <- NA   #number of PIS
      WindowStats$Num.Hap[i] <- NA
      WindowStats$SNVs[i] <- NA
      WindowStats$SNPs[i] <- NA
     # WindowStats$Singletons <- NA
    }
  }
  outfilen<-paste("DivStats", threshold, ".csv")  
  write.csv(WindowStats, outfilen, row.names = FALSE, append = TRUE)
  print(WindowStats)
}

##Split FASTA files according to row
FastabyInd<-function(aln, outputdir)
{
  fastaFile<- Biostrings::readDNAStringSet(aln)
  for (i in 1:length(fastaFile)){
    outfilen<-file.path(outputdir, file.name=paste(names(fastaFile[i]),aln, sep="_"))
    seqRFLP::dataframe2fas(fastaFile[i], file=outfilen)
  }
}

##Output FASTA file on unique sequences (option to set threshold for number of nucleotides that are N (>=) for sequence to be removed = nuc.threshold)
get.unique <-  function(inputFAS=NA, outfile="UniqueSeq.fas", nuc.threshold=1) {
  align<-read.dna(file=inputFAS, format="fasta")
  
  foo <- function(x) sum(x == 240)
  g <- apply(align, 1, foo) #apply function to all rows (1) in alignment (x)
  i <- which(g  >= nuc.threshold) 
  
  haploseq<-sidier::FindHaplo(align=align,saveFile=FALSE)
  uniqueseq<-dplyr::distinct(as.data.frame(haploseq), Haplotype.Name, .keep_all = TRUE)
  u.rm<-which(!(rownames(align) %in% uniqueseq[,1]))
  if(length(i))u.rm<-c(u.rm,i)
  u.rm<-rownames(align)[u.rm]
  
  phylotools::rm.sequence.fasta(inputFAS, outfile = outfile, to.rm=u.rm)
}


##Concatenate fasta files based on sequence ID, note uses a vector of fasta file names (chr string)
concseq<-function(fastafiles = fastafiles, outputname= "Concatenated_Seq.fas"){
  mylist<-list()
  
  for(i in 1:length(fastafiles)){
    #read in fastafile
    fastaFile<-Biostrings::readDNAStringSet(fastafiles[i])
    
    #save as df
    fastaFile.df<-as.data.frame(fastaFile)
    
    #add in column Sample ID
    fastaFile.df$SampleID<-rownames(fastaFile.df)
    
    mylist[[i]]<-fastaFile.df
    
  }
  
  joindf<-plyr::join_all(mylist, by = "SampleID")  ##works on lists
  
  # Remove SampleID column
  rownames(joindf)<-joindf$SampleID
  drops<-"SampleID"
  joindf<-joindf[ , !(names(joindf) %in% drops)]
  
  
  #unite final columns
  df_final <-tidyr::unite(joindf,col = "ConcSeq", sep = "", remove = TRUE)  #unites into single column ("NewSeq")
  outfilen<-file.path(getwd(), file.name=paste(outputname))
  seqRFLP::dataframe2fas(df_final, file=outfilen)
}


### Col deletion with frameshift ammendment to gene coords 
genome.coldel.only<-function(fasta, threshold = 1, csv = "TestData_coords.csv", prefix = ""){
  ## csv should have columns for Start and Stop (coordinates), GeneID (Gene duplications need to be uniquely labelled), GeneType ("PCG" or "NPGC" only), input files (csv and fasta = file path), threshold = proportion of samples allowed to not have data (e.g. 3/30 samples missing is maximum allowed = 0.1 = 90% cut-off) 
  
  ### Data
  coords<-read.csv(csv, header = T)
  coords[is.na(coords$GeneID),"GeneID"]<-LETTERS[seq(1,nrow(coords[is.na(coords$GeneID),]))] ##if there are any non-coding regions that are not named - label alphabetically
  
  
  x<-ape::read.dna(fasta, format="fasta")
  
  
  ### Identify which columns have N's in any of the samples (threshold = 0)
  if (!inherits(x, "DNAbin")) x <- ape::as.DNAbin(x)    #object-oriented method based on class of argument
  if (!is.matrix(x)) stop("DNA sequences not in a matrix")   #make sure DNAbin is in matrix format
  foo <- function(x) sum(x == 240)   #same as x[x==240], i.e. number of N's in alignment
  missing <- apply(x, 2, foo)   #apply the function foo to all the columns (2) in alignment (x); g <- vector of number of N's in each column (length = number of columns)
  i <- as.vector(which(missing / nrow(x) > threshold)) # i is vector of column numbers
  
  
  ### Identify whether it is a PCG or NPCG
  df2rm<-data.frame(i = as.numeric(i), region = NA)
  
  for(a in 1:nrow(df2rm)){
    df2rm$Start[a]<-coords[which(spatstat.utils::inside.range(df2rm$i[a], coords[c("Start","Stop")])), "Start"] 
  }
  
  for(b in 1:nrow(df2rm)){
    df2rm$Stop[b]<-coords[which(spatstat.utils::inside.range(df2rm$i[b], coords[c("Start","Stop")])), "Stop"] 
  }
  
  for(c in 1:nrow(df2rm)){
    df2rm$region[c]<-coords[which(spatstat.utils::inside.range(df2rm$i[c], coords[c("Start","Stop")])), "GeneType"] 
  }
  
  
  ### For PCG's remove in triplicate
  if(nrow(df2rm[df2rm$region == "PCG", ]) > 0){
    
    df2rmPCG <- df2rm[df2rm$region == "PCG", ]
    df2rmNPCG <- df2rm[df2rm$region == "NPCG", ]
    
    for(d in 1:nrow(df2rmPCG)){
      df2rmPCG$index[d]<-match(df2rmPCG$i[d], seq(df2rmPCG$Start[d],df2rmPCG$Stop[d]))
    }
    
    removed <- as.numeric(df2rmPCG$i) ##only PCG values
    indexed<-as.numeric(df2rmPCG$index)
    
    for (k in 1:length(removed)) {
      if (indexed[k] %% 3 == 0) {
        funthird <- function(x)
          c(x - 1, x - 2)
        removed <- append(removed, funthird(removed[k]))
        
      } else if ((indexed[k] + 1) %% 3 == 0) {
        funsecond <-  function(x)
          c(x - 1, x + 1)
        removed <- append(removed, funsecond(removed[k]))
        
      } else {
        funfirst <- function(x)
          c(x + 1, x + 2)
        removed <- append(removed, funfirst(removed[k]))
      }
    }
    
    removed<-sort(removed)
    removed<-unique(removed)
    removed<-na.omit(removed)
    df2<-data.frame(i = removed, region = "PCG")
    
    
    # Add to rest of columns to be deleted
    df2rm<-rbind(df2rmNPCG[1:2], df2)
  } 
  
  df2rm<-unique(df2rm)#make sure there's no duplicates
  df2rm<-df2rm[order(df2rm$i),]
  
  rownames(df2rm)<-NULL
  
  ### ID Genes
  df2rm$GeneID<-NA
  
  for(e in 1:nrow(df2rm)){
    df2rm$GeneID[e]<-coords[which(spatstat.utils::inside.range(df2rm$i[e], coords[c("Start","Stop")])), "GeneID"] 
  }
  
  
  ### Counts
  countdf<-plyr::count(df2rm, vars="GeneID") #output = dt with label + freq columns
  coords<-merge(coords, countdf, by="GeneID", all.x = T, sort=F)
  coords[is.na(coords$freq), "freq"]<- 0
  coords<-coords[order(coords$Start),]
  
  ### Ammend coordinates
  coords$Length<- coords$Stop - coords$Start +1
  coords$cumsum<-cumsum(coords$freq)
  coords$lenf<- coords$Length - coords$freq 
  
  coords$Start2<-coords$Start
  coords$Stop2<-coords$Stop
  
  for(f in 1:nrow(coords)){
    coords$Stop2[f]<-coords$Stop[f]-coords$cumsum[f]
  }
  
  for(g in 1:(nrow(coords)-1)){
    coords$Start2[g+1]<-coords$Start[g+1] - coords$cumsum[g] 
  }
  
  
  ### double check results:  
  #coords$len2<- coords$Stop2 - coords$Start2 + 1
  #coords$len2 == coords$lenf
  #which(coords$len2 != coords$lenf)
  
  rows_to_remove <- coords$lenf > 0  ###length of these should be added to cumsum
  coords<-coords[rows_to_remove,]   
  ###
  
  ### Remove columns
  i<-df2rm$i
  if (length(i)) x <- x[, -i] 
  
  cutoff<-(1-threshold)*100
  coords$Stop2[nrow(coords)]<-ncol(x)
  
  ape::write.dna(x, file=paste(prefix, paste(cutoff, "cutoff.fas", sep=""), sep="_"), format="fasta",  colw=1000000)
  write.csv(coords, paste(prefix, "coords", paste(cutoff, "cutoff.csv", sep=""), sep="_"), row.names = F)
  write.csv(df2rm, paste(prefix, "removed", paste(cutoff, "cutoff.csv", sep=""), sep="_"), row.names = F)
  
  coords
  
}
