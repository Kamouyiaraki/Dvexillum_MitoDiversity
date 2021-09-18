################# FUNCTION hap.custom()  v0.3
library(pegas)
library(haplotypes)

##haplotype formats from package haplotypes are not compatible with that from package pegas
##creates haplotype DNAbin compatible with pegas, but allows for indels to be dealt differently as with haplotypes::haplotype()
##Used to read in FASTA only at the moment
##Also allows for haplotype diversity estimates (hapdiv, based on pegas::hap.div() ) to be made on all samples in specified populations (pops) as well as +/- variance. 
##Populations are specified as a vector of length = number of sequences in your FASTA
##output is haplotype DNAbin & table of results with two rows added to the bottom of a table structured: cols = defined populations, rows = frequencies of each haplotype, last two rows = hap.div estimate and variance
##output of hap.custom() can be saved to be used for other pegas functions (e.g. hap)

###Note: haplotype diversity (pegas::hap.div()) calculated as: 
#h = n(1- sum (p^2))/(n-1)
#where n = sample size
#where p = the frequency (prop) of each haplotype
#Originally presented in Nei 1975, Nei 1978, and Nei and Tajima, 1980 


hap.custom<- function(fastafile, indel="5th", labels=NULL, hapdiv=FALSE, pops=NULL, variance=FALSE){
  #indels = argument in haplotypes::distance() on how to treat indels
  #labels = what to use to label each haplotype ID (pegas::haplotype), default = roman numerals
  #haplotype diversity stats
  #populations defined for each population
  #variance

  ### Read FASTA file and calculate distance matrix ("N")  
  fas<-haplotypes::read.fas(fastafile)
  
  d<-distance(fas,indels=indel)
  if(length(d)==0) stop("at least two DNA sequences are required") #total number of comparisons (i.e. do we have only one haplotype? If true then can't move forward)
  
  
  ### Run modified haplotype function based on combined haplotypes and pegas packages
  
  ## Part 1=Haplotypes::haplotypes source code:
  
  x <- as.matrix(d) #turn distance (class dist()) into matrix
  diag(x)<-0   #fill diagonal
  nseq<-nrow(x)  #total number of sequences
  whap<-x[1,]==0  #T/F of which in the top row have zero distance (same haplotype as the first sample = True, different haplotype = FALSE)
  haploind<-list(which(whap))   #column ID of which are the which are the same haplotype as the first sample (as a list) i.e. only 1 haplotype displayed
  haplovec<-which(whap)  ##column ID of which are the which are the same haplotype as the first sample (as a vector)
  
  x<-x[-haplovec,]
  nseq<-nrow(x)
  
  
  if(length(x)>1){   #total number of cells in matrix
    
    for(i in 1:nseq)  #repeat pairwise comparison for the rest of the rows (samples)
    { 
      whap<-x[i,]==0      #pairwise comparison (row of matrix)
      whap[haplovec]<-FALSE   #keep only comparisons that haven't been used so far
      haploind<-c(haploind,list(which(whap)))  #combines list of haploind, with ID's that are identical to the rest of the samples (list of all samples containing ID matches)
      haplovec<-unique(c(haplovec,which(whap)))  #c(haplovec, which(whap)) = combines vector for first row with all other sample comparisons (length= number of samples) and makes sure there are no duplicates
      haplovec<- unique(haplovec)
      }	
  }
  
 

  
   
  empthaplo<-sapply(haploind,length) #vector of frequencies (i.e. number of ID's in each item of the list)
  haploind<-unique(haploind[empthaplo>0])  #makes sure that you are only keeping unique comparisons (i.e. full list of all haplotypes)
  haplolistint<-lapply(haploind, as.integer)  #row numbers as integer
  haplolist<-lapply(haploind,names)   #list originally contains names + row number; this replaces that with just names
  
  uniqhapindex<-sapply(haploind,"[",1)  #unlists into the first value of every haplotype (Na's for remaining samples)
  

  
  freq<-sapply(haplolist,length) #vector of frequencies with number of ID's in each item of the list > 0
  hapnum<-length(freq)
  
  ##This part is removed from haplotype::haplotypes() 
  #hapdistmat<-as.matrix(x[uniqhapindex,uniqhapindex])  #turns into matrix of pairwise comparisons 
  #names(haploind)<-paste("haplotype",1:hapnum,sep="")
  #names(haplolist)<-names(haploind)
  #hapobj<-new("Haplotype",haplist=haplolist,hapind=haploind,uniquehapind=uniqhapindex,d=hapdistmat,freq=freq,nhap=hapnum)
  #hapobj
  
  
  #############################################  haplotype diversity: variable pops, variable variance (T/F), variable hapdiv(T/F)
  if(hapdiv){
    if(length(pops)>0){
      
      
      factors<-pops
      
      flevels<-levels(factor(factors))    #population levels as factor (specify in function)
      hapmat<-matrix(0, hapnum, length(flevels))   #create a matrix (columns = populations, rows = each haplotype)
      hapvec<-vector("numeric",length(factors)) #an empty vector for the number of sequences
      rownames(hapmat)<-1:hapnum    #name rownames
      colnames(hapmat)<-flevels  #name column names
      
      for(m in 1:hapnum)   ### for each row (i.e. haplotype), add in number of samples in each column (i.e. population)
      {
        hind<- haploind[[m]]      #all individuals of each type (haploind)
        fstab<-table(factor(factors[hind]))    #subsets factors (i.e. population IDs specified) for the numbers that match rownumbers in haploind and counts
        hapmat[m,names(fstab)]<-fstab     #puts frequencies into matrix
      }
      
      
      hapdf<-data.frame(Pop = colnames(hapmat), Hap.Div = rep(NA, ncol(hapmat)), Variance = rep(NA, ncol(hapmat)))
      
      for(b in 1:ncol(hapmat)){  
        f<-hapmat[,b]
        n <- sum(f)  #total number of samples
        p <- f/n   #proportion of samples in each haplotype (frequencies)
        sump2 <- sum(p^2)  #sum of squared frequencies
        n1 <- n - 1L   #N-1
        hapdf[b,"Hap.Div"]  <- (1 - sump2) * n / n1
        if (variance) {
          tmp <- sump2^2
          sump3 <- sum(p^3)
          hapdf[b,"Variance"] <- (sump2 - tmp + 4 * n1 * (sump3 - tmp))/(n * n1)
        }
      }
      res<-hapdf
    }else{
      n <- sum(freq)
      p <- freq/n
      sump2 <- sum(p^2)
      n1 <- n - 1L
      res <- (1 - sump2) * n / n1
      if (variance) {
        tmp <- sump2^2
        sump3 <- sum(p^3)
        var <- (sump2 - tmp + 4 * n1 * (sump3 - tmp))/(n * n1)
        res <- c(res, var)
      }        
      namevec<- c("Hap.Div", "Variance")
      names(res)<- (namevec[1:length(res)])
      
    }
    print(res)
  }  
  #############################################  
  
  
  ## Part 2= change for pegas format:  class(obj) = "haplotype" "DNAbin"  ##turns alignment into alignment of only the two haplotypes + relative frequencies
  
  dnabinFAS<-read.dna(fastafile, format="fasta")
  
  obj <- dnabinFAS[uniqhapindex, ] #extract only unique sequences
  
  if (is.null(labels)) labels <- as.character(as.roman(seq_along(haploind)))  #creates labels for each haplotype identified
  rownames(obj) <- labels  #label sequences by haplotype ID
  class(obj) <- c("haplotype", "DNAbin")
  #attr(obj, "index") <- lapply(i, function(x) which(h == x))
  attr(obj, "index") <- haplolistint
  
  nms.x <- deparse(fastafile)   #name of fastafile used 
  attr(obj, "from") <- nms.x    #add the name of the fastafile as an attribute to object

  obj
}

