---
title: "hap.custom() Example"
author: "Maria Kamouyiaros"
date: "19/07/2021"
---


## A short intro

The package `pegas` can be used to create haplotype networks in R. However, this package does not deal well with any indels or missing data. The package `haplotypes` has the functionality to deal with that issue.

Object format "haplotype" that are generated using the package `haplotypes` is not compatible with the package `pegas`. 
`hap.custom()` creates a "haplotype,DNAbin" compatible with `pegas`, but allows for indels to be dealt differently as with `haplotypes::haplotype()` (for example specifying deletions as a 5th character state using `indel = "5th" `)


From the `haplotypes` R documentation: 

>"the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.
>
>sic: Treating gaps as a missing character and coding them separately following the simple indel coding method.
>
>5th: Treating gaps as a fifth state character.
>
>missing: Treating gaps as a missing character."


Additionally the `hap.custom()` function has the ability to estimate haplotype diversity (+/- variance) of all samples or haplotype diversity (+/- variance) of specified populations at the same time. 


## `hap.custom()`

The syntax is as follows: 

`hap.custom(fastafile, indel="5th", labels=NULL, hapdiv=FALSE, pops=NULL, variance=FALSE)`


- fastafile = file path of your FASTA file (at the moment this reads in FASTA files only)
- indel = indel coding method used by `haplotypes` package ("sic", "missing" or "5th")
- labels = haplotype ID label coding method used by `pegas::haplotype` (default = roman numerals) 
- hapdiv = T/F whether to include haplotype diversity estimates
- pops = population IDs for haplotype diversity estimates; a vector of length = number of sequences in your FASTA (only works if hapdiv = TRUE)
- variance = T/F whether to include +/- variance (only works if hapdiv = TRUE)



**Additional information:**

Haplotype diversity (`pegas::hap.div()`) is calculated as: `h = n(1- sum (p^2))/(n-1)`

- where n = sample size
- where p = the frequency (prop) of each haplotype

Originally presented in Nei 1975, Nei 1978, and Nei and Tajima, 1980. 


## Output 

Output is a **haplotype,DNAbin** object which can be used in other `pegas` functions, e.g. plotting haplotype networks. 

<br>
Additional results are printed out in the console: 

- Haplotype ID labels and their frequencies as a printed table (columns = defined populations, rownames = haplotype ID, rows = frequencies of each haplotype)
- Haplotype Diversity estimates (if `hapdiv = TRUE`)

<br>

## Warning messages 

There is only one potential warning message: 

>"You haven't defined your populations!"

if `pops=NULL`. 


```{r source, include=FALSE}
setwd("../Desktop")
source("hap_cust.R")
```

## Example use:

In these examples I'm using a fake dataset: 

- 9 samples 
- from 3 populations
- a total of 4 haplotypes



```{r exampleFAS, echo=FALSE}
exampleFAS <- "hap4_pop3.fas"
dnabinFAS<-read.dna(exampleFAS, format="fasta")
dnabinFAS
```


- 2 of the samples (pop2) in this set have a codon deletion


### Example 1: I just want the DNAbin format and I don't want to treat dels as a 5th base

```{r example1, echo=T}
hap.custom(fastafile = exampleFAS, indel="missing")
```

### Example 2: I just want the DNAbin format and I want to treat dels as a 5th base

```{r example2, echo=T}
hap.custom(fastafile = exampleFAS, indel="5th")

```


### Example 3: I want to treat dels as a 5th base and I want haplotype diversity 

```{r example3, echo=T}
hap.custom(fastafile = exampleFAS, indel="5th", hapdiv = T)

```


### Example 4: I want to treat dels as a 5th base and I want haplotype diversity + variance

```{r example4, echo=T}
hap.custom(fastafile = exampleFAS, indel="5th", hapdiv = TRUE, variance = T)

```




### Example 5: I want to treat dels as a 5th base and I want haplotype diversity with specified populations

```{r example5, echo=T}
seqnames<-c("pop1", "pop1", "pop1", "pop2", "pop2", "pop3", "pop3", "pop3", "pop3")
hap.custom(fastafile = exampleFAS, indel="5th", hapdiv = T, variance=T, pops=seqnames)

```


### Example 6: I want to use the DNAbin output to plot a haplotype network using pegas (including v1.0-1) 

(I still want to treat dels as a 5th base and I want haplotype diversity with specified populations)


```{r example6, echo=TRUE, warning=FALSE}

seqnames<-c("pop1", "pop1", "pop1", "pop2", "pop2", "pop3", "pop3", "pop3", "pop3")
ExHaps <- hap.custom(fastafile = "hap4_pop3.fas", indel = "5th", hapdiv = T, pops = seqnames, variance=T )

ExNet <- pegas::rmst(ExHaps)
plot(ExNet)

```

