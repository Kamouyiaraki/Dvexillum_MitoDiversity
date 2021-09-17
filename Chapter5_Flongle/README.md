## Code used in "Chapter 5:  Barcoding and metabarcoding of the non-native invasive tunicate, Didemnum vexillum, using the Flongle by Oxford Nanopore Technology."

For the 6 experiments in this chapter, the code needed can be adapted from one of the following (all available in `Ch5_full_code.yml`):

1. **Single sample haplotyping analysis using a single Sanger reference sequence (Experiment 1)**
2. ***In silico* pseudorandom sampling of starting HAC reads (Experiment 1)**
    - includes **Example test**: mapping a pseudorandom subset of reads to a single incorrect reference to test consensus sequence accuracy
4. **Multiplexed sample haplotyping (Experiments 2-6)**
5. **Analysis of mapping results from a multiplexed run (haplotype classification) in R (Experiments 2-4)**
6. **Analysis of mapping results from a metabarcoding multiplexed run (Experiments 5 & 6)**


Everything in this chapter was done on a unix machine (Bio-Linux), except for guppy re-basecalling was done on a Windows machine. 

### Version details: 

- Bio-Linux 8.0.7
  - Linux 3.13.0-170-generic
- R version 4.0.2 (2020-06-22) -- "Taking Off Again"
- NanoPlot 1.35.1
- NanoStat 1.1.2
- NanoFilt 2.8.0
- MinIONQC version 1.4.1
- conda 4.8.3
- minimap2 2.17-r941
- bcftools 1.9
- samtools 1.7
- Python 3.6.13
- Guppy - Version 4.4.2+9623c1626 
