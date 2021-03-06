## Code used in "Chapter 4: Patterns of Didemnum vexillum haplotype diversity at the mitochondrial NAD1 marker across four sites on the west coast of Scotland."

`hap.custom.R` : a custom function that uses both package `haplotypes` and `pegas` (<v1.0-1) to output a 'haplotype,DNAbin' object that can be used to make haplotype networks using `pegas` and generate estimates of haplotype diversity. 
This function has been validated for further use, and example usage and full function details are presented in: `hap.custom_doc.md`


Note: to recreate the haplotype maps in the full analysis code you will need an API Google key. 

Package minimum requirements: 

- dplyr *v1.0.2*
- ggmap *v3.0.0*
- ggplot2 *v3.3.3*
- ggrepel *v0.8.2*
- ggthemes *v4.2.0*
- haplotypes *v1.1.2*
- measurements *v1.4.0*
- pegas *v0.11*   **2021 update (v1.0-1) does not function the same - this only impacts plotting the haplotype network in this case**
- phylotools *v0.2.2*
- purrr *v0.3.4*
- rgdal *v1.5-18*
- tidyr *v1.1.2*
