# operOn.and.On
Leanie Williams

2025 01-07

## Introduction

- Can we determine bacterial operon structure using the GFF and genome information?
- Can operons be used as the features for transcriptomic experiments?
- How do we find promoters?
- What are the regulators?

## Methods

*Work in progress*

1. Features in an operon must be co-oriented. 
	- co oriented would mean that the features are on the same strand (+ or -)
2. Features must be found serially.
	- serially means that the number of bases between features should not exceed 500 bases
(*note* we will examine the threshold of bases and how it impacts the output)


**input**: gff
**output** : csv



---

Testing:

**input**: RdKW20.annot.gff
- sequence-region 1 1830138
- species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=727

## Results

### GFF to Operons
We define functions for parsing the GFF, finding operons based on our rules (see Methods), and returning the results as a csv. The csv contains a header  `Feature, Locus_tag, Gene_type, Operon_ID, Operon_Start, Operon_End, Strand` and each feature (defined as genes or pseudogenes) is specified within an operon. 

The user can include a separation threshold (`--sep_thresh`) which will specify how many bases between features. This parameter has a range between 0 and 500 bases, and a default of 500. This parameter informs the results in the columns `oriented` and `oriented_nearby`. If a feature is located on the '+' strand (sense), they receive a 1 and the '-' strand (antisense) is 0. The results of `oriented_nearby` are `TRUE` if the features are both on the same strand and located within the defined separation threshold, defining the operon. 

### BED formatted operons
The user may return an optional [[BED]](http://useast.ensembl.org/info/website/upload/bed.html) formatted file (`--bed`).
The user may return a bed file for intergenic regions using `--inter`. 

**To do**: 

- argument parser
- define feature name, BED attributes (chromosome, start, stop, name, score, strand, ID, feature type, description etc etc...)


