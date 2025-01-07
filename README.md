# operOn.and.On
Leanie Williams

2025 01-07

## Introduction

- Can we determine bacterial operon structure using the GFF and genome information?
- Can operons be used as the features for transcriptomic experiments? 

## Methods

*Work in progress*
rules:
1. same strand
2. serially located
	- overlap of less than 500
3. no fight club


1. features in an operon must be co-oriented. 
	- co oriented would mean that the features are on the same strand (+ or -)
2. features must be found serially.
	- serially means that the number of bases between features should not exceed 500 bases
(*note* we will examine the threshold of bases and how it impacts the output)

** input **: gff
** output ** : csv
---

Testing:
**input**: RdKW20.annot.gff
- sequence-region 1 1830138
- species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=727

## Results

- argument parser
- define feature name, BED attributes (chromosome, start, stop, name, score, strand, ID, feature type, description etc etc...) 
