# rsomics-bed-fisher

Fisher's exact test for interval overlap significance between two BED files.

## Usage

```
rsomics-bed-fisher -a peaks.bed -b promoters.bed -g hg38.genome
```

## Origin

Equivalent to `bedtools fisher`. Implementation follows the bedtools algorithm:
builds a 2×2 contingency table from interval overlap counts and genome size,
then applies Fisher's exact test (hypergeometric distribution, three alternatives).

Upstream: [bedtools](https://github.com/arq5x/bedtools2) (MIT).  
Reference: Quinlan & Hall (2010). BEDTools. Bioinformatics 26(6): 841–842.
DOI: [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)

License: MIT OR Apache-2.0.
