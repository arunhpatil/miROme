# Stepwise processing tutorial to generate microRNAOme data

Here we describe the protocols implemented to gather cellular information and process a large scale expression analysis of cellular miRNA in humans. 

#### General Protocol to Achieve 2,077 samples (2,406 runs) of human cellular microRNA

**Step 1**.  Query at NCBI Sequence Read Archive (https://www.ncbi.nlm.nih.gov/sra)<br>
```
((miRNA[All Fields] OR microRNA[All Fields] OR (small[All Fields] AND RNA[All Fields]) AND ("Homo sapiens"[Organism] OR ("Homo sapiens"[Organism] OR human[All Fields]))) AND "Homo sapiens"[Organism] AND (cluster_public[prop] AND "library layout single"[Properties] AND 1900[MDAT] : 2900[MDAT] NOT "strategy epigenomic"[Filter] NOT "strategy genome"[Filter] NOT "strategy exome"[Filter] AND "filetype fastq"[Properties]))
``` 
Then download metadata through Run Selector as explained [here](https://github.com/NCBI-Hackathons/ncbi-cloud-tutorials/blob/master/SRA%20tutorials/tutorial_SRA_run_selector.md)
<br>
***Step 2***.
## Citation
A curated human cellular microRNAome based on 196 primary cell types. GigaScience 2022

