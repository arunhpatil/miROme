# General Protocol to Achieve 2,077 samples (2,406 runs) of human cellular microRNA

Here we describe the protocols implemented to gather cellular information and process a large scale expression analysis of cellular miRNA in humans. 

### Obtaining metadata
**Step 1**.  Query at NCBI Sequence Read Archive (https://www.ncbi.nlm.nih.gov/sra)<br>
```
((miRNA[All Fields] OR microRNA[All Fields] OR (small[All Fields] AND RNA[All Fields]) AND ("Homo sapiens"[Organism] OR ("Homo sapiens"[Organism] OR human[All Fields]))) AND "Homo sapiens"[Organism] AND (cluster_public[prop] AND "library layout single"[Properties] AND 1900[MDAT] : 2900[MDAT] NOT "strategy epigenomic"[Filter] NOT "strategy genome"[Filter] NOT "strategy exome"[Filter] AND "filetype fastq"[Properties]))
``` 
Then download metadata through Run Selector as explained [here](https://github.com/NCBI-Hackathons/ncbi-cloud-tutorials/blob/master/SRA%20tutorials/tutorial_SRA_run_selector.md).
<br>
<br>

### Clean the metadata and selecting SRA runs
**Step 2**. Manually curate the run list (58,117) to positively select samples that appeared to be from human primary cells.
<br>
<br>

### Downloading SRA runs
**Step 3**. Use fasterq-dump from the [NCBI SRA-toolkit](https://hpc.nih.gov/apps/sratoolkit.html) on each run.  This was performed by using a Python script (`createDownload.py`) to create a shell script to sequentially download each to a computer cluster.  The fastq files were converted to fastq.gz using gzip command (as part of the shell script).<br><br>
Example to download a single file is shown below:<br> 
> `fasterq-dump -e 40 -t temp DRR041393` <br>
> where `-e` specifies number of parallel executions and `-t` specifies the temporary folder. <br>
<br>

To download bulk files, you need to run the python script, `createDownload.py` as shown below: <br>
The script, test input and output shell script can be found [here](https://github.com/mhalushka/miROme/tree/main/Processing_Tutorial/).

> `python createDownload.py SRR_Accessions.txt`
> 
The SRR_Accessions.txt is a text input file which contains SRR accessions that need to be downloaded and an output file `download_runs.sh` is created. This download_runs.sh is a shell script and needs to be executed as shown below:
> `bash download_runs.sh` <br>
> The output of the command should look as shown below:
```
spots read      : 2,301,241
reads read      : 2,301,241
reads written   : 2,301,241
spots read      : 3,782,536
reads read      : 3,782,536
reads written   : 3,782,536
spots read      : 2,906,968
reads read      : 2,906,968
reads written   : 2,906,968
(These messages will continue until the end of the file, only subset is shown above)
```
> Note: The header lines are specific to the server used during this project, if this does not apply to your work station, please remove them before execution;

<br>

### Fetching adapter sequences 
**Step 4**. A python script (`adapter_detect.py`) was created to identify the first 5 rows of let-7a with adaptor sequence across all downloaded files. These were manually identified for which adaptor sequence type was used.  All samples with 4N, 5’, or UMI-based adaptors were excluded as they would not work with the optional miREC step of miRge3.0. All other adaptor types were identified and adaptor sequence was supplemented into the miRge3.0 parameters for accurate processing of the fastq.gz files. 

Move all the downloaded folders in a new directory `SRR_folder`. <br>
```
mkdir SRR_folder
mv *.fastq.gz ./SRR_folder
```

Execute the python script `adapter_detect.py` script, as shown below: <br>
`python adapter_detect.py SRR_folder > new_Adapters_nextbatch.txt` <br>
The python script takes in a folder name as input and iterates through all of the fastq.gz files and fetches let-7a sequences. The output is redirected to a file called new_Adapters_nextbatch.txt. The output is a tab delimited file with two columns. The first column contains the let-7a + adapter sequence and the second column is SRR runs.
```
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036697.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036697.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036697.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036697.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036697.fastq.gz
CTGGTGAGGTAGTAGGTTGTATAGTTCTGTAGGCAC    SRR2038610.fastq.gz
CTGGTGAGGTAGTAGGTTGTATAGTTCTGTAGGCAC    SRR2038610.fastq.gz
CTGGTGAGGTAGTAGGTTGTATAGTTCTGTAGGCAC    SRR2038610.fastq.gz
CTGGTGAGGTAGTAGGTTGTATAGTTCTGTAGGCAC    SRR2038610.fastq.gz
CTGGTGAGGTAGTAGGTTGTATAGTTCTGTAGGCAC    SRR2038610.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036709.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036709.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036709.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCCGGGTG    DRR036709.fastq.gz
TGAGGTAGTAGGTTGTATAGTTTGGAATTCTCGGGT    DRR036709.fastq.gz
```
From the example above, DRR036697 and DRR036709 represent runs with Illumina adapters (`TGGAATTCTCGGGT`). However, SRR2038610 has a 4N nucleotide adapter on either end of let-7a and is therefore excluded in this analysis. 

### Executing miRge3.0 on the SRA runs
**Step 5**. miRge3.0 was performed on 4,184 runs using this general command: <br>
`miRge3.0 -s SRAS-file.fastq.gz -a <adapter_sequence> -gff -bam -trf -lib miRge3_Lib -on human -db miRBase -o OutputDir -mEC -ks 20 -ke 20`. <br>
Generally, 11 runs were within a single miRge3.0 run. An example run is shown below:<br>

```
miRge3.0 -s SRR2061941.fastq.gz,SRR2061942.fastq.gz,SRR2061943.fastq.gz,SRR2061944.fastq.gz,SRR2061945.fastq.gz,SRR2061946.fastq.gz,SRR2061947.fastq.gz,SRR2061948.fastq.gz,SRR2061949.fastq.gz -a TGGAATTCTCGGGTGCCAAGGAACTCCAG  -gff -bam -trf -lib miRge3_Lib -on human -db mirbase -o epiC_out -mEC -ks 20 -ke 20
```

## Citation
A curated human cellular microRNAome based on 196 primary cell types. GigaScience 2022

## Resources
1. Arun H Patil, Marc K Halushka. **miRge3.0: a comprehensive microRNA and tRF sequencing analysis pipeline**. [NAR Genomics and Bioinformatics]( <https://academic.oup.com/nargab/article/3/3/lqab068/6325159>). 2021.
2. NCBI Sequence Read Archive (https://www.ncbi.nlm.nih.gov/sra)
3. [NCBI SRA-toolkit](https://hpc.nih.gov/apps/sratoolkit.html) 

