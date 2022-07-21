# Stepwise processing tutorial to generate microRNAOme data

Here we describe the protocols implemented to gather cellular information and process a large scale expression analysis of cellular miRNA in humans. 

#### General Protocol to Achieve 2,077 samples (2,406 runs) of human cellular microRNA

**Step 1**.  Query at NCBI Sequence Read Archive (https://www.ncbi.nlm.nih.gov/sra)<br>
```
((miRNA[All Fields] OR microRNA[All Fields] OR (small[All Fields] AND RNA[All Fields]) AND ("Homo sapiens"[Organism] OR ("Homo sapiens"[Organism] OR human[All Fields]))) AND "Homo sapiens"[Organism] AND (cluster_public[prop] AND "library layout single"[Properties] AND 1900[MDAT] : 2900[MDAT] NOT "strategy epigenomic"[Filter] NOT "strategy genome"[Filter] NOT "strategy exome"[Filter] AND "filetype fastq"[Properties]))
``` 
Then download metadata through Run Selector as explained [here](https://github.com/NCBI-Hackathons/ncbi-cloud-tutorials/blob/master/SRA%20tutorials/tutorial_SRA_run_selector.md)
<br>
<br>
**Step 2**. Manually curate the run list (58,117) to positively select samples that appeared to be from human primary cells.
<br>
<br>
**Step 3**. Use fasterq-dump from the [NCBI SRA-toolkit](https://hpc.nih.gov/apps/sratoolkit.html) on each run.  This was performed by using a Python script to create a shell script to sequentially download each to a computer cluster.  The fastq files were converted to fastq.gz using gzip command (as part of the shell script).<br><br>
Example to download a single file is shown below:<br> 
> `fasterq-dump -e 40 -t temp DRR041393` <br>
> where `-e` specifies number of parallel executions and `-t` specifies temporary folder. <br>
<br>

To download bulk files, you need to run the python script, `createDownload.py` as shown below: 
The script, test input and output shell script can be found [here](https://github.com/mhalushka/miROme/tree/main/Processing_Tutorial/)

> `python createDownload.py SRR_Accessions.txt`
> 
The SRR_Accessions.txt is a text input file which contains SRR accessions that need to be downloaded and an output file `download_runs.sh` is created. This download_runs.sh is a shell script and need to be executed as shown below:
> `bash download_runs.sh` <br>
> Note: The header lines are specific to server used during this project, if this doesnot apply to your work station, please remove them before execution;

<br>

**Step 4**. A python script (`adapter_detect.py`) was created to identify the first 5 rows of let-7a with adaptor sequence across all downloaded files. These were manually identified for which adaptor sequence type was used.  All samples with 4N, 5’ or UMI-based adaptors were excluded as they would not work with the miREC step of miRge3.0. All other adaptor types were identified and adaptor sequence was supplemented into the miRge3.0 parameters for accurate processing of the fastq.gz files. 

Assuming all your downloaded runs  execute the adapter_detect.py script, please 
Move all the donwloaded folders in a new directory `SRR_folder` <br>
```
mkdir SRR_folder
mv *.fastq.gz ./SRR_folder
```


## Citation
A curated human cellular microRNAome based on 196 primary cell types. GigaScience 2022

