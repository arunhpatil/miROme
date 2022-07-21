#!/usr/bin/python
import sys
# This python script will create shell script to queue downloads sequentially 
# The commands to download the runs are listed in the following file, user may change this name according to the needs or use sys.argv arguments to call it from the command line.
out = open("download_runs.sh", "w+")
out.write("#!/bin/bash") # Indicating that we are calling shell script.
out.write("#SBATCH --partition=preempt --time=2-00:00:00 --output=h_rundata2.log -c 40 --mem=120G") # This is specific to BlueHive cluster at RU, this is optional for your needs. 
out.write("module load sratoolkit/2.9.2") # loading module sratoolkit to call fasterq-dump, this is specific to the BlueHive cluster, user should skip this line if this is not similar 
out.write("")
userinfile = sys.argv[1]
# read each accession and write the command to download these accessions from NCBI
with open(userinfile, "r") as inf: 
	infile = inf.readlines()
	for i in infile:
		x = i.strip()
		out.write("/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp " + str(x)) # This will specify the path to fasterq-dump, if you have path set globally, you could skip the absolute path and retain only fasterq-dump onwards 
		out.write("gzip "+str(x)+".fastq") # This will zip the files to save storage space, if you have enough room in the server/workstation, you may skip this as well by commenting it out. 
		#print("rm -rf /home/apatil8/ncbi/public/sra/*") # optional if the user wants to clear the temporary space manually then they can uncomment this line.
