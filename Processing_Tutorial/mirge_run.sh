
# To execute miRge3.0, change the directory to SRR_folder (Alternatively, one could provide absolute paths)
miRge3.0 -s SRR3996364.fastq.gz,SRR3996365.fastq.gz,SRR3996366.fastq.gz,SRR3996367.fastq.gz,SRR3996368.fastq.gz,SRR3996369.fastq.gz -gff -bam -trf -lib /mnt/d/Halushka_lab/Arun/GTF_Repeats_miRge2to3/miRge3_Lib/revised_hsa -on human -db mirbase -o OutputDir -mEC -ks 20 -ke 20

miRge3.0 -s DRR041356.fastq.gz,DRR041386.fastq.gz,DRR041393.fastq.gz,DRR041399.fastq.gz,DRR041408.fastq.gz,DRR041422.fastq.gz,DRR041430.fastq.gz,DRR041436.fastq.gz,DRR041443.fastq.gz  -a illumina  -gff -bam -trf -lib miRge3_Lib -on human -db mirbase -o OutputDir -mEC -ks 20 -ke 20

miRge3.0 -s DRR041450.fastq.gz,DRR041470.fastq.gz,DRR041530.fastq.gz,DRR041578.fastq.gz,DRR041592.fastq.gz,DRR041620.fastq.gz,SRR5121485.fastq.gz,SRR5121486.fastq.gz,SRR5121487.fastq.gz -a illumina  -gff -bam -trf -lib miRge3_Lib -on human -db mirbase -o OutputDir -mEC -ks 20 -ke 20
