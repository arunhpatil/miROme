#!/bin/bash
#SBATCH --partition=preempt --time=2-00:00:00 --output=h_rundata2.log -c 40 --mem=120G
module load sratoolkit/2.9.2

/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041393
gzip DRR041393.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041530
gzip DRR041530.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041578
gzip DRR041578.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041592
gzip DRR041592.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041620
gzip DRR041620.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041356
gzip DRR041356.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041386
gzip DRR041386.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041430
gzip DRR041430.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041399
gzip DRR041399.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041408
gzip DRR041408.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041422
gzip DRR041422.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041436
gzip DRR041436.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041443
gzip DRR041443.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041450
gzip DRR041450.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp DRR041470
gzip DRR041470.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR3996364
gzip SRR3996364.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR3996365
gzip SRR3996365.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR3996366
gzip SRR3996366.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR3996367
gzip SRR3996367.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR3996368
gzip SRR3996368.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR3996369
gzip SRR3996369.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR5121485
gzip SRR5121485.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR5121486
gzip SRR5121486.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR5121487
gzip SRR5121487.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR1012333
gzip SRR1012333.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR1012334
gzip SRR1012334.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR7585373
gzip SRR7585373.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR6433711
gzip SRR6433711.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR6433712
gzip SRR6433712.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR6433713
gzip SRR6433713.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR6433714
gzip SRR6433714.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR6433715
gzip SRR6433715.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027315
gzip SRR4027315.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027317
gzip SRR4027317.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027318
gzip SRR4027318.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027319
gzip SRR4027319.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027326
gzip SRR4027326.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027328
gzip SRR4027328.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027329
gzip SRR4027329.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027372
gzip SRR4027372.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027373
gzip SRR4027373.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027392
gzip SRR4027392.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027395
gzip SRR4027395.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027396
gzip SRR4027396.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027397
gzip SRR4027397.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027398
gzip SRR4027398.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027422
gzip SRR4027422.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027423
gzip SRR4027423.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027424
gzip SRR4027424.fastq
/gpfs/fs1/sfw2/sratoolkit/2.9.2/bin/fasterq-dump -e 40 -t temp SRR4027425
gzip SRR4027425.fastq
