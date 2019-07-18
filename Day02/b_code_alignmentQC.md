## General QC of alignment

We will make one directory per sample and collect the different QC outputs in the respective folder. 
`MultiQC` will eventually combine all the different QC stats into one interactive document that will use the naming scheme of the sub-directories.

### Prepare for alignment QC

```
mkdir alignment_qc
cd alignment_qc
mkdir SNF2_1 WT_1
```

### STAR log files

Just need to be copied from the alignment folder:

```
ln -s ~/mat/precomputed/results_alignment/ bam_files
for SAMPLE in SNF2_1 WT_1
do
	for i in bam_files/${SAMPLE}/*Log*final.out
	do
		cp $i ${SAMPLE}/STAR_log.final.out
	done
done
```

### samtools flagstat

For each read, this tool will look at the FLAG (field 4) and count the types
of settings.

```
export PATH=~/mat/software/samtools-1.7/bin/:$PATH

for SAMPLE in SNF2_1 WT_1
do
	for i in bam_files/${SAMPLE}/*bam
	do
	# we need to index the BAM files first, this will create a .bai file
		samtools index $i
	# afterwards, we can run samtools flagstat
		samtools flagstat $i >  ${SAMPLE}/flagstat.out
	done
done
```

Example output via `head WT_1/flagstat.out`:

```
1182835 + 0 in total (QC-passed reads + QC-failed reads)
133369 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
1182835 + 0 mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped

```

### RSeQC

#### bam stat

General bam file statistics (not RNA-seq specific) similar to `samtools flagstat`

```
PATH=~/mat/software/anaconda2/:$PATH

for SAMPLE in WT_1 SNF2_1
do
	bam_stat.py -i bam_files/${SAMPLE}/*bam > ${SAMPLE}/rseqc_bam_stat.out
done
```

Example output via `head -n10 WT_1/rseqc_bam_stat.out`:

```
#==================================================
#All numbers are READ count
#==================================================

Total records:                          1182835

QC failed:                              0
Optical/PCR duplicate:                  0
Non primary hits                        133369
```


#### Read distribution

How many reads fall into exons, introns, etc. 
This will be based on the annotation file that you provide!

```
REF_DIR=~/mat/referenceGenomes/S_cerevisiae/

for SAMPLE in WT_1 SNF2_1
do
read_distribution.py \
    -r ${REF_DIR}/sacCer3.bed \
    -i bam_files/${SAMPLE}/*.bam > ${SAMPLE}/rseqc_read_distribution.out
done
```
`head -n10 WT_1/rseqc_read_distribution.out `:

```
Total Reads                   1049466
Total Tags                    1059871
Total Assigned Tags           992608
=====================================================================
Group               Total_bases         Tag_count           Tags/Kb             
CDS_Exons           8832031             990363              112.13            
5'UTR_Exons         0                   0                   0.00              
3'UTR_Exons         0                   0                   0.00              
Introns             69259               630                 9.10              
TSS_up_1kb          2421198             1260                0.52       
```

### Gene body coverage

RNA-seq-specific QC, assesses whether there's a skew towards preferential coverage of either 5' or 3' end.

Takes a relatively long time to compute.

```
for SAMPLE in WT_1 SNF2_1
do
geneBody_coverage.py \
	-i bam_files/${SAMPLE}/*.bam \
	-r ${REF_DIR}/sacCer3.bed \
	-o ${SAMPLE}/rseqc_geneBody_coverage.out &
done
```


### QoRTs

As an alternative to `MultiQC`, `QoRTs` can be used to summarize multiple QC checks in one file.

```
# to call up the help manual for the QC module of QoRTs:
java -Xmx4g -jar ~/mat/software/qorts.jar QC man
```

`QoRTs` produces a lot of individual files if run in `QC` mode.

```
for SAMPLE in WT_1 SNF2_1
do
java -Xmx4g -jar ~/mat/software/qorts.jar QC \
     --singleEnded  \
	 --generatePdfReport \
     bam_files/${SAMPLE}/*.bam ${REF_DIR}/sacCer3.gtf $SAMPLE

done
```

The inclusion of QoRTs doesn't fully work yet.
Look at example QoRTs output at chagall.med.cornell.edu


### MultiQC

To summarize all the QC results into one document, we will use `MultiQC`, which will generate an interactive html report based on all the QC stats we've accumulated so far.

For simplicity, it's easiest to gather all QC stats per sample.

```
mkdir multiQC
cd multiQC
```

Let's link the individual stats.

```
ln -s ../alignment_results/WT_1
ln -s ../alignment_results/SNF2_1
```

Run MultiQC:

```
~/mat/software/anaconda2/bin/multiqc . --dirs  
```
```
[INFO   ]           rseqc : Found 2 read_distribution reports
[INFO   ]           rseqc : Found 2 gene_body_coverage reports
[INFO   ]           rseqc : Found 2 bam_stat reports
[INFO   ]        samtools : Found 2 flagstat reports
[INFO   ]            star : Found 2 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```

To look at the html file, you'll have to download it to your computer.
In a separate terminal, navigate to a directory of your liking, e.g. `~/Downloads`.
Then copy the html file from the server:

```
scp <user name>@<server IP>:<path to>/multiQC/*html .
```
