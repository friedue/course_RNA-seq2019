# DAY 2

## Alignment with STAR

for alignment, two things are needed:

    1. reference genome
    2. annotation

The genome index was already generated! --> see course notes.

We're just going to do the alignment -- for all samples of `WT_1`


    runSTAR=~/mat/software/STAR-2.7.1a/bin/Linux_x86_64/STAR
    REF_DIR=~/mat/referenceGenomes/S_cerevisiae
    mkdir ~/class_Summer2019/alignment_STAR


* __list fastq.gz__ files separated with comma, but no whitespaces

		for i in ~/mat/precomputed/rawReads_yeast_Gierlinski/WT_1/*fastq.gz
		do
			FILES=`echo $i,$FILES`
		done

        $runSTAR --genomeDir ${REF_DIR}/STARindex/ \
                  --readFilesIn $FILES \
                  --readFilesCommand zcat \ # gzipped fastq files
                  --outFileNamePrefix alignment_STAR/WT_1_ \
                  --outFilterMultimapNmax 1 \ # only reads with one location will be returned
                  --outReadsUnmapped Fastx \ # extra output file for unaligned reads (Fasta)
                  --outSAMtype BAM SortedByCoordinate \
                  --runThreadN 1 \
                  --twopassMode Basic


Most important things to consider:

* handling of __multi-mapped__ reads (e.g., how the best alignment score is assigned and the number and order in which secondary alignments are reported);
* optimization for __very small genomes__
	* `--genomeSAindexNbases` must be reduced (for recomm., see STAR manual) 
* (defining the minimum and maximum __intron sizes__ that are allowed (the default setting for the maximum intron size is 1,000,000 bp)
* handling of genomes with __more than 5,000 scaffolds__ (usually reference genomes in a draft stage);
* using STAR for the detection of chimeric __(fusion) and circular tx__.

* __Result__: BAM file = compressed SAM file**

## For interaction with BAM/SAM files: samtools suite
	
 * __export__ samtools path (for convenience)
 
 	```
    $ export PATH=~/mat/software/samtools-1.2:$PATH
	
	# check out the header
    $ samtools view -H  WT_1_Aligned.sortedByCoord.out.bam
    
	# look at just the alignment section of the BAM file
    $ samtools view WT_1_Aligned.sortedByCoord.out.bam | head
   ```

### filter SAM files based on flags, optional fields

#### count unmapped reads: flag = 4

    samtools view -f 4 WT_1_Aligned.sortedByCoord.out.bam | wc -l
    
#### count mapped reads: flag != 4

	samtools view -F 4 WT_1_Aligned.sortedByCoord.out.bam | wc -l

#### MAPQ < 20

    samtools view -q 20 WT_1_Aligned.sortedByCoord.out.bam | head
    
    # to pipe into new file:
    samtools view -q 20 -h -b <BAM>


## Alignment QC

Let's do visual inspection with IGV!

## IGV

go to [https://www.broadinstitute.org/igv/IGV](https://www.broadinstitute.org/igv/IGV) --> Downloads --> email address --> go to archive --> download `2.2.13.zip` --> double click --> go to Downloads Folder in Terminal --> `sh igv.sh`  

"File" --> load from URL --> `chagall.med.cornell.edu/RNASEQcourse` --> BAM

* to __filter reads with insert size smaller than 1kb__ use CIGAR:
    * for the sake of simplicity, let's work on the SAM file:
    
            $ samtools view -h WT_1_Aligned.sortedByCoord.out.bam > WT_1_Aligned.sortedByCoord.out.sam

    * here's an example using grep, excluding lines with at least four digits followed by N

            $ egrep -v "[0-9][0-9][0-9][0-9]N" WT_1_Aligned.sortedByCoord.out.sam > smallInsert_reads.sam

    * same thing with awk, which can be used to match a regex within a specified column

            $ awk '!($6 ~ /[0-9][0-9][0-9][0-9]N/) {print $0}' WT_1_Aligned.sortedByCoord.out.sam > smallInsert_reads.sam


Homework: try to find optimal parameter settings for the yeast data set

