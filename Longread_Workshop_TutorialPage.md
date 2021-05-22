## Bioinformatics Workshop 2021 - S13 Long-read Genomics, Metagenomics and Amplicon Analysis
### Tutorial/Practical

##### *Pierre Ramond, Swan LS Sow*
##### *Thursday, 27th May 2021*

### 0. Background
Besides the ease, simplicity, speed and relatively lower cost of long-reads generated from 3rd generation sequencing tech, one of the key advantages of longer amplicons lies in its potential of increased taxonomic resolution that cannot be achieved by targeted 16S sub-region sequencing used in short-read sequencing platforms (Johnson et.al. 2019). Longer reads also allow significant improvements of genome assemblies (Koren and Phillippy, 2015).

The primary aim of this tutorial is to illustrate the benefits of longer amplicons by comparing the quality of taxonomic assignments of long versus short amplicons of the same sequences.

### 1. Dataset
We will work with long-read 16S and 18S amplicon dataset generated from samples taken from the Wadden Sea and Southern Ocean in this practical session. 16S amplicons were generated with the primer set A519F-1492R-pB-3771 (Martijn et.al., 2019), while the 18S amplicons were generated with the primer set Euk528F-U1391R (Edgcomb et.al. 2011).

First let's get some data into the project folder. Login to ADA, and *cd* to the following directory where all files required for this tutorial are located:

```
cd /export/lv4/projects/workshop_2021/S13_LongRead/
```

To be efficient with disk space, please make links from the sequence data fasta files to your working folders instead of copying them over.

```
ln -s /export/lv4/projects/workshop_2021/S13_LongRead/reads/ /export/lv3/scratch/workshop_2021/Users/<your_username>
```


### 2. Extracting specific sub-regions of the 16S & 18S rRNA gene
The original reads generated from the MinION sequencing are ~1100 bp for the 16S amplicons and ~1200 bp for the 18S amplicons. We will use *cutadapt* to trim the sequences to the desired fragment lengths and extract specific 16S and 18S rRNA gene sub-regions. For example, to extract the 18S V4 region, we use the primer sequences that were developed by Stoeck et.al. (2010) as the adapter sequence parameter in *cutadapt* as follows:

```
cutadapt -j 0 -e 0.3 -O 12 \
  --discard-untrimmed \
  -a CCAGCASCYGCGGTAATTCC...TYRATCAAGAACGAAAGT \
  -a ACTTTCGTTCTTGATYRA...GGAATTACCGCRGSTGCTGG \
  -M 600 \
  -o longread_wk/18S_sub_V4_STOECK.fasta \
18S.fastq
```

The 16S V4 region can be extracted by using the 515F-806R primer sequences developed by Caporaso et.al. (2011):

```
cutadapt -j 0 -e 0.3 -O 12 \
  --discard-untrimmed \
  -a GTGCCAGCMGCCGCGGTAA...ATTAGAWADDDBDGTAGTCC \
  -a GGACTACHVHHHTWTCTAAT...TTACCGCGGCKGCTGGCAC \
  -M 600 \
  -o longread_wk/16S_sub_V4_806R.fasta \
16S.fastq
```

Or alternatively, with the 515F-926R primer pair (Parada et.al., 2015), which targets the V4-5 region:

```
cutadapt -j 0 -e 0.3 -O 12 \
  --discard-untrimmed \
  -a GTGYCAGCMGCCGCGGTAA...AAACTYAAAKRAATTGRCGG \
  -a CCGYCAATTYMTTTRAGTTT...TTACCGCGGCKGCTGRCAC \
  -M 600 \
  -o longread_wk/16S_sub_V4_926R.fasta \
16S.fastq
```

Now, let's make a list of the reads that matched the adapter(primer) sequences from the *cutadapt* step above. You are encouraged to figure out how to do this based on unix commands covered in the previous workshop sessions, but for simplicity we have also provided the code below:

<details>
<summary>
<a class="btnfire small stroke"><em class="fas fa-chevron-circle-down"></em>&nbsp;&nbsp;Show me the code!</a>    
</summary>

```
grep ">" longread_wk2/18S_sub_V4_STOECK.fasta | sed 's/>//' | sed 's/\s.*$//' > longread_wk2/18S_reads_ID.txt
grep ">" longread_wk2/16S_sub_V4_806R.fasta | sed 's/>//' | sed 's/\s.*$//' > longread_wk2/16S_806R_reads_ID.txt
grep ">" longread_wk2/16S_sub_V4_926R.fasta | sed 's/>//' | sed 's/\s.*$//' > longread_wk2/16S_926R_reads_ID.txt
```

</details>


We'll then extract the long reads that came through the *cutadapt* pipeline based on the list of reads with *seqkit*'s *grep* function:

```
seqkit grep -f longread_wk2/18S_reads_ID.txt 18S.fastq -o longread_wk2/18S_og_reads.fastq
seqkit grep -f longread_wk2/16S_806R_reads_ID.txt 16S.fastq -o longread_wk2/16S_og_reads_806R.fastq
seqkit grep -f longread_wk2/16S_926R_reads_ID.txt 16S.fastq -o longread_wk2/16S_og_reads_926R.fastq
```

Finally, we'll remove the adapters, primers and Unique Molecular Identifiers (UMIs) from the long reads by trimming the first and last 80 bp of each sequence with *seqkit*'s *subseq* function:

```
seqkit subseq -r 80:-80  18S_og_reads.fastq > 18S_og_reads_trimm.fastq
seqkit subseq -r 80:-80  16S_og_reads_806R.fastq > 16S_og_reads_806R_trimm.fastq
seqkit subseq -r 80:-80  16S_og_reads_926R.fastq > 16S_og_reads_926R_trimm.fastq
```

### 3. Comparing taxonomy of the V4 sub-region fragments to the long reads

First, make a directory for the sub-region fragments and copy all the trimmed sequence fragments to this folder

<details>
<summary>
<a class="btnfire small stroke"><em class="fas fa-chevron-circle-down"></em>&nbsp;&nbsp;Show me the code!</a>    
</summary>

```
mkdir sub_regions

cp 16S_og_reads_806R_trimm.fastq sub_regions
cp 16S_og_reads_926R_trimm.fastq sub_regions
cp 16S_sub_V4_806R.fasta sub_regions
cp 16S_sub_V4_926R.fasta sub_regions
cp 18S_og_reads_trimm.fastq sub_regions
cp 18S_sub_V4_STOECK.fasta sub_regions
```

</details>

### X. Generating length gradients
