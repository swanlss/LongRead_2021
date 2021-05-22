## Bioinformatics Workshop 2021
### S13 Long-read Genomics, Metagenomics and Amplicon Analysis
### Tutorial/Practical

##### *Pierre Ramond, Swan LS Sow*
##### *Thursday, 27th May 2021*
<p>&nbsp;</p>

### 0. Background
Besides the ease, simplicity, speed and relatively lower cost of long-reads generated from 3rd generation sequencing tech, one of the key advantages of longer amplicons lies in its potential of increased taxonomic resolution that cannot be achieved by targeted 16S sub-region sequencing used in short-read sequencing platforms (Johnson et.al. 2019). Longer reads also allow significant improvements of genome assemblies (Koren and Phillippy, 2015).

The primary aim of this tutorial is to illustrate the benefits of longer amplicons by comparing the quality of taxonomic assignments of long versus short amplicons of the same sequences.
<p>&nbsp;</p>

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
<p>&nbsp;</p>

### 2. Extracting specific sub-regions of the 16S & 18S rRNA gene
The original reads generated from the MinION sequencing are ~1100 bp for the 16S amplicons and ~1200 bp for the 18S amplicons. We will use **cutadapt** to trim the sequences to the desired fragment lengths and extract specific 16S and 18S rRNA gene sub-regions. For example, to extract the 18S V4 region, we use the primer sequences that were developed by Stoeck et.al. (2010) as the adapter sequence parameter in **cutadapt** as follows:

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

</details><p>&nbsp;</p>

We'll then extract the long reads that came through the **cutadapt** pipeline based on the list of reads with **seqkit**'s *grep* function:

```
seqkit grep -f longread_wk2/18S_reads_ID.txt 18S.fastq -o longread_wk2/18S_og_reads.fastq
seqkit grep -f longread_wk2/16S_806R_reads_ID.txt 16S.fastq -o longread_wk2/16S_og_reads_806R.fastq
seqkit grep -f longread_wk2/16S_926R_reads_ID.txt 16S.fastq -o longread_wk2/16S_og_reads_926R.fastq
```

Finally, we'll remove the adapters, primers and Unique Molecular Identifiers (UMIs) from the long reads by trimming the first and last 80 bp of each sequence with **seqkit**'s *subseq* function:

```
seqkit subseq -r 80:-80  18S_og_reads.fastq > 18S_og_reads_trimm.fastq
seqkit subseq -r 80:-80  16S_og_reads_806R.fastq > 16S_og_reads_806R_trimm.fastq
seqkit subseq -r 80:-80  16S_og_reads_926R.fastq > 16S_og_reads_926R_trimm.fastq
```
<p>&nbsp;</p>

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
</details><p>&nbsp;</p>

We'll use **MOTHUR** to classify our sequences with the silva.nr_v138 database for the 16S sequences and pr2_version_4.13.0_18S for the 18S sequences. MOTHUR works better with sequences in the fasta format, so we'll first convert all fastq sequences to fasta format with **seqtk**:

```
cd sub_regions

seqtk seq -A 16S_og_reads_926R_trimm.fastq > `ls 16S_og_reads_926R_trimm.fastq | sed 's/\.fastq/\.fasta/'`
seqtk seq -A 16S_og_reads_806R_trimm.fastq > `ls 16S_og_reads_806R_trimm.fastq | sed 's/\.fastq/\.fasta/'`
seqtk seq -A 18S_og_reads_trimm.fastq > `ls 18S_og_reads_trimm.fastq | sed 's/\.fastq/\.fasta/'`
```

Classify the sequences with **MOTHUR** using the default *wang* method (Wang et.al., 2007). This should take ~ 5-10 minutes. The cutoff value indicates that only classified sequences with bootstrap values >= 80% will be returned.

```
# 16S
mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk2/sub_regions/16S_og_reads_806R_trimm.fasta, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"

mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk2/sub_regions/16S_og_reads_926R_trimm.fasta, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"

mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk2/sub_regions/16S_sub_V4_806R.fasta, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"

mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk2/sub_regions/16S_sub_V4_926R.fasta, reference=silva.nr_v138_1.align, taxonomy=silva.nr_v138_1.tax, cutoff=80)"
```
```
# 18S
mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk2/sub_regions/18S_og_reads_trimm.fasta, reference=pr2_version_4.13.0_18S_mothur.fasta,taxonomy=pr2_version_4.13.0_18S_mothur.tax, cutoff=80)"

mothur "#set.dir(input=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk/databases/);classify.seqs(fasta=/export/lv4/projects/NIOZ200/Data/Analysis_Bonito/6_UMI_BINNING/longread_wk2/sub_regions/18S_sub_V4_STOECK.fasta, reference=pr2_version_4.13.0_18S_mothur.fasta,taxonomy=pr2_version_4.13.0_18S_mothur.tax, cutoff=80)"
```
<p>&nbsp;</p>

### 4. Generating length gradient fragments from the trimmed long reads

To further test the theory, we'll generate sequence fragments from the original long reads with 100bp length variations. Make a new folder for the length gradient fragments, and then subfolders within that folder for the 16S and 18S fragments. Make a copy of the fasta file containing the trimmed long reads.

<details>
<summary>
<a class="btnfire small stroke"><em class="fas fa-chevron-circle-down"></em>&nbsp;&nbsp;Show me the code!</a>    
</summary>

```
mkdir Length_gradients
mkdir Length_gradients/18S
mkdir Length_gradients/16S

cp 18S_og_reads_trimm.fastq Length_gradients/18S/18S_trim_original.fastq
cp 16S_og_reads_806R_trimm.fastq Length_gradients/16S/16S_trim_original.fastq
```

</details>
