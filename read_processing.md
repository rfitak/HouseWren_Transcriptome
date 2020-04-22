# Cleaning the Raw Sequencing Reads
This section will start with the raw sequencing data and perform a series a cleaning steps to prepare the sequences for the transcriptome assembly.  The various steps include:
1.  Filtering low-quality reads, trimming low-quality bases, adapter identification and removal
    - Program: [fastp](https://github.com/OpenGene/fastp) for paired-end reads 
2.  Correcting for sequencing errors
    - Program: [rcorrector](https://github.com/mourisl/Rcorrector)
3.  Removing read pairs with an unfixable error.
    - Program: [FilterUncorrectabledPEfastq.py](https://github.com/harvardinformatics/TranscriptomeAssemblyTools).
4.  Remove reads originating from rRNA.
    - Program:  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) alignment to the [SILVA](https://www.arb-silva.de) rRNA database.

### Raw Data Summary:

| Name | Type | # Paired Reads | # Bases | Q20 bases | Q30 Bases |
| --- | --- | --- | --- | --- | --- |
| HOWR-1 | 150 bp; Paired-end | 350,315,944 | 105,094,783,200 | 97.6% | 94.0% |
| HOWR-2 | 150 bp; Paired-end | 149,747,741 | 44,924,322,300 | 97.1% | 93.6% |
| __Total__ | 150 bp; Paired-end | 500,063,685 | 150,019,105,500 | 97.4% | 93.9% |


## Step 1: Read trimming and filtering
Here the new software [fastp v0.20.0](https://github.com/OpenGene/fastp) was used to trim all reads. It combines a QC (similar to FastQC) along with various trimming and filtering functions. The publication can be found here:  
Chen S, Zhou Y, Chen Y, Gu (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. _Bioinformatics_ 34(17):i884–i890. https://doi.org/10.1093/bioinformatics/bty560

_Installation:_
```bash
# Install fastp using git
git clone https://github.com/OpenGene/fastp.git
cd fastp
make
```

_Run fastp_  
Only sample HOWR-1 is shown, the exact same code was repeated for HOWR-2 (`name="HOWR-2"`).
```bash
# Assign names to each sample
name="HOWR-1"

# Run fastp cleaning
fastp \
   --in1=${name}_R1_001.fastq.gz \
   --in2=${name}_R2_001.fastq.gz \
   --out1=${name}_cleaned.R1.fastq.gz \
   --out2=${name}_cleaned.R2.fastq.gz \
   --adapter_fasta=adapters.fa \
   --cut_front \
   --cut_tail \
   --cut_window_size=4 \
   --cut_mean_quality=20 \
   --qualified_quality_phred=20 \
   --average_qual=20 \
   --unqualified_percent_limit=30 \
   --n_base_limit=5 \
   --length_required=50 \
   --low_complexity_filter \
   --complexity_threshold=30 \
   --overrepresentation_analysis \
   --trim_poly_x \
   --poly_x_min_len=10 \
   --html=${name}.html \
   --json=${nname}.json \
   --report_title="$name" \
   --thread=16
```
_Parameters Explained:_
- --in1/--in2 :: input forward and reverse read files, recognizes gzip
- --out1/-out2 :: output forward and reverse read files, recognizes gzip
- --adapter_fasta :: a file of known Illumina adapters to trim, see [adapters.fa](./Data/adapters.fa)
- --cut_front :: enable a 5' sliding window trimmer, like trimmomatic
- --cut_tail :: enable a 3' sliding window trimmer, like trimmomatic
- --cut_window_size=4 :: window size for the trimming
- --cut_mean_quality=20 :: mean base score across the window required, or else trim the last base
- --qualified_quality_phred=20 :: minimum base quality score to keep
- --average_qual=20 :: remove read of the average quality across all bases is < 20
- --unqualified_percent_limit=30 :: Percent of bases allowed to be less than q in a read
- --n_base_limit=5 :: if one read's number of N bases is >5, then this read pair is discarded
- --length_required=50 :: minimum read length to keep after trimming
- --low_complexity_filter :: filter sequences with a low complexity
- --complexity_threshold=30 :: threshold for sequence complexity filter
- --overrepresentation_analysis :: look for overrepresented sequences, like adapters
- --trim_poly_x :: trim strings of homopolymers at the 3' end of reads
- --poly_x_min_len 10 :: minimum length of homopolymer ot trim
- --json=${name}.json :: output file name, JSON format
- --html=${name}.html :: output file name, HTML format
- --report_title="$name" :: output report tile
- --thread=16 :: number of cpus to use

_See the Output HTML/PDF Files from ```fastp``` below:_
- [HOWR-1](./Data/HOWR-1.pdf)
- [HOWR-2](./Data/HOWR-2.pdf)

### Output Summary
_HOWR-1 Reads_
```
# HOWR-1 PE reads
Read1 before filtering:
total reads: 350315944
total bases: 52547391600
Q20 bases: 51683521994(98.356%)
Q30 bases: 50218703328(95.5684%)

Read2 before filtering:
total reads: 350315944
total bases: 52547391600
Q20 bases: 50838173238(96.7473%)
Q30 bases: 48611492675(92.5098%)

Read1 after filtering:
total reads: 340474459
total bases: 47631199046
Q20 bases: 47085051064(98.8534%)
Q30 bases: 45950789198(96.472%)

Read2 aftering filtering:
total reads: 340474459
total bases: 47569274139
Q20 bases: 46572539164(97.9047%)
Q30 bases: 44820344133(94.2212%)

Filtering result:
reads passed filter: 680948918
reads failed due to low quality: 15303720
reads failed due to too many N: 47110
reads failed due to too short: 3058662
reads failed due to low complexity: 1273478
reads with adapter trimmed: 261704514
bases trimmed due to adapters: 7048775313
reads with polyX in 3' end: 3158499
bases trimmed in polyX tail: 132541574

Duplication rate: 12.3609%

Insert size peak (evaluated by paired-end reads): 150

JSON report: HOWR-1.json
HTML report: HOWR-1.html

fastp v0.20.0, time used: 7621 seconds
```

_HOWR-2 Reads_
```
# HOWR-2 PE reads
Read1 before filtering:
total reads: 149747741
total bases: 22462161150
Q20 bases: 21896887478(97.4834%)
Q30 bases: 21187597776(94.3257%)

Read2 before filtering:
total reads: 149747741
total bases: 22462161150
Q20 bases: 21719331733(96.693%)
Q30 bases: 20860532350(92.8697%)

Read1 after filtering:
total reads: 141658829
total bases: 19271163888
Q20 bases: 18992644933(98.5547%)
Q30 bases: 18513125388(96.0665%)

Read2 aftering filtering:
total reads: 141658829
total bases: 19283102833
Q20 bases: 18921388552(98.1242%)
Q30 bases: 18312653502(94.9674%)

Filtering result:
reads passed filter: 283317658
reads failed due to low quality: 7017686
reads failed due to too many N: 20682
reads failed due to too short: 6737414
reads failed due to low complexity: 2402042
reads with adapter trimmed: 144401586
bases trimmed due to adapters: 4514773965
reads with polyX in 3' end: 4433478
bases trimmed in polyX tail: 344676081

Duplication rate: 15.4203%

Insert size peak (evaluated by paired-end reads): 150

JSON report: HOWR-2.json
HTML report: HOWR-2.html
 
fastp v0.20.0, time used: 3220 seconds
```

### After read cleaning:

| Name | Type | # Paired Reads | # Bases | Q20 bases | Q30 Bases |
| --- | --- | --- | --- | --- | --- |
| HOWR-1 | Raw sequences | 350,315,944 | 105,094,783,200 | 97.6% | 94.0% |
| HOWR-1 | Cleaned sequences | 340,474,459 | 95,200,473,185 | 98.4% | 95.3% |
| HOWR-2 | Raw sequences | 149,747,741 | 44,924,322,300 | 97.1% | 93.6% |
| HOWR-2 | Cleaned sequences | 141,658,829 | 38,554,266,721 | 98.3% | 95.5% |
| __Total__ | Cleaned sequences | 482,133,288 | 133,754,739,906 | 98.4% | 95.4% |

<br>

## Step 2: Correcting for sequencing errors
[Rcorrector v1.0.4 ce5d06b](https://github.com/mourisl/Rcorrector) is a kmer-based, error-correction method for RNA-seq data. The error-correction context differs substantially in RNA-seq data from whole-genome sequencing (WGS) data, thus requiring specific correction procedures.  Methods like _rcorrector_ can significantly improve downstream RNA-seq data analysis, especially if used to produce a _de novo_ transcriptome assembly.  The publication describing _rcorrector_ can be found here:  
Song L, Florea L (2015) Rcorrector: Efficient and accurate error correction for Illumina RNA-seq reads. _GigaScience_ 4:48. https://doi.org/10.1186/s13742-015-0089-y

_Installation:_
```bash
# Install Rcorrector using git
git clone https://github.com/mourisl/rcorrector.git
cd rcorrector/
make
```

_Run rcorrector_
```bash
run_rcorrector.pl \
   -1 HOWR-1_cleaned.R1.fastq.gz,HOWR-2_cleaned.R1.fastq.gz \
   -2 HOWR-1_cleaned.R2.fastq.gz,HOWR-2_cleaned.R2.fastq.gz \
   -k 23 \
   -t 16
```
_Parameters Explained:_
- -1 :: comma separated list of forward read files (PE reads), recognizes gzip
- -2 :: comma separated list of reverse read files (PE reads), recognizes gzip
- -k :: k-mer length, default=23
- -t :: number of threads

_Output from Rcorrector_
```bash
Stored 713812866 kmers
Weak kmer threshold rate: 0.027540 (estimated from 0.950/1 of the chosen kmers)
Processed 964266576 reads
	Corrected 319704705 bases.
```
This means 319704705 / 133754739906 = __0.24%__ of bases were corrected

## Step 3: Remove read pairs with an unfixable error.
The step above uses `rcorrector` to identify erroneous k-mers, then correct them if possible.  However, sometimes `rcorrector` cannot correct ann erroneous k-mer, because it cannot identify a change that doesn't result in another rare k-mer.  As a result, with so much data, it is simply better to just eliminate these reads.  This saves not only disk space, but will speed up computational time and result in fewer assembly errors.  We will use the python script `FilterUncorrectabledPEfastq.py` which can be obtained from Harvard Informatics' [TranscriptomeAssemblyTools](https://github.com/harvardinformatics/TranscriptomeAssemblyTools).

_Run FilterUncorrectabledPEfastq.py_
```bash
# Remove unfixable errors from HOWR-1
FilterUncorrectabledPEfastq.py \
   -1 HOWR-1_cleaned.R1.cor.fq.gz \
   -2 HOWR-1_cleaned.R2.cor.fq.gz  \
   -s HOWR-1 

# Remove unfixable errors from HOWR-2
FilterUncorrectabledPEfastq.py \
   -1 HOWR-2_cleaned.R1.cor.fq.gz \
   -2 HOWR-2_cleaned.R2.cor.fq.gz \
   -s HOWR-2

# Change all output read file names and compress them
gzip unfixrm_HOWR-1_cleaned.R1.cor.fq
mv unfixrm_HOWR-1_cleaned.R1.cor.fq.gz HOWR-1_cleaned.R1.cor.unfixrm.fq.gz
gzip unfixrm_HOWR-1_cleaned.R2.cor.fq
mv unfixrm_HOWR-1_cleaned.R2.cor.fq.gz HOWR-1_cleaned.R2.cor.unfixrm.fq.gz
gzip unfixrm_HOWR-2_cleaned.R1.cor.fq
mv unfixrm_HOWR-2_cleaned.R1.cor.fq.gz HOWR-2_cleaned.R1.cor.unfixrm.fq.gz
gzip unfixrm_HOWR-2_cleaned.R2.cor.fq
mv unfixrm_HOWR-2_cleaned.R2.cor.fq.gz HOWR-2_cleaned.R2.cor.unfixrm.fq.gz
```

_Output from FilterUncorrectabledPEfastq.py_
```bash
# Output file: rmunfixable_HOWR-1.log
total PE reads:340474459
removed PE reads:6605775
retained PE reads:333868684
R1 corrected:75851187
R2 corrected:68220420
pairs corrected:103973135
R1 unfixable:2036835
R2 unfixable:2646289
both reads unfixable:1922651

# Output file: rmunfixable_HOWR-2.log
total PE reads:141658829
removed PE reads:3757435
retained PE reads:137901394
R1 corrected:32881899
R2 corrected:29190921
pairs corrected:43493718
R1 unfixable:1298352
R2 unfixable:1360213
both reads unfixable:1098870
```
| Name | Type | # Paired Reads | # Bases |
| --- | --- | --- | --- |
| HOWR-1 | Raw sequences | 350,315,944 | 105,094,783,200 |
| HOWR-1 | Cleaned sequences | 340,474,459 | 95,200,473,185 |
| HOWR-1 | R1 corrected | 75,851,187 | NA |
| HOWR-1 | R2 corrected | 68,220,420 | NA |
| HOWR-1 | Pairs corrected | 103,973,135 | NA |
| HOWR-1 | Retained | 333,868,684 | 93,343,653,959 |
| HOWR-2 | Raw sequences | 149,747,741 | 44,924,322,300 |
| HOWR-2 | Cleaned sequences | 141,658,829 | 38,554,266,721 |
| HOWR-2 | R1 corrected | 32,881,899 | NA |
| HOWR-2 | R2 corrected | 29,190,921 | NA |
| HOWR-2 | Pairs corrected | 43,493,718 | NA |
| HOWR-2 | Retained | 137,901,394 | 37,542,362,474 |
| __Total__ | Cleaned + Corrected + Retained | 471,770,078 | 130,886,016,433 |

## Step 4: Remove rRNA contamination
The library preparation protocol we used ([NEBNext® Ultra™ II Directional RNA Library Prep Kit for Illumina®](https://www.neb.com/products/e7760-nebnext-ultra-ii-directional-rna-library-prep-kit-for-illumina#Product%20Information)) is based off a dUTP protocol, which is essentially a poly(A) RNA enrinchment with strand specificity.  Depsite enriching for mRNA, numerous rRNA sequences easily make their way into the library, at times at high abundance.  Here, we will screen for rRNA sequences by comparing (mapping) our transcriptome read pairs to the [SILVA](https://www.arb-silva.de) database.  Reads that match rRNA sequences will be removed, since we do not want them as part of our transcriptome and subsequent SNP calling.

_Download and prepare the SILVA database_
```bash
# Download the SILVA SSU Parc database (release 138)
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138_SSUParc_tax_silva.fasta.gz

# Download the SILVA LSU Parc database (release 132)
wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_LSUParc_tax_silva.fasta.gz

# Concatenate the two together, change uracil (U) to thymine (T)
zcat \
   SILVA_138_SSUParc_tax_silva.fasta.gz \
   SILVA_132_LSUParc_tax_silva.fasta.gz | \
   sed '/^[^>]/ y/uU/tT/' > SILVA.fa

# Remove files not needed
rm SILVA_138_SSUParc_tax_silva.fasta.gz SILVA_132_LSUParc_tax_silva.fasta.gz
```

_Index the database and map reads - only retaining the unmapped reads_
```bash
# Index SILVA reference
bowtie2-build --threads 16 SILVA.fa SILVA.fa

# Map HOWR-1
bowtie2 \
   --very-sensitive-local \
   --phred33 \
   -x SILVA.fa \
   -1 HOWR-1_cleaned.R1.cor.unfixrm.fq.gz \
   -2 HOWR-1_cleaned.R2.cor.unfixrm.fq.gz \
   --nofw \
   --threads 16 \
   --un-conc-gz HOWR-1_cleaned.R1.cor.unfixrm.rmrRNA.fq.gz

# Map HOWR-2
bowtie2 \
   --very-sensitive-local \
   --phred33 \
   -x SILVA.fa \
   -1 HOWR-2_cleaned.R1.cor.unfixrm.fq.gz \
   -2 HOWR-2_cleaned.R2.cor.unfixrm.fq.gz \
   --nofw \
   --threads 16 \
   --un-conc-gz HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.gz
```
_Parameters explained_
- -very-sensitive-local :: a present for bowtie2 to perform sensitive, local mapping of reads
- --phred33 :: standard quality score encoding of fastq files (default)
- -x SILVA.fa :: reference genome prefix
- -1 :: forward read file
- -2 :: reverse read file
- -nofw :: reads are in the 'RF' orientation as opposed to 'FR', since we used the dUTP stranded protocol
- --threads 16 :: use 16 cpus
- --un-conc-gz :: store unaligned, properly paired reads in this gzipped output file
