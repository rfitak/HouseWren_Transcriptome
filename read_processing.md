# Cleaning the Raw Sequencing Reads
This section will start with the raw sequencing data and perform a series a cleaning steps to prepare the sequences for the transcriptome assembly.  The various steps include:
1.  Filtering low-quality reads, trimming low-quality bases, adapter identification and removal
    - Program: [fastp](https://github.com/OpenGene/fastp) for paired-end reads 
2.  Correcting for sequencing errors
    - Program: [rcorrector](https://github.com/mourisl/Rcorrector)
3.  Removing reads that map conclusively to the mitochondrial genome

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
Detecting adapter sequence for read1...
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

Detecting adapter sequence for read2...
No adapter detected for read2

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
total reads: 340602122
total bases: 47652261110
Q20 bases: 47104879440(98.8513%)
Q30 bases: 45969609495(96.4689%)

Read2 aftering filtering:
total reads: 340602122
total bases: 47601307266
Q20 bases: 46601380520(97.8994%)
Q30 bases: 44846464372(94.2127%)

Filtering result:
reads passed filter: 681204244
reads failed due to low quality: 15370542
reads failed due to too many N: 47112
reads failed due to too short: 2630734
reads failed due to low complexity: 1379256
reads with adapter trimmed: 258405625
bases trimmed due to adapters: 6837694734
reads with polyX in 3' end: 5836677
bases trimmed in polyX tail: 279591204

Duplication rate: 12.3608%

Insert size peak (evaluated by paired-end reads): 150

JSON report: .json
HTML report: HOWR-1.html

fastp --in1 HOWR-1_R1_001.fastq.gz --in2 HOWR-1_R2_001.fastq.gz --out1 HOWR-1_cleaned.R1.fastq.gz --out2 HOWR-1_cleaned.R2.fastq.gz --detect_adapter_for_pe --adapter_fasta /home/rfitak/.bin/adapters.fa --cut_front --cut_tail --cut_window_size=4 --cut_mean_quality=20 --qualified_quality_phred=20 --average_qual 20 --unqualified_percent_limit=30 --n_base_limit=5 --length_required=50 --low_complexity_filter --complexity_threshold=30 --overrepresentation_analysis --trim_poly_x --poly_x_min_len 10 --html=HOWR-1.html --json=.json --report_title=HOWR-1 --thread=16 
fastp v0.20.0, time used: 7718 seconds
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

JSON report: .json
HTML report: HOWR-2.html

fastp --in1=HOWR-2_R1_001.fastq.gz --in2=HOWR-2_R2_001.fastq.gz --out1=HOWR-2_cleaned.R1.fastq.gz --out2=HOWR-2_cleaned.R2.fastq.gz --adapter_fasta=/home/rfitak/.bin/adapters.fa --cut_front --cut_tail --cut_window_size=4 --cut_mean_quality=20 --qualified_quality_phred=20 --average_qual=20 --unqualified_percent_limit=30 --n_base_limit=5 --length_required=50 --low_complexity_filter --complexity_threshold=30 --overrepresentation_analysis --trim_poly_x --poly_x_min_len=10 --html=HOWR-2.html --json=.json --report_title=HOWR-2 --thread=16 
fastp v0.20.0, time used: 3220 seconds
```

## Step 2: Correcting for sequencing errors
[Rcorrector vce5d06b](https://github.com/mourisl/Rcorrector) is a kmer-based, error-correction method for RNA-seq data. The error-correction context differs substantially in RNA-seq data from whole-genome sequencing (WGS) data, thus requiring specific correction procedures.  Methods like _rcorrector_ can significantly improve downstream RNA-seq data analysis, especially if used to produce a _de novo_ transcriptome assembly.  The publication describing _rcorrector_ can be found here:
Song L, Florea L (2015) Rcorrector: Efficient and accurate error correction for Illumina RNA-seq reads. _GigaScience_ 4:48. https://doi.org/10.1186/s13742-015-0089-y

_Installation:_
```bash
# Install fastp using git
git clone https://github.com/mourisl/rcorrector.git
cd rcorrector/
make
```
