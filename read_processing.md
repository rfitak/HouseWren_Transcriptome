# Cleaning the Raw Sequencing Reads
This section will start with the raw sequencing data and perform a series a cleaning steps to prepare the sequences for the transcriptome assembly.  The various steps include:
1.  Filtering low-quality reads, Trimming low-quality bases, adapter identification and removal
    - Program: [fastp](https://github.com/OpenGene/fastp) for paired-end reads   
2.  Removing reads that map conclusively to the American Shad mitochondrial genome
    - A mitogenome is already available, so we want to minimize their presence
    - Program [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for the mapping
6.  Kmer counting and Error-correcting the sequencing reads
    - Program: [musket v1.1](http://musket.sourceforge.net/homepage.htm)

### Raw Data Summary:

| Name | Type | Insert Size | # Paired Reads | # Bases | Q20 bases | Q30 Bases |
| --- | --- | --- | --- | --- | --- | --- |
| HOWR-1 | 151 bp; Paired-end | TBD | 313,465,270 | 94,666,511,540 | 88.2% | 81.4% |
| HOWR-2 | 151 bp; Mate-pair | 5-7 kb | 226,692,460 | 68,461,122,920 | 93.1% | 85.3% |
| __Total__ | n/a | n/a | 763,516,406 | 230,581,954,612 | 93.0% | 85.5% |


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
An example run is shown below, please see the scripts [fastp_PE500.sh](./Data/fastp_PE500.sh), [fastp_MP5k.sh](fastp_MP5k.sh), and [fastp_MP10k.sh](./Data/fastp_MP10k.sh) for more details on Job information.
```bash
# Assign names to each sample
name="HOWR-1"

# Trim PE reads
echo "Trimming the PE reads"
fastp \
   -i ${fwd1} \
   -I ${rev1} \
   -o ${name1}_F.trimmed.fq.gz \
   -O ${name1}_R.trimmed.fq.gz \
   --detect_adapter_for_pe \
   --cut_front \
   --cut_tail \
   --cut_window_size=4 \
   --cut_mean_quality=20 \
   --qualified_quality_phred=20 \
   --unqualified_percent_limit=30 \
   --n_base_limit=5 \
   --length_required=50 \
   --low_complexity_filter \
   --complexity_threshold=30 \
   --overrepresentation_analysis \
   --json=${name1}.json \
   --html=${name1}.html \
   --report_title="$name1" \
   --thread=8
```
_Parameters Explained:_
- --in1/--in2 :: input forward and reverse read files, recognizes gzip
- --out1/-out2 :: output forward and reverse read files, recognizes gzip
- --detect_adapter_for_pe :: enable PE adapter trimming
- --cut_front :: enable a 5' sliding window trimmer, like trimmomatic
- --cut_tail :: enable a 3' sliding window trimmer, like trimmomatic
- --cut_window_size=4 :: window size for the trimming
- --cut_mean_quality=20 :: mean base score across the window required, or else trim the last base
- --qualified_quality_phred=20 :: minimum base quality score to keep
- --unqualified_percent_limit=30 :: Percent of bases allowed to be less than q in a read
- --n_base_limit=5 :: if one read's number of N bases is >5, then this read pair is discarded
- --length_required=50 :: minimum read length to keep after trimming
- --low_complexity_filter :: filter sequences with a low complexity
- --complexity_threshold=30 :: threshold for sequence complexity filter
- --overrepresentation_analysis :: look for overrepresented sequences, like adapters
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
Detecting adapter sequence for read1...
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

Detecting adapter sequence for read2...
>Illumina TruSeq Adapter Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

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
total reads: 141847827
total bases: 19297745476
Q20 bases: 19017542785(98.548%)
Q30 bases: 18536747306(96.0565%)

Read2 aftering filtering:
total reads: 141847827
total bases: 19324744083
Q20 bases: 18958583251(98.1052%)
Q30 bases: 18346289679(94.9368%)

Filtering result:
reads passed filter: 283695654
reads failed due to low quality: 7153124
reads failed due to too many N: 20694
reads failed due to too short: 5971254
reads failed due to low complexity: 2654756
reads with adapter trimmed: 138880654
bases trimmed due to adapters: 3948422123
reads with polyX in 3' end: 9485728
bases trimmed in polyX tail: 806371712

Duplication rate: 15.4203%

Insert size peak (evaluated by paired-end reads): 150

JSON report: .json
HTML report: HOWR-2.html

fastp --in1 HOWR-2_R1_001.fastq.gz --in2 HOWR-2_R2_001.fastq.gz --out1 HOWR-2_cleaned.R1.fastq.gz --out2 HOWR-2_cleaned.R2.fastq.gz --detect_adapter_for_pe --adapter_fasta /home/rfitak/.bin/adapters.fa --cut_front --cut_tail --cut_window_size=4 --cut_mean_quality=20 --qualified_quality_phred=20 --average_qual 20 --unqualified_percent_limit=30 --n_base_limit=5 --length_required=50 --low_complexity_filter --complexity_threshold=30 --overrepresentation_analysis --trim_poly_x --poly_x_min_len 10 --html=HOWR-2.html --json=.json --report_title=HOWR-2 --thread=16 
fastp v0.20.0, time used: 3031 seconds
```
