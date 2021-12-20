# Assembling the House Wren Transcriptome using Trinity
This section focuses on generating the _de novo_ transcriptome assembly using the processed sequencing reads completed in the [previous section](./read_processing.md).  Afterwards, the assembled transcripts, which represent various isoforms, fragments, and chimeras of transcripts are combined/collapsed into Super Transcripts prior to identifying SNPs.

## Step 1:  Build transcriptome assembly
_Download Trinity singularity container_
```bash
wget https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/trinityrnaseq.v2.9.1.simg

# Load singularity
module load singularity

# Run trinity
singularity \
   exec \
   -e ~/PROGRAMS/trinityrnaseq.v2.9.1.simg \
   Trinity \
      --seqType fq \
      --SS_lib_type RF \
      --left `pwd`/HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.1.gz,`pwd`/HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.1.gz \
      --right `pwd`/HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.2.gz,`pwd`/HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.2.gz \
      --max_memory 512G \
      --CPU 32 \
      --output `pwd`/trinity1
```

_Get statistics from the final assembly_  
_All Trinity perl scripts from https://github.com/trinityrnaseq/trinityrnaseq_
```bash
# Get assembly stats
TrinityStats.pl Trinity.fasta
```

_Output_
```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	293691
Total trinity transcripts:	564671
Percent GC: 46.16

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 9873
	Contig N20: 6500
	Contig N30: 4685
	Contig N40: 3430
	Contig N50: 2457

	Median contig length: 462
	Average contig: 1082.80
	Total assembled bases: 611426756


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 7704
	Contig N20: 4642
	Contig N30: 2906
	Contig N40: 1842
	Contig N50: 1177

	Median contig length: 364
	Average contig: 725.45
	Total assembled bases: 213058014
```

## Step 2: Assembly assessment using Transrate
From the [Transrate website](http://hibberdlab.com/transrate/):  
"Transrate is software for _de-novo_ transcriptome assembly quality analysis. It examines your assembly in detail and compares it to experimental evidence such as the sequencing reads, reporting quality scores for contigs and assemblies. This allows you to choose between assemblers and parameters, filter out the bad contigs from an assembly, and help decide when to stop trying to improve the assembly"

_Run transrate v1.0.3_
```bash
transrate \
   --assembly=Trinity.fasta \
   --left=HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.1.gz,HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.1.gz \
   --right=HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.2.gz,HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.2.gz \
   --threads=32 \
   --output=. \
   --loglevel=info
```

_Output_
```
[ INFO] 2020-06-24 08:32:14 : Contig metrics:
[ INFO] 2020-06-24 08:32:14 : -----------------------------------
[ INFO] 2020-06-24 08:32:14 : n seqs                       564671
[ INFO] 2020-06-24 08:32:14 : smallest                        177
[ INFO] 2020-06-24 08:32:14 : largest                       44898
[ INFO] 2020-06-24 08:32:14 : n bases                   611426756
[ INFO] 2020-06-24 08:32:14 : mean len                    1082.76
[ INFO] 2020-06-24 08:32:14 : n under 200                     118
[ INFO] 2020-06-24 08:32:14 : n over 1k                    140335
[ INFO] 2020-06-24 08:32:14 : n over 10k                     4349
[ INFO] 2020-06-24 08:32:14 : n with orf                   110569
[ INFO] 2020-06-24 08:32:14 : mean orf percent              35.54
[ INFO] 2020-06-24 08:32:14 : n90                             376
[ INFO] 2020-06-24 08:32:14 : n70                            1032
[ INFO] 2020-06-24 08:32:14 : n50                            2457
[ INFO] 2020-06-24 08:32:14 : n30                            4686
[ INFO] 2020-06-24 08:32:14 : n10                            9877
[ INFO] 2020-06-24 08:32:14 : gc                             0.46
[ INFO] 2020-06-24 08:32:14 : bases n                           0
[ INFO] 2020-06-24 08:32:14 : proportion n                    0.0
[ INFO] 2020-06-24 08:32:14 : Contig metrics done in 110 seconds
[ INFO] 2020-06-24 08:32:14 : Calculating read diagnostics...
[ INFO] 2020-06-24 13:43:05 : Read mapping metrics:
[ INFO] 2020-06-24 13:43:05 : -----------------------------------
[ INFO] 2020-06-24 13:43:05 : fragments                 426002474
[ INFO] 2020-06-24 13:43:05 : fragments mapped           42544746
[ INFO] 2020-06-24 13:43:05 : p fragments mapped              0.1
[ INFO] 2020-06-24 13:43:05 : good mappings              26280824
[ INFO] 2020-06-24 13:43:05 : p good mapping                 0.06
[ INFO] 2020-06-24 13:43:05 : bad mappings               16263922
[ INFO] 2020-06-24 13:43:05 : potential bridges                 0
[ INFO] 2020-06-24 13:43:05 : bases uncovered           445334069
[ INFO] 2020-06-24 13:43:05 : p bases uncovered              0.73
[ INFO] 2020-06-24 13:43:05 : contigs uncovbase            420242
[ INFO] 2020-06-24 13:43:05 : p contigs uncovbase            0.74
[ INFO] 2020-06-24 13:43:05 : contigs uncovered            564671
[ INFO] 2020-06-24 13:43:05 : p contigs uncovered             1.0
[ INFO] 2020-06-24 13:43:05 : contigs lowcovered           564671
[ INFO] 2020-06-24 13:43:05 : p contigs lowcovered            1.0
[ INFO] 2020-06-24 13:43:05 : contigs segmented             40990
[ INFO] 2020-06-24 13:43:05 : p contigs segmented            0.07
[ INFO] 2020-06-24 13:43:05 : Read metrics done in 18651 seconds
[ INFO] 2020-06-24 13:43:05 : No reference provided, skipping comparative diagnostics
[ INFO] 2020-06-24 13:44:18 : TRANSRATE ASSEMBLY SCORE     0.0047
[ INFO] 2020-06-24 13:44:18 : -----------------------------------
[ INFO] 2020-06-24 13:44:18 : TRANSRATE OPTIMAL SCORE      0.0177
[ INFO] 2020-06-24 13:44:18 : TRANSRATE OPTIMAL CUTOFF     0.1739
[ INFO] 2020-06-24 13:44:19 : good contigs                 259399
[ INFO] 2020-06-24 13:44:19 : p good contigs                 0.46
[ INFO] 2020-06-24 13:44:19 : Writing contig metrics for each contig to /home/fitaklab/data/HOUSE_WREN/TRANSRATE/Trinity/contigs.csv
[ INFO] 2020-06-24 13:45:57 : Writing analysis results to assemblies.csv
```


## Step 3: Generate Super Transcripts
"SuperTranscripts provide a gene-like view of the transcriptional complexity of a gene. SuperTranscripts were originally defined by Nadia Davidson, Anthony Hawkins, and Alicia Oshlack as described in their publication [SuperTranscripts: a data driven reference for analysis and visualisation of transcriptome. Genome Biology, 2017](https://doi.org/10.1186/s13059-017-1284-1). SuperTranscripts are useful in the context of genome-free _de novo_ transcriptome assembly in that they provide a genome-like reference for studying aspects of the gene including differential transcript usage (aka. differential exon usage) and as a substrate for mapping reads and identifying allelic polymorphisms.

A SuperTranscript is constructed by collapsing unique and common sequence regions among splicing isoforms into a single linear sequence."  
From: https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts

_Generate Super Transcripts file_
```bash
Trinity_gene_splice_modeler.py \
   --trinity_fasta Trinity.fasta \
   --out_prefix Trinity.SuperTrans \
   --incl_malign \
   --verbose
```

_Summary of the Super Transcripts_
```bash
# Get assembly stats
TrinityStats.pl Trinity.SuperTrans.fasta
```

_Output_
```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	293691
Total trinity transcripts:	293691
Percent GC: 45.01

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 9738
	Contig N20: 5680
	Contig N30: 3496
	Contig N40: 2177
	Contig N50: 1408

	Median contig length: 379
	Average contig: 798.98
	Total assembled bases: 234652073


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 9738
	Contig N20: 5680
	Contig N30: 3496
	Contig N40: 2177
	Contig N50: 1408

	Median contig length: 379
	Average contig: 798.98
	Total assembled bases: 234652073
```

## Summary Table
