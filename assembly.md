# Assembling the House Wren Transcriptome using Trinity
This section focuses on generating the _de novo_ transcriptome assembly using the processed sequencing reads completed in the [previous section](./read_processing.md).  Afterwards, a series of cleaning, quality control, and assessment steps are performed prior to identifying SNPs.

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

## Step 2: Generate Super Transcripts
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
