# Identifying SNPs
In this section we first identify SNPs from the reads aligned to the Super Transcripts reference fasta file. the full pipeline is described here https://github.com/trinityrnaseq/trinityrnaseq/wiki/Variant-Calling. This pipeline is well established, and other work has actually found that the Trinity + STAR + GATK pipeline produces the best variant calls (e.g., see [Zhao et al. 2019](https://doi.org/10.1186/s12864-019-5533-4)). After calling SNPs, we perform a series of steps in order to generate a set of high-quality candidate SNPs and flanking regions to be submitted to the Sequenom genotyping facility.

## Step 1:  Running Trinity's built-in GATK pipeline on the Super Transcripts
See the full pipeline here: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Variant-Calling

Dependencies that need to be installed:  
  - [Picard v2.23.0](https://broadinstitute.github.io/picard/)
    - `wget https://github.com/broadinstitute/picard/releases/download/2.23.0/picard.jar`
  - [STAR v2.7.3a](https://github.com/alexdobin/STAR)
  - [samtools v1.9](https://samtools.github.io)
  - [GATK v4.1.7.0](https://gatk.broadinstitute.org/hc/en-us)
    - `wget https://github.com/broadinstitute/gatk/releases/download/4.1.7.0/gatk-4.1.7.0.zip`
    - `unzip gatk-4.1.7.0.zip`

_Setup Picard and GATK environmental variables_
```bash
export PICARD_HOME=/home/rfitak/.bin/
export GATK_HOME=/home/rfitak/PROGRAMS/gatk-4.1.7.0
```

_Run SNP-calling Pipeline_
```bash
# Merge read files together
cat \
   HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.1.gz \
   HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.1.gz > left.fq.gz
cat \
   HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.2.gz \
   HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.2.gz > right.fq.gz

# Run Trinity SNP Pipeline
run_variant_calling.py \
   --st_fa ./trinity1/Trinity.SuperTrans.fasta \
   --st_gtf ./trinity1/Trinity.SuperTrans.gtf \
   -p left.fq.gz right.fq.gz \
   -o ./SNPS \
   -l 150 \
   -t 32 \
   -m 512000000000

rm -rf left.fq.gz right.fq.gz
```

**See various output and log files from the code above:**
- [snp-calling.err](./Data/snp-calling.err)
- [variant_calling.log](./Data/variant_calling.log)
- [Log.out](./Data/Log.out)
- [Log.final.out](./Data/Log.final.out)

## SNP output summary
The final set of called SNPs are in the output file `filtered_output.vcf`.  In total, **2,200,840** variants were identified (including SNPs, indels, etc).  However, not all these SNPs are of high quality.  The pipeline above does filter SNPs using GATK's `VariantFiltration` engine and the following criteria:
- `-window 35`
  - The window size (in bases) in which to evaluate clustered SNPs
- `-cluster 3`
  - The number of SNPs which make up a cluster.
- `--filter-name FS -filter "FS > 30.0"`
  - Strand bias estimated using Fisher's Exact Test. Strand bias is a type of sequencing bias in which one DNA strand is favored over the other, which can result in incorrect evaluation of the amount of evidence observed for one allele vs. the other. The FisherStrand (FS) annotation is one of several methods that aims to evaluate whether there is strand bias in the data. It uses Fisher's Exact Test to determine if there is strand bias between forward and reverse strands for the reference or alternate allele. The output is a Phred-scaled p-value. The higher the output value, the more likely there is to be bias. More bias is indicative of false positive calls.
  - FS > 30 means a P-value for bias < 0.001
- `--filter-name QD -filter "QD < 2.0"`
  - The QD is the QUAL score normalized by allele depth (AD) for a variant
  - This annotation puts the variant confidence QUAL score into perspective by normalizing for the amount of coverage available. Because each read contributes a little to the QUAL score, variants in regions with deep coverage can have artificially inflated QUAL scores, giving the impression that the call is supported by more evidence than it really is. To compensate for this, we normalize the variant confidence by depth, which gives us a more objective picture of how well supported the call is.

All variants in `filtered_output.vcf` are labeled as either `PASS` (passing all the above criteria), or has one or more of the labels `QD`, `SnpCluster`, and `FS` when they don't pass that specific filter above. Below is a summary of the variants based on the filtering criteria.  In total, **1,299,784** (59.1%) variants passed the filter.
```bash
# Find the number of SNPs matching each filter label
grep -v "^#" filtered_output.vcf | cut -f7 | sort | uniq -c

# Output
   2340 FS
    603 FS;QD
    308 FS;QD;SnpCluster
   1674 FS;SnpCluster
1299784 PASS
 200579 QD
  93293 QD;SnpCluster
 602259 SnpCluster
```

_Here are some additional summary stats and plots for "PASS" SNPs_
```bash
# Get summary stats using bcftools
bcftools stats -f 'PASS' filtered_output.vcf > out.ck

# Plot results
plot-vcfstats -p SNP-STATS out.ck
```
![summary.png](./images/summary.png)

|![depth.0.png](./images/depth.0.png)|![indels.0.png](./images/indels.0.png)|
| --- | --- |
|![tstv_by_qual.0.png](./images/tstv_by_qual.0.png)|![substitutions.0.png](./images/substitutions.0.png)|

## Step 2:  Additional SNP filtering
This section passes all the variants through additional checks to identify a final set of robust SNPs to submit for Sequenom genotyping.
1. Only retain biallelic SNPs with a heterozygous genotype that PASS the initial filter and have adequate depth of coverage
2. Generate 200 bp of flanking sequence for each SNP
    - remove SNPs with <80 bp flanking sequence on either side
    - these SNPs are poor candidates for Sequenom primer development
3. Remove any Super Transcripts (i.e. Genes) with a TPM count <1
    - after mapping reads back to the Super Transcripts, sequences with littel read coverage lack support
4. BLAST to a reference bird genome.
    - remove any loci that cross splice junctions
5. Remove any loci that overlap repeat masked region

#### Step 2.1: Only retain biallelic SNPs with a heterozygous genotype that PASS the initial filter and have adequate depth of coverage

_**2.1.1:** Label Homozygous sites or too little coverage in the filter column_

```bash
# Add filters "HOM" and "LowCov"
bcftools filter -s "HOM" -e 'COUNT(GT="hom")>0' filtered_output.vcf | \
   bcftools filter -s "LowCov" -e 'FORMAT/AD[*:*]<2' -O z > filtered_output2.vcf.gz
   # still 2200840 variants in the file

# Get ONLY the SNPs (exclude MNPs, indels, etc) with PASS filter:
bcftools view --min-alleles 2 --max-alleles 2 --types snps -f PASS filtered_output2.vcf.gz -O z > clean.vcf.gz
tabix -p vcf clean.vcf.gz
   # 1,018,044 SNPs retained
```

#### Step 2.2: Generate 200 bp of flanking sequence and remove SNPs with inadequate flanking sequence

_**2.2.1:** Get flanking sequence using `vcfprimers` in [vcflib](https://github.com/vcflib/vcflib) package_

```bash
module load anaconda3/97
source activate samtools # vcflib installed here via conda

# Get left and right flanking sequences
vcfprimers -l 200 -f ../trinity1/Trinity.SuperTrans.fasta clean.vcf.gz > flanks.tmp.fa

# Add flanking sequences
   # fill-fs -l 100 -r ../trinity1/Trinity.SuperTrans.fasta clean.vcf.gz > clean.flank.vcf  # failed, cant figure out why it wants to take negative values
```

_**2.2.2:** Convert flanking sequence file into a clean, fasta file_

```bash
# Convert to tab format, then merge
paste <(fasta2tab flanks.tmp.fa | sed -n 1~2p) <(fasta2tab flanks.tmp.fa | sed -n 2~2p) > tmp1

# Get SNP alleles
grep -v "^#" <(zcat clean.vcf.gz) | cut -f4,5 | sed -e "s/^./[&/g" -e "s/.$/&]/g" -e "s_\t_/_g" > tmp2

# Merge to final table of flanking sequences
paste <(cut -f1 tmp1 | sed "s/_LEFT//g") <(paste <(cut -f2 tmp1) tmp2 <(cut -f4 tmp1) | tr -d "\t") > flanking-seqs.tsv

# Next, remove any SNPs where one of the flanking sequences is less than 80 bp. Add a tab space at the end for searching later
seqtk comp flanks.tmp.fa | cut -f1-2 | awk '$2 <= 80' | cut -f1 | sed "s/_LEFT\|_RIGHT//g" | sort | uniq | sed 's_$_\t_g' > remove.IDs
   # A total of 205,221 SNP loci need to be removed

# Convert to a fasta file to use for BLASTING to zebra finch genome reference. Keep reference allele (first allele)
grep -v -f remove.IDs flanking-seqs.tsv | tab2fasta | sed "s_\[\(.\)/.\]_\1_g" > flanking-subset.fasta
grep -v -f remove.IDs flanking-seqs.tsv > flanking-subset.tsv
   # A total of 812,823 SNP loci remain
```

_fasta2tab_
```perl
#!/usr/bin/env perl
# fasta2tab
# Johan Nylander
local $/ = '>';
while(<>) {
    chomp;
    next if($_ eq '');
    my ($h, @S) = split /\n/;
    my $s = join('', @S);
    print STDOUT "$h\t$s\n" unless (!$h);
}
```

_tab2fasta_
```perl
#!/usr/bin/env perl
# tab2fasta
# Johan Nylander
while(<>) {
    chomp;
    my ($h, $s) = split /\t+/;
#    $s =~ s/\S{60}/$&\n/sg;   # don't wrap to 60 bases per line
    print STDOUT ">$h\n$s\n";
}
```

#### Step 2.3: Remove any Super Transcripts (i.e. Genes) with a TPM count <1
In this section, the trimmed/cleaned sequencing reads used in the Trinity assembly are mapped back to the Super Transcipts using [bowtie2 v2.4.2](https://github.com/BenLangmead/bowtie2). After mapping, [rsem v1.3.3](https://deweylab.github.io/RSEM/) is used to estimate the expression level for each Super Transcipt in **TPM** (_transcripts per million_). The metric **TPM** scales expression for each transcript by the gene length and sequencing depth. In general, sequences with a **TPM**<1 are often considered to lack support for downstream analyses.  In [Kaiser et al. 2017](https://doi.org/10.1111/1755-0998.12589), the authors omitted contigs (i.e., trasncripts) with a **TPM**<2. The `align_and_estimate_abundance.pl` perl script used is available in the `util` folder as part of the Trinity package.

_**2.3.1:** align reads and estimate abundance (TPM)_
```bash
# Marge together the forward and reverse read files
cat HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.1.gz HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.1.gz > left.fq.gz
cat HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.2.gz HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.2.gz > right.fq.gz

# Align and estimate abundance (bowtie2 and rsem must be in the PATH)
 align_and_estimate_abundance.pl \
	 --transcripts Trinity.SuperTrans.fasta \
	 --seqType fq \
	 --left left.fq.gz \
	 --right right.fq.gz \
	 --est_method RSEM \
	 --aln_method bowtie2 \
	 --prep_reference \
	 --thread_count 32 \
	 --SS_lib_type RF \
	 --output_dir rsem_outdir 

# Remove the read files at the end.
rm -rf left.fq.gz right.fq.gz
```

_**2.3.2:** Generate list of transcripts with TPM≥1_
```bash
# Move counts/expression file to 'SNPs' folder
mv rsem_outdir/RSEM.genes.results ../../SNPs

# How to sort a file with headers (e.g., sort by TPM)
cat RSEM.genes.results | (sed -u 1q; sort -nr -k6,6) | less

# Generate list of transcripts with TPM≥1
   # Add "_" to end of each trasncript ID to improve search
awk '$6 >= 1' RSEM.genes.results | cut -f1 | sed '1d' | sed 's/$/_/g' > Super.Transcripts_TPM1.list
   # 17,682 transcripts in total selected
```

_**2.3.3:** Generate a final subset of SNP loci for transcripts with TPM≥1_
```bash
# Subset FASTA and TSV files of flanking sequences for just these Super Transcripts
   # sed command is to remove the extra '--' inserted by grep
grep -f Super.Transcripts_TPM1.list -A1 flanking-subset.fasta | sed '/^--/d' > flanking-subset-TPM1.fasta
grep -f Super.Transcripts_TPM1.list flanking-subset.tsv | sed '/^--/d' > flanking-subset-TPM1.tsv

# How many sequences are left?
grep -c "^>" flanking-subset-TPM1.fasta
   # 262,871
```
To summarize, at the end of this step (**Step 2.3**), a total of 262,871 putative SNPs remain.  These SNPs have been filtered for all sorts of identification criteria, and are also restricted to transcripts with reasonably high support from the assembly for being true transcripts.

#### Step 2.4: BLAST to the zebrafinch reference genome
In this section, we BLAST each flanking region and remove any loci that match multiple locations.  These are indicative of potential paralogs and not true SNPs. As a reference, we are using the [zebra finch reference genome](https://www.ncbi.nlm.nih.gov/genome/?term=txid59729[orgn]).  In particular, the most recent assembly [GCA_003957565.4 bTaeGut1.4.pri](https://www.ncbi.nlm.nih.gov/genome/367?genome_assembly_id=1613864).  This assembly was recently built by the Vertebrate Genomes Project and contains a chromosome-level assembly composed of 199 scaffolds.  Although the Carolina wren (_Thryothorus ludovicianus_) has a [genome available](https://www.ncbi.nlm.nih.gov/genome/?term=txid74200[Organism:noexp]) and is more closely related to the house wren, the zebra finch genome is a much higher quality assembly.

_**2.4.1:** Gather zebra finch genome and setup BLAST database_
```bash
# Download and build BLAST DB for the Zebrafinch bTaeGut1.4
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz
mv GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz bTaeGut1.4.fna.gz
gunzip bTaeGut1.4.fna.gz

# Build BLAST database
makeblastdb -in bTaeGut1.4.fna -dbtype nucl -title bTaeGut1.4 -out bTaeGut1.4
gzip bTaeGut1.4.fna
```

_**2.4.2:** Perform the BLAST search_
```bash
# Load BLAST module (v2.11.0+)
module load blast+

# Set env variable to blast database locations
export BLASTDB=$PWD:$BLASTDB

# Blast Abyss contigs using 32 cores (only took a few minutes to run)
blastn \
   -query flanking-subset-TPM1.fasta \
   -db  bTaeGut1.4 \
   -max_target_seqs 5 \
   -max_hsps 1 \
   -evalue 1e-3 \
   -outfmt '6' \
   -out flanking-subset-TPM1.blastout \
   -num_threads 32
```

_**2.4.3:** Retain loci with only a single BLAST hit_
```bash
# Get list of loci with a single hit (for searching the tsv file, add tab character at the end)
cut -f1 flanking-subset-TPM1.blastout | \
   sort | \
   uniq -c | \
   sed "s/^ *//g" | \
   awk '$1 == 1' | \
   cut -d' ' -f2 | \
   sed 's_$_\t_g' > 1-hit-seqs.list

# Filter SNPs and FASTA to retain these IDs.
seqtk subseq -l0 ../flanking-subset-TPM1.fasta 1-hit-seqs.list > ../flanking-subset-TPM1-1hit.fasta
grep -f 1-hit-seqs.list flanking-subset-TPM1.tsv | sed '/^--/d' > flanking-subset-TPM1-1hit.tsv
   # 208,986 loci with a single blast hit
```

#### Step 2.5: Remove SNPs overlapping a splice junction or masked region
In this section, we will remove every SNP that overlaps either 1) a region in the reference genome that has been soft masked or 2) a region that contains an annotate splice junction. For the formers, the masked bases in the reference genome (made lower case letters) are indicative of either low quality bases or bases that contain an annotated repetive element. For the latter, a splice junction is an exon-intron junction, as we do want a SNP near a junction since the fragment size amplified by PCR may differ from that predicted/expected for the Sequenom analysis. We will combine these two "intervals" (GFF/BED files of genomic locations) into a single file then retains SNPs not overlapping these regions based on the blast results.

_**2.5.1:** Make a BED file of masked regions in the reference genome (long perl command)_
```bash
# This only works if each FASTA sequence is one-per-line (not wrapped)
seqtk seq -l0 bTaeGut1.4.fna | \
   perl -lne 'if(/^>(\S+)/){ $n=$1} else { while(/([a-z]+)/g){ printf("%s\t%d\t%d\n",$n,pos($_)-length($1),pos($_)) } }' > bTaeGut1.4.masked.bed
````
_**2.5.2:** Make a BED file of splice junctions in the reference genome_
```bash
# Get annotations
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.gff.gz
mv GCF_003957565.2_bTaeGut1.4.pri_genomic.gff.gz bTaeGut1.4.gff.gz
gunzip bTaeGut1.4.gff.gz

# Convert annotations to GTF format
gffread -E bTaeGut1.4.gff -T -o bTaeGut1.4.gtf

# Get splice junctions/introns
python3 hisat2_extract_splice_sites.py bTaeGut1.4.gtf > bTaeGut1.4.sj.bed
````

_**2.5.3:** Make a BED file of BLAST Hit locations in the reference genome_
```bash
# Get BLAST 1hit subset
grep -f 1-hit-seqs.list flanking-subset-TPM1.blastout > flanking-subset-TPM1-1hit.blastout

# Convert to bed format
blast2bed flanking-subset-TPM1-1hit.blastout
   # output is flanking-subset-TPM1-1hit.blastout.bed
```
_**2.5.4:** Retain SNPs without any overlaps_
```bash
# Get non-overlapping SNPs
bedtools intersect -v -a flanking-subset-TPM1-1hit.blastout.bed -b bTaeGut1.4.masked.bed | \
   bedtools intersect -v -a - -b bTaeGut1.4.sj.gtf | \
   cut -f4 | \
   sed 's_$_\t_g' > TPM1-1Hit-noOverlaps.list
  # starting with 208,986 loci, 31,268 did not overlap masked regions, and 8,257 of these did not overlap a splice junction.
  
# Filter SNPs and FASTA to retain these IDs.
seqtk subseq -l0 flanking-subset-TPM1-1hit.fasta TPM1-1Hit-noOverlaps.list > flanking-subset-TPM1-1hit-noOverlaps.fasta
grep -f TPM1-1Hit-noOverlaps.list flanking-subset-TPM1-1hit.tsv | sed '/^--/d' > flanking-subset-TPM1-1hit-noOverlaps.tsv
   # 8,257 loci with a single blast hit
```
_blast2bed_
```bash
#!/bin/bash

INPUT="$1"

#check if input file is provided
if [ -f "$INPUT" ];
then
    OUTPUT=$(echo $INPUT | sed -e s/.bls$//)'.bed'
else
    echo "ERROR: No input file prvided or file does not exist.
Usage: ./blast2bed <blastoutput.bls>
The blast file should be in blast outfmt 6 or 7.
See Readme.org for more details."
    exit 1
fi

#converting blast to bed
echo "converting $INPUT in $OUTPUT"

# echo "#This file was generated from $INPUT using blast2bed" > $OUTPUT
#grep -v '^#' "$INPUT" | awk '{print $2,"\t",$9-1,"\t",$10,"\t",$1}' | sort >> $OUTPUT
grep -v '^#' "$INPUT" | perl -ane 'if($F[8]<=$F[9]){print join("\t",$F[1],$F[8]-1,$F[9],$F[0],"0","+"),"\n";}else{print join("\t",$F[1],$F[9]-1,$F[8],$F[0],"0","-"),"\n";}' | sort >> $OUTPUT


echo "done"
```

#### Summary of the SNP filtering procecure
1. Started with **2,200,840** called variants
2. Of these, **1,299,784** variants passed the initial quality filter
3. Of these, **1,018,044** were biallelic SNPs that passed the coverage filter and were called heterozygous.
4. Of these, **812,823** SNPs had adequate flanking sequence (≥80 bp) for primer design
5. Of these, **262,871** SNPs were found on Super Transcripts with a TPM > 1 (17,682 Super Transcripts)
6. Of these, **208,986** SNPs (including flanking region) had a single BLAST hit to the zebra finch reference genome
7. Of these **31,268** SNPs did not overlap a masked region in the reference genome
8. Of these, **8,257** SNPs did not overlap a splice junction (intron-exon boundary) in the referene genome.
    - 2,677 Super Transcripts represented

_Build Final output table_
```bash
# Build final tab-separated table
paste \
   <(cut -f1 flanking-subset-TPM1-1hit-noOverlaps.tsv | sed "s/_[0-9]*$//g") \
   <(cut -f1 flanking-subset-TPM1-1hit-noOverlaps.tsv | cut -d"_" -f5) \
   <(cut -f1 flanking-subset-TPM1-1hit-noOverlaps.tsv) \
   <(cut -f1 flanking-subset-TPM1-1hit-noOverlaps.tsv | sed 's_$_\t_g' | grep -f - BLAST/flanking-subset-TPM1-1hit.blastout | cut -f 2,9,10) \
   <(cut -f2 flanking-subset-TPM1-1hit-noOverlaps.tsv) | \
   cat <(echo -e "SuperTranscript\tPosition\tSNP_ID\tbTaeGut1.4_Chrom\tbTaeGut1.4_start\tbTaeGut1.4_end\tSequence") - > final.table.tsv
```
