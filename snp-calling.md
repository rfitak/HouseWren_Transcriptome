# Identifying SNPs
TBD

## Step 1:  Running Trinity's built-in GATK pipeline on the Super Transcripts
TBD

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

See various output and log files from the code above:
- [snp-calling.err](./Data/snp-calling.err)
- [variant_calling.log](./Data/variant_calling.log)
- [Log.out](./Data/Log.out)
- [Log.final.out](./Data/Log.final.out)


