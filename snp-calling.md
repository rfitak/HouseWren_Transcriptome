# Identifying SNPs
TBD

## Step 1:  Running Trinity's built-in GATK pipeline on the Super Transcripts
TBD

Dependencies that need to be installed:  
  - [Picard v2.23.0](https://broadinstitute.github.io/picard/)
    - `wget https://github.com/broadinstitute/picard/releases/download/2.23.0/picard.jar`
  - [STAR v2.7.3a](https://github.com/alexdobin/STAR)
  - [samtools v1.9](https://samtools.github.io)
  - [GATK v3.8]()

_Setup Picard and GATK environmental variables_
```bash
export PICARD_HOME=/home/rfitak/.bin/
export GATK_HOME=/home/rfitak/PROGRAMS/GATK
```

_Run SNP-calling Pipeline_
```bash
# Run Trinity SNP Pipeline
run_variant_calling.py \
   --st_fa Trinity.SuperTrans.fasta \
   --st_gtf Trinity.SuperTrans.gtf \
   -p /home/rfitak/fitaklab/HOUSE_WREN/HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.1.gz,/home/rfitak/fitaklab/HOUSE_WREN/HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.1.gz /home/rfitak/fitaklab/HOUSE_WREN/HOWR-1_cleaned.cor.unfixrm.rmrRNA.fq.2.gz,/home/rfitak/fitaklab/HOUSE_WREN/HOWR-2_cleaned.cor.unfixrm.rmrRNA.fq.2.gz \
   -o SNPS \
   -l 150 \
   -t 32
```
