# Assembling the House Wren Transcriptome using Trinity
This section focuses on generating the _de novo_ transcriptome assembly using the processed sequencing reads completed in the [previous section](./read_processing.md).  Afterwards, a series of cleaning, quality control, and annotation steps are performed prior to identifying SNPs.

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
