# Annotating the House Wren Transcriptome
This section perform several quality control assessments annd annotation of the _de novo_ transcriptome assembly completed in the [previous section](./assembly.md).  Although a fulll annotation is not necessary to identify SNPs for parentage analysis, a fully annotated transcriptome cann be submitted to the NCBI's TSA]() database and will be a useful resource for other researchers.

## Step 1: Assess assembly using Transrate

## Step 2: Identify Open Reading Frames (ORFs)

_Install [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)_
```bash
# Grab from GitHub
git clone https://github.com/TransDecoder/TransDecoder.git

# Install (verify hmmer v3.3 is installed and executables are in your PATH)
cd TransDecoder
make test
```

