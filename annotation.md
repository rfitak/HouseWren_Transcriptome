# Annotating the House Wren Transcriptome
This section perform several quality control assessments annd annotation of the _de novo_ transcriptome assembly completed in the [previous section](./assembly.md).  Although a fulll annotation is not necessary to identify SNPs for parentage analysis, a fully annotated transcriptome cann be submitted to the NCBI's [TSA](https://www.ncbi.nlm.nih.gov/genbank/tsa/) database and will be a useful resource for other researchers.

## Step 1: Assess assembly using Transrate
This step will perform an assessment of the _de novo_ assembly using [Transrate v1.03](http://hibberdlab.com/transrate/). According to the Transrate website:  
"Transrate is software for de-novo transcriptome assembly quality analysis. It examines your assembly in detail and compares it to experimental evidence such as the sequencing reads, reporting quality scores for contigs and assemblies. This allows you to choose between assemblers and parameters, filter out the bad contigs from an assembly, and help decide when to stop trying to improve the assembly."

_Install Transrate_
```bash
# These are prebuilt linux x86_64 binaries so nothing more is really needed.
wget https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
tar -zxvf transrate-1.0.3-linux-x86_64.tar.gz

# Set correct permissions
chmod -R go-w /home/rfitak/ 

# Verify dependencies
./transrate --install-deps=ref
```


## Step 2: Identify Open Reading Frames (ORFs)

_Install [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)_
```bash
# Grab from GitHub
git clone https://github.com/TransDecoder/TransDecoder.git

# Install (verify hmmer v3.3 is installed and executables are in your PATH)
cd TransDecoder
make test
```

