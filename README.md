<h3><p align="center">Supplementary Methods for:</p></h3>
<h2><p align="center">The <i>de novo</i> assembly of the house wren (<i>Troglodytes aedon</i>) transcriptome</p></h2>

<I><h5>Rachael DiSciullo<sup>1</sup>, Anna Forsman<sup>2,3</sup>, Robert R. Fitak<sup>2</sup>, John Hunt<sup>4</sup>, Pirmin Nietlisbach<sup>1</sup>, Charles F. Thompson<sup>1</sup>, Scott K. Sakaluk<sup>1</sup></h5></I>

1. Behavior, Ecology, Evolution, and Conservation Section, School of Biological Sciences, Illinois State University, Normal, IL, USA
2. Department of Biology, Genomics and Bioinformatics Cluster, University of Central Florida, Orlando, FL 32816 USA.
3. Current location: Department of Biology, Colby College, Waterville, ME, USA
4. School of Science, Western Sydney University, Penrith, Australia

<br>

<p align="center">
  <img src="images/DiSciullo_Image of Research at Illinois State 2020_image.jpg" alt="House Wren" width="500">
</p>
<p align="center"><sup>Photo of a male house wren (<i>Troglodytes aedon</i>) singing.  Photo copyrighted by Rachael DiSciullo, please do not use or distribute without permission.</sup>
</p>

<br>

***
___This GitHub repository contains a summary of the various code, software, and data analysis pipelines used for the aforementioned study of the house wren. The contents represented here are only to be used as an example and not intended to be comeprehensive. The authors make no representation about the suitability or accuracy of this code, software, or data for any purpose, and make no warranties, either expressed or implied, for a particular purpose or that the use of this software or data will not infringe any third party patents, copyrights, trademarks, or other rights. The code, software and data are provided "as is". All content is hereby registered under the GNU General Public License v3.0, see [LICENSE](./LICENSE). Any publication that significantly relies upon the use of the content generated herein shall appropriately cite:___

<p align="center">DiSciullo et al. (2024) Male song structure predicts offspring recruitment to the breeding population in a migratory bird. <i>Evolution</i> 78(6): 1054–1066. [https://doi.org/10.1093/evolut/qpae036](https://doi.org/10.1093/evolut/qpae036)</p>

***
  
<h2><p align="center">Table of Contents</p></h2>
<div align="center">
 
[Read processing, cleaning, and error-correction](./read_processing.md)

[de novo assembly using Trinity](./assembly.md)

[SNP Calling](./snp-calling.md)

</div>

***

<h2><p align="center">Summary of the Programs/Software Used in this Study</p></h2>  

| Program | Version | Citation |
| --- | --- | --- |
| FASTP | 0.20.0 | [Chen et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884–i890 (2018).](https://doi.org/10.1093/bioinformatics/bty560) |
| RCORRECTOR | 1.0.4 ce5d06b | [Song L, Florea L. Rcorrector: Efficient and accurate error correction for Illumina RNA-seq reads. GigaScience 4, 48 (2015).](https://doi.org/10.1186/s13742-015-0089-y) |
| FilterUncorrectabledPEfastq.py | N/A | https://github.com/harvardinformatics/TranscriptomeAssemblyTools |
| BOWTIE2 | 2.3.5.1 | [Langmed, B. & Salzberg, SL. Fast gapped-read alignment with Bowtie 2. Nature Methods 9(4), 357–359 (2012](https://dx.doi.org/10.1038%2Fnmeth.1923) |
| SILVA database | SSU:138, LSU:132 | [Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Research 41, D590-D596 (2013)](http://nar.oxfordjournals.org/content/41/D1/D590) |
| TRINITY | 2.9.1 | [Grabherr, M., Haas, B., Yassour, M. et al. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat Biotechnol 29, 644–652 (2011)](https://doi.org/10.1038/nbt.1883) |
| SAMTOOLS | 1.9 | [Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079 (2009).](https://doi.org/10.1093/bioinformatics/btp352) |
| PICARD | 2.23.0 | [Picard Toolkit.” 2018. Broad Institute, GitHub Repository. http://broadinstitute.github.io/picard/; Broad Institute](https://broadinstitute.github.io/picard/) |
| GATK | 4.1.7.0 | [Van der Auwera, G.A., et al. From FastQ data to high‐confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Current Protocols in Bioinformatics, 43: 11.10.1-11.10.33 (2013).](https://doi.org/10.1002/0471250953.bi1110s43) |
| STAR | 2.7.3a | [Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29(1), 15‐21 (2013).](https://doi.org/10.1093/bioinformatics/bts635) |
| VCFTOOLS | 0.1.12b | Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et al. 2011. The variant call format and VCFtools. Bioinformatics. 27:2156–8. |

