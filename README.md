<h3><p align="center">Supplementary Methods for:</p></h3>
<h2><p align="center">The <i>de novo</i> assembly of the house wren (<i>Troglodytes aedon</i>) trascriptome</p></h2>

<I><h5>Rachael DiSciullo<sup>1</sup>, Robert R. Fitak<sup>2</sup>, Anna Forsman<sup>2</sup>, Scott Sakaluk<sup>2</sup>, Charles Thompson<sup>1</sup></h5></I>

1. Department of Biology, Illinois State University, Normal, IL 61790 USA.
2. Department of Biology, Genomics and Bioinformatics Cluster, University of Central Florida, Orlando, FL 32816 USA.

***
___This GitHub repository contains a summary of the various code, software, and data analysis pipelines used for the aforementioned study of the house wrenn. The contents represented here are only to be used as an example and not intended to be comeprehensive. The authors make no representation about the suitability or accuracy of this code, software, or data for any purpose, and make no warranties, either expressed or implied, for a particular purpose or that the use of this software or data will not infringe any third party patents, copyrights, trademarks, or other rights. The code, software and data are provided "as is". All content is hereby registered under the GNU General Public License v3.0, see [LICENSE](./LICENSE). Any publication that significantly relies upon the use of the content generated herein shall appropriately cite:___

<p align="center">DiSciullo et al. (in prep.) The <i>de novo</i> assembly of the house wren (<i>Troglodytes aedon</i>) trascriptome. TBD</p>

***
  
<h2><p align="center">Table of Contents</p></h2>
<div align="center">
 
[Read processing, cleaning, and error-correction](./read_processing.md)

[SNP Calling - including base and variant quality score recalibrations](./SNP-calling.md)

</div>

***

<h2><p align="center">Summary of the Programs/Software Used in this Study</p></h2>  

| Program | Version | Citation |
| --- | --- | --- |
| FASTP | 0.20.0 | [Chen et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884–i890 (2018).](https://doi.org/10.1093/bioinformatics/bty560) |
| RCORRECTOR | cccc | [Song L, Florea L. Rcorrector: Efficient and accurate error correction for Illumina RNA-seq reads. GigaScience 4, 48 (2015).](https://doi.org/10.1186/s13742-015-0089-y) |
| BWA | 0.6.2 | Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25, 1754-1760 (2009). |
| SAMTOOLS | 0.1.19 | Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079 (2009). |
| VCFTOOLS | 0.1.12b | Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, et al. 2011. The variant call format and VCFtools. Bioinformatics. 27:2156–8. |

