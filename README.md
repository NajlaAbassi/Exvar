# exvar: An R Package for Gene Expression and Genetic Variation Data Analysis and Visualization

<img src="inst/logo.png" align="right" alt="" width="180" />

[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic_2.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13459761.svg)](https://doi.org/10.5281/zenodo.13459761)

<<<<<<< HEAD
=======
<p align="center">
  <img src="https://raw.githubusercontent.com/omicscodeathon/exvar/main/logo.png" alt="exvar logo" width="300"  />
</p>

>>>>>>> upstream/main

## Table of Contents
1. [Motivation](#Motivation)
2. [About exvar](#About-exvar)
3. [Supported Species](#Supported-Species)
4. [Installation](#Installation)
5. [Functions](#Functions)
6. [License](#License)
7. [Contact](#Contact)
8. [Citation](#Citation)
9. [Contributors](#Contributors)


## Motivation

RNA sequencing data manipulation workflows are complex and require various skills and tools. This creates the need for user-friendly and integrated genomic data analysis and visualization tools.

We developed a novel R package using multiple R packages to perform gene expression analysis and genetic variant calling from RNA sequencing data. The package could be used to analyze eight species’ data. Multiple public datasets were analyzed using the developed package to validate the pipeline for all the supported species.

## About exvar

<<<<<<< HEAD
The Exvar R package performs gene expression and  genetic variants (SNPs, Indels, and CNVs) data analysis and  visualization.
=======
The exvar R package performs gene expression and  genetic variats (SNPs, Indels, and CNVs) data analysis and  visualization.
>>>>>>> upstream/main

It is developed to facilitate and minimize the skills required for the genetic expression and variants calling. It includes 9 functions, providing an easy biologist-friendly workflow.

## Supported Species

It could be used to analyze data from eight species including:

    * Homo Sapiens (hg19, 37, and 38)

    * Mus Musculus

    * Arabidopsis Thaliana

    * Drosophila Melanogaster

    * Danio rerio

    * Rattus norvegicus

    * Caenorhabditis elegans  

    * and Saccharomyces Cerevisiae.

### Limitations:

| Function/Species           | processfastq | count | geneExpression | callsnp | callcnv | callindels | vizexp | vizsnp | vizcnv |
|-----------------------------|--------------|--------|------------|---------|---------|-----------|--------|--------|--------|
| **Homo Sapiens**           | ✓            | ✓      | ✓          | ✓       | ✓       | ✓         | ✓      | ✓      | ✓      |
| **Mus Musculus**           | ✓            | ✓      | ✓          | ✓       | ✗       | ✓         | ✓      | ✓      | ✓      |
| **Arabidopsis Thaliana**   | ✓            | ✓      | ✓          | ✓       | ✗       | ✓         | ✓      | ✓      | ✗      |
| **Drosophila Melanogaster**| ✓            | ✓      | ✓          | ✓       | ✗       | ✓         | ✓      | ✓      | ✓      |
| **Danio rerio**            | ✓            | ✓      | ✓          | ✓       | ✗       | ✓         | ✓      | ✓      | ✓      |
| **Rattus norvegicus**      | ✓            | ✓      | ✓          | ✓       | ✗       | ✓         | ✓      | ✓      | ✓      |
| **Caenorhabditis elegans** | ✓            | ✓      | ✓          | ✓       | ✗       | ✓         | ✓      | ✓      | ✗      |
| **Saccharomyces Cerevisiae**| ✓           | ✓      | ✓          | ✓       | ✗       | ✓         | ✓      | ✓      | ✓      |

## Installation

The package could be installed as follows:

    install.packages("devtools")

    library(devtools)

<<<<<<< HEAD
    devtools::install_github("omicscodeathon/Exvar")
=======
    devtools::install_github("omicscodeathon/exvar/Package")
>>>>>>> upstream/main

    library(exvar)


## Functions

The package consists of 6 data analysis functions (processfastq(); count(); geneExpression(); callsnp(); callcnv(); and callindels()), 3 data visualization functions (vizexp(), vizsnp(), and vizcnv()), and a requirement() function is also provided to install all the dependencies. The list of functions and their roles are represented in Table 1.

| Function        | Role                                   | Input                    | Output                     |
|-----------------|----------------------------------------|--------------------------|-----------------------------|
| `requirement()` | Install required packages              | -                        | -                           |
| `processfastq()`| Preprocess fastq files                 | Fastq files              | BAM files                   |
| `count()`      | Gene count analysis                    | BAM files                | CSV file                    |
| `geneExpression()`  | Identify DEGs                          | BAM files                | CSV file                    |
| `callsnp()`     | SNP calling                            | BAM files                | VCF files                   |
| `callcnv()`     | CNV calling                            | BAM files                | CSV file                    |
| `callindels()`   | Indel calling                          | BAM files                | VCF files                   |
| `vizexp()`      | Analyze and visualize gene expression  | CSV file                 | Interactive interface       |
| `vizsnp()`      | Analyze and visualize SNP data          | VCF files                | Interactive interface       |
| `vizcnv()`      | Analyze and visualize CNV data          | CSV file                 | Interactive interface       |


#### Operating system :
- Linux for the processfastq(), count(), geneExpression(), callsnp(), callcnv(), and callindels() functions

- The visualization functions vizexp(), vizsnp(), and vizcnv() are platform-independent.

#### Programming language: R 4.2.1

#### Requirements:
The RAM and CPU requirements differ according to the analyzed data size. For Human whole genome data 16  CPUs and  65 GB RAM are required.

## License  

Artistic license 2.0
Any restrictions to use by non-academics: None.

## Citation

Ben Aribi, H., Dixon, I., Abassi, N., & Awe, O. I. (2024). exvar: An R Package for Gene Expression and Genetic Variation Data Analysis and Visualization. Zenodo. https://doi.org/10.5281/zenodo.13459761

## Contributors

   - Hiba Ben Aribi, Faculty of Sciences of Tunis, University of Tunis El Manar, Tunis, Tunisia.

   - Imraan Dixon, Faculty of Health Sciences, University of Cape Town, Cape Town, South Africa.

   - Najla Abassi, Higher Institute of Biotechnology Sidi Thabet, Manouba University, Tunisia.

   - Olaitan I. Awe, Ph.D., African Society for Bioinformatics and Computational Biology, Cape Town, South Africa.
