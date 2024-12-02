# Differential expression analysis workflow using newtflow 

This project provides a simple pipeline that aims to identify genes that are differently expressed in two distinct groups of RNAseq samples. 
The workflow is designed to take in raw paired end read data in the form of fastq files as input. The packages that were used are designed for short read data, but in theory this workflow will also run on long read data.
After a successful run the user will be provided with a quality control report for each of the input samples, a MA-plot that visualizes the log2-fold change of the observed genes and a hitlist of genes that are appear to be differentiallu expressed after Benjamini Hochberg correction. Details on how to interpret these results can be found below.

## Set up

To clone this repository
1. Change the working directory to the location where you want the cloned repository
2. Clone the repository using
```console
git clone https://github.com/maubermann/nf-differential-expression
```

## Dependencies

To ensure reporoducibilities containerized versions of all required packages have been employed. In order to run the workflow you will need the ***docker engine*** and ***nextflow*** installed on your system.  
I encourage the installation of [Docker Desktop](https://docs.docker.com/get-started/get-docker/) which provides a graphical user interface.  
For a more minimal installation of the docker engine please [go here](https://docs.docker.com/engine/install/) and follow the instructions specific to your operating system.
For a guide on how to install nextflow, plese [go here](https://www.nextflow.io/docs/latest/install.html).

In order to test this workflow I used:  
* Docker Desktop 4.34.3  
* Nextflow 24.10.0

## Input files:
Examples for all files and directories that are required to run this workflow are provided in this repository. Below I specify what files need to be provided and how the examples in this srepository were obtained. 

* ### A directory containing all paired end FASTQ files to be analyzed
  The fastq files all have to lie in the same parent directory. For any given RNAseq sequencing experiment the two files containing the paired reads must follow the naming convention `<EXPERIMENT_NAME>_1.fastq` and `<EXPERIMENT_NAME>_2.fastq`.  
  I provided example [FASTQ files](FASTQ_sub) that contain RNAseq data from mice that were chronically exposed to morphine in contrast to a control group that were trated with saline. The experiments were conducted by [Rao et al (2023)](https://pubmed.ncbi.nlm.nih.gov/37543246/).
  The data was found through [gene expression omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239387), a database repository of high throughput gene expression data and hybridization arrays, chips, microarrays. The data was then downloaded through the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) using the [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).
  ```console
  fasterq-dump SRR25436327 SRR25436328 SRR25436329 SRR25436330 SRR25436331 SRR25436332 SRR25436333 SRR25436334 --verbose
  ```
  In order to provide a more compact files to test this workflowlow, for every sequencing experiment I sampled 100,000 random paired reads using [seqtk](https://github.com/lh3/seqtk).
  ```console
  seqtk sample -s100 <EXPERIMENT_NAME>_1.fastq 100000 > sub_<EXPERIMENT_NAME>_1.fastq
  seqtk sample -s100 <EXPERIMENT_NAME>_2.fastq 100000 > sub_<EXPERIMENT_NAME>_2.fastq
  ```

  
  
* ### Reference transcriptome as fasta file
  The workflow can automatically handle compressed `.gz` fasta files.  
  I provided a reference transcriptome of a mouse strain that is closely related to the strain used by Rao et al. It was downloaded from [here](https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus_c57bl6nj/cdna/) using the Ensembl genome browser on 24.11.2024.
* ### GTF file containing annotations for the reference transcriptome
  The workflow can automatically handle compressed `.gz` GTF files.
  I provided a GTF file for the provided reference transcriptome. It was downloaded from [here](https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus_c57bl6nj/) using the Ensembl genome browser on 24.11.2024
* ### Text file with a table specifying which experiments belong to which group
  The first line of the file has to be `sample  condition`, every other line needs to contain an experiment name followed by a tab and the name of the group the experiment belongs to (for example control, treatment). 
  There must only be two distinct group names. Refer to the provided [experiment_info.txt](experiment_info.txt).

* ### CSV file containing all experiment names
  The file specifying the experiment names needs to contain one line that contains the names of the experiments seperated by `, `.  
