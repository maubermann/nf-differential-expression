# Differential expression analysis workflow with newtflow 

This project provides a simple pipeline that aims to identify genes that are differently expressed in two distinct groups of RNAseq samples. 
The workflow is designed to take in raw paired end short read RNAseq data in the form of fastq files as input.
After a successful run the user will be provided with a quality control report for each of the input samples, a MA-plot that visualizes the log2-fold change of all observed genes and a hitlist of genes that appear to be differentially expressed after Benjamini Hochberg correction. Details on how to interpret these results can be found in the **Workflow** section of this README. Further investigating the differentially expressed genes can give valuable insight into the mechanisms behind the underlying differences between the two groups of interest. A next step to gain further interpretable insight could be GO enrichment analysis - i consider implementing it in a future version of this workflow.


## Set up

To clone this repository
1. Change the working directory to the location where you want the cloned repository
2. Clone the repository using
```console
git clone https://github.com/maubermann/nf-differential-expression
```

## Dependencies

To ensure reporoducibility, containerized versions of the required packages have been employed. In order to run the workflow you will only need the ***docker engine*** and ***nextflow*** installed on your system.  
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
  The first line of the file has to be `sample  condition`, every other line needs to contain an experiment name followed by a tab and the name of the group the experiment belongs to (for example control, treatment). Note: The names of the experiments need to match the prefixes of the FASTQ files.
  There must only be two distinct group names. Refer to the provided [experiment_info.txt](experiment_info.txt).

* ### CSV file containing all experiment names
  The file specifying the experiment names needs to contain one line only that contains the names of the experiments seperated by `, `. Refer to the provided [allExperiments.csv](allExperiments.csv).

## Running the workflow
In order to run the workflow you need to ensure that the docker daemon is running. When using docker desktop, you will simply need to open docker desktop.  
Running the workflow on the provided example data does not not require additional command line arguments. Simply change your working directory to the root directory of this repository and run
```console
nextflow run deWorkflow.nf
```
When using your own data you will need to provide the reqired input files via the command line. Run
```console
nextflow run deWorkflow.nf --input_dir "<FASTQ_DIRECTORY>" --referenceTranscriptome "<FASTA_FILE>" --geneAnnotations "<GTF_FILE>" --experimentInfo "<GROUP_TABLE>" --allExperiments "<CSV_FILE>"
```
Here you will need to replace the arguments in the quotation marks with the required input paths (do not remove the quotation marks).  
Note: The listing of the input files from the previous section is in the same order as the command line arguments. 

## Workflow 

![plot](https://github.com/maubermann/nf-differential-expression/blob/main/supplementary/nextflow_dag.png)

### Fastp quality control and preprocessing  
This step is performed for every pair of FASTQ files. Fastp aims to remove bad quality bases from both ends of the reads and to trim off adapter sequences - which are remainings from primer sequences and don't provide biological signal. Further Fastp removes bad quality reads. Fastp provides a detailed report html report that can be interactively examined through a web browser. Further, a less detailed json report is provided that is more suitable for any downstream analysis of the quality control results. This step not only requires the FASTQ files but also a list containing all experiment names. This list is required to identify for every given experiment name in the list the related pair of FASTQ files. For further details on the employed methods refer to the fastp paper by [Chen et al](https://doi.org/10.1093/bioinformatics/bty560)

### Salmon indexing and quantification
In order to align (or rather quasi-map) the raw reads from the FASTQ files, Salmon requires building a trascriptome index from a reference transcriptome first. The reference trascriptome is composed of cDNA sequences that are complimentary to all known or predicted mRNA transcripts in the reference organism. Building an index from the transcriptome enables fast lookup for k-mers in the transcriptome.  
After building the index, for every sequencing experiment salmon quantifies the expressed genes using from paired reads. Salmon exploits the fast lookup of kmers to perform so called quasi-mapping of the reads which is a heuristic for performing traditional alignment against the reference. Salmon utilizes statistical models to correcto for different biases in the data, more speciffically varying read abundance from the ends of transcripts and from regions with high or low CG content. It should be noted, that some randoms components of the algorithm do not allow for perfect reproducibility of results. For further details refer to the Salmon paper by [Patro et al](https://www.nature.com/articles/nmeth.4197).

### DESeq2 differential expression analysis  
As DESeq2 requires gene-level count matrices the transcript abundance estimates are first mapped to their corresponding genes. The transcript to gene mapping is obtained from the annotations in the GTF input file. The gene-level count matrix contains conts for every gene for every sample. The input table, specifying the two distinct groups of samples is required to perform differential expression analysis. DESeq2 utilizes the median of ratios method to account for differences in library size (sequencing depth) by estimate size factors to obtain normalized read counts. which is very robust to outliers and differentially expressed genes. Under assumed negative biomial distribution of the read counds, while accounting for dispersion between samples, DESeq2 fits a log linear regression to estimate the log2 fold change between the two samples for every gene.

$$ \log _2 FC(gene_i) = \log \left(\frac{\text{expression of gene$_i$ in group 1}}{\text{expression of gene$_i$ in group 2}} \right)$$

Whether the log2 fold change is significantly different from zero 0 (in other words whether the appears to be differentially expressed in the two groups) is then determined using the Wald test, providing the user with a p-value for every gene. To adjust for multiple testing and reduce the false discovery rate, Benjamini-Hochberg procedure is applied. For further details refer to the DESeq2 paper by [Love et al](https://link.springer.com/article/10.1186/s13059-014-0550-8). 
After completion of differential expression analyisis the user is provided whith a MA-plot that visualizes for every gene its mean expression value and the estimated log2 fold change. Genes for which the adjusted p-value is below a threshold of 0.1 are highlighted in the plot. Further, genes with an adjusted p-value below a significance threshold of 0.05 are provided to the user in a hit list. 

## Output files 
Below i will give an overview over the output of the workflow and briefly discuss a reference output for the provided example data. Please keep in mind that the FASTQ files used for this analysis were subsampeled and therefoere strongly reduced in size. Due to this limitation, the results for this dataset might not be biologically meaningful.  
After a successful run of the workflow, the root directory of this repository will contain an additional directory named `results`. In `results/salmon` you will find the a directory for every sequencing experiment containing the quantification files that were created by Salmon (along with some other data). I will not discuss these files here, as they are not ment for interpretation by the user. In `results/fastp` you will find html and json reports for every sequencing experiments and the preprocessed FASTQ files. In `results/DESEQ2` you will find a MA-plot along with the hitlist genes that appear to be differentially expressed (adjusted p < 0.05).

### Fastp report 

The fastp reports are similar to FastQC reports, however they also provide information on the performed preprocessing. Here some details of the fastp html report for the sequencing experiment `SRR25436327` are discussed. The report indicates that 197,204 reads passed quality control, the remaining 2,796 reads were filtered out due to low quality, a high amount of N (bases that could not get called) or due to short length. 
![read_quality](https://github.com/maubermann/nf-differential-expression/blob/main/supplementary/read_quality_SRR25436327_1.png)
Above the read quality for `SRR25436327_1` are shown before (left) and after (right) preprocessing. The plot indicates for any given position a quality score for the basecalling of the four nucleotides. A score above 20 is considered good. One can observe that the quality scores slightly improve through processing, likely by excluding low quality reads. Generally speaking, the results show good read quality, even before filtering.

![base_contents](https://github.com/maubermann/nf-differential-expression/blob/main/supplementary/Base_contents_SRR25436327_1.png)
Above the base contents per position of the reads in `SRR25436327_1` are shown before (left) and after (right) preprocessing. The fact that certain bases are over- or underrepresented at the 5' end of the reads might indicate the presence of certain adapter sequences or primer residuals. While fastp attempted to detect and remove adapter sequences, the base contents at the 5' end don't change. Very slight improvements can be observed at the 3' ends of the reads, likely due to the filtering of low quality reads.

### DESeq2 
The outputs of DESeq2 provide indicate which genes appear to be differentially expressed. 
![MA](https://github.com/maubermann/nf-differential-expression/blob/main/supplementary/MA-plot.png)


