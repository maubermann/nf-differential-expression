params.input_dir = "${projectDir}/example_data/FASTQ_sub"

params.allExperiments = "${projectDir}/example_data/allExperiments.csv"
params.experimentInfo = "${projectDir}/example_data/experiment_info.txt"

params.referenceTranscriptome ="${projectDir}/example_data/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.cdna.all.fa.gz"
params.geneAnnotations = "${projectDir}/example_data/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.113.gtf.gz"



//The fastQC process is unused 
process fastQC{
    publishDir 'results', mode: 'symlink'
    container 'staphb/fastqc:0.12.1'
    
    input:
        path fastq_file

    output:
        path "${fastq_file.baseName}_fastqc.html"
        path "${fastq_file.baseName}_fastqc.zip"

    script:
    """
    fastqc ${fastq_file}
    """ 

}



process FASTP {
    publishDir 'results/fastp', mode: 'symlink'
    container 'biocontainers/fastp:v0.19.6dfsg-1-deb_cv1'

    input:
        path fastq_dir
        val experiment_name

    output:
        tuple val("${experiment_name}"), path("t_${experiment_name}_1.fastq"), path( "t_${experiment_name}_2.fastq") ,emit: trimmedTuple
        path "${experiment_name}_fastp_report.html"    ,emit:  htmlRep
        path "${experiment_name}_fastp_report.json"    ,emit:  jsonRep

    script:
    """
    fastp -i ${fastq_dir}/${experiment_name}_1.fastq \
            -I ${fastq_dir}/${experiment_name}_2.fastq \
            -o t_${experiment_name}_1.fastq \
            -O t_${experiment_name}_2.fastq \
            -h ${experiment_name}_fastp_report.html\
            -j ${experiment_name}_fastp_report.json\
            --detect_adapter_for_pe --cut_front --cut_tail
    
    """
    //outputs trimmed fastq files and quality reports
}

process SALMON_INDEX{
    container 'combinelab/salmon:1.10.3'

    input:
    path input_transcriptome

    output:
    path "index_file"

    script:
    """
    salmon index -t ${input_transcriptome} -i index_file
    """
    //indexing step for quasi-mapping in SALMON_QUANTIFY
}

process SALMON_QUANTIFY{

    publishDir 'results/salmon', mode: 'symlink'
    container 'combinelab/salmon:1.10.3'
    input:
    tuple val(experiment_name), path(fastq_file_1), path(fastq_file_2)
    path index_file

    output:
    path "${experiment_name}_quant"

    script:
    """
    salmon quant -i ${index_file} -l A \
        -1 ${fastq_file_1} \
        -2 ${fastq_file_2} \
        --validateMappings \
        -o ${experiment_name}_quant
    """
    //recommended additional flag: --gcBias
    //Transcript quantification for a pair of FASTQ files

}


process DESEQ2{

    container 'resolwebio/rnaseq:5.11.0'
    publishDir 'results/DESEQ2', mode: 'symlink'

    input:
    path all_quant_dir
    path experiment_info
    
    path deseq2_r_script 
    path annotations

    output:
    path "genes_padj_significant.csv"
    path "MA_plot.pdf"

    script:
    """
    chmod +x ${deseq2_r_script}
    ./${deseq2_r_script} . ${experiment_info} ${annotations} 
    """
    //the first line of the script ensures, that the script is made executable
    // R script performs basic differential expression analysis with the DESeq2 R package

}

workflow  {

    //declaration of all external files
    fastq_directory = file(params.input_dir)
    
    transcriptome = file(params.referenceTranscriptome)
    annotation_file = file(params.geneAnnotations)

    experiment_info = file(params.experimentInfo)
    all_filenames_ch = Channel.fromPath(params.allExperiments).splitCsv().flatten()

    rscript_path = "${projectDir}/DESeq2_DEanalysis.R"

    //fastp process
    FASTP(fastq_directory, all_filenames_ch)

    //salmon processes
    SALMON_INDEX(transcriptome)
    SALMON_QUANTIFY(FASTP.out.trimmedTuple, SALMON_INDEX.out)

    //DESeq2 analysis
    DESEQ2(SALMON_QUANTIFY.out.collect(), experiment_info, rscript_path, annotation_file)

}
