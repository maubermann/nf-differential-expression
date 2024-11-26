params.input_dir = "FASTQ_sub"

params.allExperiments = "allExperiments.csv"
params.experimentInfo = "experiment_info.txt"

params.tx2gene = "tx2gene.csv"
params.referenceTranscriptome ="Mus_musculus_c57bl6nj.C57BL_6NJ_v1.cdna.all.fa.gz"
params.rscript = "DESeq2_DEanalysis.R"

process fastQC{

    publishDir 'results', mode: 'symlink'
    //container 'maubermann/test:1.0'
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
    //315 mds and 151


}


process fastp{

    publishDir 'results/fastp', mode: 'symlink'
    container 'maubermann/test:1.0'
    //container 'biocontainers/fastp:v0.20.1_cv1'
    //trims off adapters for paired read FASTQ files
    
    input:
        path fastq_file_1
        path fastq_file_2
        //tuple path( fastq_file_1), path(fastq_file_2)
        //other files that wouldn't be called directly but that are needed in the wd would go here

    output:
        tuple path("${fastq_file_1.baseName}_trim.fastq"), path("${fastq_file_2.baseName}_trim.fastq") ,emit: trim
        //same tuple could go in input
        path "${fastq_file_1.baseName}_fastp_report.html"    ,emit:  htmlRep
        path "${fastq_file_1.baseName}_fastp_report.json"    ,emit:  jsonRep

    script:
    """
    fastp -i ${fastq_file_1} \
            -I ${fastq_file_2} \
            -o ${fastq_file_1.baseName}_trim.fastq \
            -O ${fastq_file_2.baseName}_trim.fastq \
            -h ${fastq_file_1.baseName}_fastp_report.html\
            -j ${fastq_file_1.baseName}_fastp_report.json\
            --detect_adapter_for_pe --cut_front --cut_tail
    """


}


process FASTP {
    publishDir 'results/fastp', mode: 'symlink'
    container 'staphb/fastp:0.23.4'

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
    //recommended: --gcBias

 
}

process DESEQ2{

    container 'resolwebio/rnaseq:5.11.0'
    publishDir 'results/DESEQ2', mode: 'symlink'

    input:
    path all_quant_dir
    path experiment_info
    path tx2gene
    path deseq2_r_script 

    output:
    path "genes_padj_significant.csv"
    path "MA_plot.pdf"

    script:
    """
    chmod +x ${deseq2_r_script}
    ./${deseq2_r_script} . ${experiment_info} ${tx2gene} 
    """

}

workflow  {
    //accessory files that arent main data inputs dont need own channel:   xfile = file(params.xfile)
    fastq_directory = file(params.input_dir)
    all_filenames_ch = Channel.fromPath(params.allExperiments).splitCsv().flatten()
    transcriptome = file(params.referenceTranscriptome)


    FASTP(fastq_directory, all_filenames_ch)
    SALMON_INDEX(transcriptome)
    SALMON_QUANTIFY(FASTP.out.trimmedTuple, SALMON_INDEX.out)

    
    tx2gene = file(params.tx2gene)
    rscript = file(params.rscript)
    experiment_info = file(params.experimentInfo)
    DESEQ2(SALMON_QUANTIFY.out.collect(), experiment_info,tx2gene,rscript)

    
}