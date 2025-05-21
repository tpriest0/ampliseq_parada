process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.6--py39hf95cd2a_1' :
        'biocontainers/cutadapt:4.6--py39hf95cd2a_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed_step1_r1  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.temp1.fastq.gz"
    def trimmed_step1_r2  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-p ${prefix}_2.temp1.fastq.gz"
    def trimmed_step2_r1  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.temp2.fastq.gz"
    def trimmed_step2_r2  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-p ${prefix}_2.temp2.fastq.gz"
    def trimmed_final  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    """
    cutadapt $args $trimmed_step1_r1 $trimmed_step1_r2 $reads > ${prefix}.cutadapt.log
    cutadapt $args2 $trimmed_step2_r1 $trimmed_step1_r1 >> ${prefix}.cutadapt.log
    cutadapt $args3 $trimmed_step2_r2 $trimmed_step1_r2 >> ${prefix}.cutadapt.log
    cutadapt $trimmed_final $trimmed_step1_r1 $trimmed_step1_r2 >> {log.log} 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "${prefix}.trim.fastq.gz" : "${prefix}_1.trim.fastq.gz ${prefix}_2.trim.fastq.gz"
    """
    touch ${prefix}.cutadapt.log
    touch ${trimmed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
