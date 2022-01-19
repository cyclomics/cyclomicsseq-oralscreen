#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process TrimFasta {
    // sed call to remove all info after first comma in every line.
    // used to reduce the length of the readname in fasta files
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    
    input:
        each path(fasta)

    output:
        path "*_trimmed.fasta"

    script:
        """
        sed 's/,.*//' $fasta > ${fasta.SimpleName}_trimmed.fasta
        """
}


process ConcatenateFasta {
    // Call `cat` on all files (eg fasta) that enter the process
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    
    input:
        path(fasta)

    output:
        path "concat.fasta"

    script:
        """
        cat *.fasta > concat.fasta
        """
}

process splitSequences {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path(fasta)

    output:
        path("${fasta.SimpleName}_seq_*.fa")

    script:
        """
        awk '/^>/{f="${fasta.SimpleName}_seq_"++d".fa"} {print > f}' < $fasta
        """
}

process ConcatenateFastq {
    // Call `cat` on all files (eg fasta) that enter the process
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    
    input:
        path(fasta)

    output:
        path "concat.fastq"

    script:
        """
        cat *.fastq > concat.fastq
        """
}

process FilterBams{
    // given a tuple, filter on the fact that the second argument is bigger or equal to 1
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), val(filter_value), path(bam), path(bai)

    output:
        tuple val(X), path(bam), path(bai)

    when:
        Float num = "$filter_value" as Float
        num >= 1
    
    script:
    """
    """
}


process ConcatBams{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple val(X), path(bam), path(bai)

    
    script:
    """
    """
}

