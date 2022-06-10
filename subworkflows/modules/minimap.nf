#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process MinimapAlign{
    // Sams are large and should never be the final data point
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'minimap_large'

    input:
        tuple val(X), path(fasta), path(reference_genome)

    output:
        tuple val(X), path("${X}.sam") 

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fasta > ${X}.sam
        """
}

process MinimapAlignMany{
    // Same 
    // Sams are large and should never be the final data point
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'minimap_large'

    input:
        each path(fasta)
        path(reference_genome)

    output:
        tuple val("${fasta.simpleName}"), path("${fasta.simpleName}.sam") 

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fasta > ${fasta.simpleName}.sam
        """
}

process Minimap2AlignAdaptive{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    //  apply at least 1 Gb of memory to the process, otherwise we apply two-four times the size of the reference genome mmi
    memory {reference_genome.size() < 500_000_000 ? "1.5GB" : "${reference_genome.size() * (1 + task.attempt * task.attempt)} B"}
    // memory "32 GB"
    // small jobs get 4 cores, big ones 8
    cpus (params.ci_run == true ? 2 :{reference_genome.size() < 500_000_000 ? 4 : 8 })
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3


    input:
        each path(fasta)
        path(reference_genome)
    
    output:
        tuple val("${fasta.simpleName}"), path("${fasta.simpleName}.bam") 

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fasta > tmp.sam 
        samtools sort -o ${fasta.simpleName}.bam tmp.sam
        rm tmp.sam
        """
}

process Minimap2AlignAdaptiveParameterized{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    //  apply at least 1 Gb of memory to the process, otherwise we apply two-four times the size of the reference genome mmi
    memory {reference_genome.size() < 500_000_000 ? "1.5GB" : "${reference_genome.size() * (1 + task.attempt * task.attempt)} B"}
    // memory "32 GB"
    // small jobs get 4 cores, big ones 8
    cpus (params.ci_run == true ? 2 :{reference_genome.size() < 500_000_000 ? 4 : 8 })

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        each path(fasta)
        path(reference_genome)
    
    output:
        tuple val("${fasta.simpleName}"), path("${fasta.simpleName}.bam") 

    script:
    // Lower parameters to increase data available to cycas
        """
        minimap2 -ax map-ont -t ${task.cpus} -m ${params.minimap2parameterized.min_chain_score} -n ${params.minimap2parameterized.min_chain_count} -s ${params.minimap2parameterized.min_peak_aln_score} $reference_genome $fasta > tmp.sam 
        samtools sort -o ${fasta.simpleName}.bam tmp.sam
        rm tmp.sam
        """
}

process Minimap2Index{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'few_memory_intensive'

     input:
        path(reference_genome)
    
    output:
        path("${reference_genome.simpleName}.mmi")

    script:
        """
        minimap2  -ax map-ont -t ${task.cpus} -d ${reference_genome.simpleName}.mmi $reference_genome
        """
}
