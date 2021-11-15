
nextflow.enable.dsl=2


process SamtoolsIndex{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'

    input:
        each path(bam)

    output:
        tuple path(bam), path("*.bai")

    script:
        """
        samtools index $bam
        """

}

process SamtoolsSort{
    // Given a sam or bam file, make a sorted bam
    // Does break when sam is not `propper` eg no @SQ tags
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'

    input:
        path(input_sam)

    output:
        path "${input_sam.SimpleName}_sorted.bam"

    script:
        """
        samtools sort $input_sam -O bam -o "${input_sam.SimpleName}_sorted.bam"
        """

}

