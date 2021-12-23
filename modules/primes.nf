
nextflow.enable.dsl=2

// TODO: remove since its depreciated by the functions in seqkit
// Steps to extract prime end of fastq files
process Extract5PrimeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    container "staphb/seqtk:1.3"
    
    input:
        path fasta
        val length

    output:
        path "*_5p.fasta"

    script:
        println("Extract5PrimeFasta has been depreciated!")
        """
        seqtk trimfq -L $length $fasta > ${fasta.simpleName}_${length}_5p.fasta
        """
}

process Extract3PrimeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "staphb/seqtk:1.3"

    input:
        path fasta
        val length

    output:
        path "*_3p.fasta"

    script:
        println("Extract3PrimeFasta has beendepreciated!")
    // flip, trim, flip back
        """
        seqtk seq -r $fasta | seqtk trimfq -L $length - | seqtk seq -r - > ${fasta.simpleName}_${length}_3p.fasta
        """
}

