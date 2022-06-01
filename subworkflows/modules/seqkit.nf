#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FastqToFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path(fastq)

    output:
        path("${fastq.SimpleName}.fa")

    script:
        """
        seqkit fq2fa $fastq -o ${fastq.SimpleName}.fa
        """
}

process Extract5PrimeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path fasta
        val length

    output:
        path "*_5p.fasta"

    script:
        """
        seqtk trimfq -L $length $fasta > ${fasta.simpleName}_${length}_5p.fasta
        """
}

process MergeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path fasta1
        path fasta2
    
    output:
        path "${fasta1.simpleName}_${fasta2.simpleName}.fasta"
    
    script:
        """
        cat $fasta1 > ${fasta1.simpleName}_${fasta2.simpleName}.fasta
        cat $fasta2 >> ${fasta1.simpleName}_${fasta2.simpleName}.fasta
        """

}

process Extract3PrimeFasta {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path fasta
        val length

    output:
        path "*_3p.fasta"

    script:
    // flip, trim, flip back
        """
        seqtk seq -r $fasta | seqtk trimfq -L $length - | seqtk seq -r - > ${fasta.simpleName}_${length}_3p.fasta
        """
}

process ExtractSpecificRead{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path fasta
        val readname

    output:
        path "*_${readname}.fasta"

    script:
        """
        seqkit grep -p $readname  $fasta  > ${fasta.simpleName}_${readname}.fasta
        """
}

process CountFastqInfo{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        path(fastq)

    output:
        path ("read_count.txt")
        path ("base_count.txt")
        path ("overview.txt")

    
    script:
        """
        seqkit stats -T $fastq | tee overview.txt | awk 'BEGIN{fs = "\t"} { sum+=\$4} END{print sum}' > read_count.txt
        cat overview.txt | tr -d \\, | awk 'BEGIN{fs = "\t"} { sum+=\$5} END{print sum}' > base_count.txt
        """
}
