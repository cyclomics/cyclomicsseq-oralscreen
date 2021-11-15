
nextflow.enable.dsl=2

process Tidehunter{
    // _tide_consensus.fasta in ConCall
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'tidehunter'
    
    input:
        path fasta

    output:
        path "${fasta.SimpleName}_tide_consensus.fasta", emit: consensus

    script:
        """
        TideHunter -t ${task.cpus} $fasta > ${fasta.SimpleName}_tide_consensus.fasta
        """
}

process TidehunterFullLength{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'tidehunter'
    
    input:
        path fasta
        path prime_3
        path prime_5

    output:
        path "${fasta.SimpleName}_tide_consensus_full_length.fasta", emit: consensus

    script:
        """
        TideHunter -t ${task.cpus} -p 20 -a 0.7 $fasta > ${fasta.SimpleName}_tide_consensus_full_length.fasta
        """
}

process Tidehunter53{
    // Parameterized, but unused right now
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'tidehunter'
    
    input:
        path fasta
        path prime_3
        path prime_5

    output:
        path "${fasta.SimpleName}.consensus.tsv", emit: consensus

    script:
        """
        echo "$params.tidehunter.headerlines" > ${fasta.SimpleName}.consensus.tsv
        TideHunter -f 2 -t ${task.cpus} -5 $prime_5 -3 $prime_3 -p $params.tidehunter.minimum_period -a $params.tidehunter.minimum_match_ratio $fasta >> ${fasta.SimpleName}.consensus.tsv
        """
}

process TideHunterTrimmmer {
    // Plug tidehunter fasta into this
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        path fasta_full_length
        path fasta_all
    
    output:
        path "${original_fasta.SimpleName}.full_length.fasta", emit: fasta_full_length 

    script:
        """
        sed 's/,.*//' $original_fasta > ${original_fasta.SimpleName}.full_length.fasta
        """
}   
