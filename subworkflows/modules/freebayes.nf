#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Freebayes {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(input_bam_file), path(input_bai_file), path(reference)
        tuple val(X), path(roi)

    output:
        path("${input_bam_file.SimpleName}.vcf")

    script:
        ref = reference
        """
        freebayes \
        -f $ref \
        -t $roi \
        --min-alternate-fraction $params.freebayes.min_alt_fraction \
        --min-alternate-count $params.freebayes.min_alt_count \
        --min-base-quality $params.freebayes.min_base_qual \
        --vcf ${input_bam_file.SimpleName}.vcf \
        $input_bam_file
        """
}

process SeparateMultiallelicVariants{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        path(vcf)

    output:
        path("${vcf.SimpleName}_norm.vcf")

    script:
        """
        bcftools norm -m -any $vcf > ${vcf.SimpleName}_norm.vcf
        """
}

process FilterFreebayesVariants{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple path(vcf), val(X), path(perbase_table)

    output:
        path("${vcf.simpleName}_filtered.vcf")
    
    script:
        """
        vcf_filter_freebayes.py -i $vcf -o ${vcf.simpleName}_filtered.vcf \
        --perbase_table $perbase_table \
        --dynamic_vaf_params $params.dynamic_vaf_params_file \
        --min_dpq $params.var_filters.min_dpq \
        --min_dpq_n $params.var_filters.min_dpq_n \
        --min_dpq_ratio $params.var_filters.min_dpq_ratio \
        --max_sap $params.var_filters.max_sap \
        --min_rel_ratio $params.var_filters.min_rel_ratio \
        --min_abq $params.var_filters.min_abq        
        """
}