
nextflow.enable.dsl=2

include {
    PYCOQC
} from "./modules/pycoqc.nf"

include {
    MinionQc
    MinionQcToJson
} from "./modules/minion_qc.nf"

include {
    CollectClassificationTypes
} from "./modules/bin"

include {
    CountFastqInfo as FastqInfoRaw
    CountFastqInfo as FastqInfoConsensus
} from "./modules/seqkit"

include {
    SamtoolsQuickcheck
    SamtoolsFlagstats
} from "./modules/samtools"

include {
    FindRegionOfInterest
} from "./modules/mosdepth"

workflow  QC_pycoqc{
    take:
        read_directory
    main:
        PYCOQC(read_directory)
    emit:
        PYCOQC.output
}

workflow QC_MinionQc {
    take:
        read_directory
    main:
        MinionQc(read_directory)
        MinionQcToJson(MinionQc.out.summary)
    emit:
        MinionQcToJson.output
}


workflow PostQC {
    take:
        read_info
        fastq_raw
        fastq_consensus
        consensus_bam
    main:
        // dont add the ID to the process
        // CollectClassificationTypes(read_info.map(it -> it[1]).collect())
        FastqInfoRaw(fastq_raw.collect())
        FastqInfoConsensus(fastq_consensus.map(it -> it[1]).collect())

        SamtoolsQuickcheck(consensus_bam)
        SamtoolsQuickcheck.out.view()
        SamtoolsFlagstats(consensus_bam)
        SamtoolsFlagstats.out.view()
        
        FindRegionOfInterest(consensus_bam)
}
