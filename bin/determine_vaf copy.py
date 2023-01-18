#!/usr/bin/env python
import logging
import time
from pathlib import Path

import pysam
from detect_indel import extract_indel_evidence
from detect_snp import extract_snp_evidence
from tqdm import tqdm
from vcf_tools import (
    create_bed_positions,
    create_bed_lines,
    initialize_output_vcf,
    write_vcf_entry,
)

from concurrent.futures import ProcessPoolExecutor, wait, ALL_COMPLETED
from multiprocessing import Pool


def write_variants(
    reference: pysam.FastaFile,
    contig: str,
    pos: int,
    bam_af: pysam.AlignmentFile,
    snp_vcf: Path,
    indel_vcf: Path,
    pileup_depth: int,
    minimum_base_quality: int = 10,
    high_base_quality_cutoff: int = 80,
    end_of_amplicon: bool = False,
):
    """
    Identify SNPs and Indels in a given pileup position.
    Results are written to VCF output files.

    Args:
        reference: Reference genome or sequence, FASTA.
        contig: Reference sequence name.
        pos: Position in the alignment pileup to check for variants.
        bam_af: Read alignments, BAM.
        snp_vcf: Output path for SNPs, VCF.
        indel_vcf: Output path for Indels, VCF.
        pileup_depth: Maximum pileup depth, integer (Default=1_000_000).
        minimum_base_quality: Minimum base quality to consider for
            variant calling.
        high_base_quality_cutoff: Cutoff to calculate ratios on high base
            quality nucleotides only (Default = 80).
        end_of_amplicon: (Not used) Flag to determine if we are reaching
            the end of the amplicon, based on positions in BED file.
    """

    positional_pileup = bam_af.pileup(
        contig,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
        stepper="all",
    )

    # list of jobs TODO: contig:position
    # based on contig:position, do snp + indel
    # multiprocessing.pool, run all positions in parallel
    for pileupcolumn in positional_pileup:

        # positional_pileup only has 1 element
        snv_evidence = extract_snp_evidence(
            pileupcolumn, contig, high_base_quality_cutoff
        )

        # TODO: filter here, write immediately to vcf
        # later merge with indels
        write_vcf_entry(snp_vcf, contig, pos, snv_evidence)

        indel_evidence = extract_indel_evidence(
            pileupcolumn,
            contig,
            reference,
            pos,
            high_base_quality_cutoff,
            end_of_amplicon,
        )
        if indel_evidence[1][0] != ".":
            # if reference allele is not '.', then an indel was found
            write_vcf_entry(indel_vcf, contig, pos, indel_evidence)


def process_pileup_column(
    contig,
    pos,
    bam,
    fasta,
    amplicon_ending=False,
    pileup_depth=1_000_000,
    minimum_base_quality=10,
):
    logging.debug(f"process column {contig} | {pos}")
    bam_af = pysam.AlignmentFile(bam, "r")
    reference = pysam.FastaFile(fasta)

    iteration = 0
    for pu_column in bam_af.pileup(
        contig,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
        stepper="all",
    ):
        snv_evidence = extract_snp_evidence(pu_column, contig, minimum_base_quality)
        indel_evidence = extract_indel_evidence(
            pu_column, contig, reference, pos, minimum_base_quality, amplicon_ending
        )
        if indel_evidence[1][0] == "." and indel_evidence[1][1] == ".":
            indel_evidence = None
        if iteration > 0:
            logging.critical("Single thread worker recieved multiple locations.")
    logging.debug(f"done process column {contig} | {pos}")
    return (contig, pos, snv_evidence, indel_evidence)


def main(
    fasta: Path,
    bed_path: Path,
    bam_path: Path,
    snp_output_path: Path,
    indel_output_path: Path,
    pileup_depth: int = 1_000_000,
    minimum_base_quality=10,
    end_of_amplicon_warn_limit=4,
    threads=16,
):
    """
    Run variant detection over given positions.

    Args:
        fasta: Reference genome or sequence, FASTA.
        bed_file_path: Genomic locations, BED.
        bam_path: Read alignments, BAM.
        snp_output_path: Output path for SNPs, VCF.
        indel_output_path: Output path for Indels, VCF.
        pileup_depth: Maximum pileup depth, integer (DEFAULT=1_000_000).
    """

    logging.debug("started main")
    pileup_parameters = []
    promisses = []

    with ProcessPoolExecutor(max_workers=threads) as executor:
        logging.debug("creating location promisses in ProcessPoolExecutor")
        for bed_file_line in create_bed_lines(bed_path):
            contig = bed_file_line[0]
            contig_region_start = int(bed_file_line[1])
            contig_region_stop = int(bed_file_line[2])

            for pos in range(contig_region_start, contig_region_stop):
                amplicon_ending = (
                    True
                    if contig_region_stop - pos > end_of_amplicon_warn_limit
                    else False
                )
                logging.debug("appending to ProcessPoolExecutor")
                promisses.append(
                    executor.submit(
                        process_pileup_column,
                        contig,
                        pos,
                        bam_path,
                        fasta,
                        amplicon_ending,
                        pileup_depth,
                        minimum_base_quality,
                    )
                )
        logging.debug("Asking for results from ProcessPoolExecutor")
        done, not_done = wait(promisses, return_when=ALL_COMPLETED)

        results = [x.result() for x in promisses]

    bam_af = pysam.AlignmentFile(bam_path, "r")

    snp_vcf = initialize_output_vcf(snp_output_path, bam_af.references)
    indel_vcf = initialize_output_vcf(indel_output_path, bam_af.references)

    for result in results:
        contig = result[0]
        position = result[1]
        snp_results = result[2]
        indel_results = result[3]

        write_vcf_entry(snp_vcf, contig, position, snp_results)
        if indel_results:
            write_vcf_entry(indel_vcf, contig, position, indel_results)

    delay_for_pysam_variantfile = 0.5
    time.sleep(delay_for_pysam_variantfile)
    snp_vcf.close()
    indel_vcf.close()


if __name__ == "__main__":
    import argparse

    dev = True
    if not dev:
        parser = argparse.ArgumentParser(description="")

        parser.add_argument("fasta", type=Path)
        parser.add_argument("bed", type=Path)
        parser.add_argument("bam", type=Path)
        parser.add_argument("snp_vcf_out", type=Path)
        parser.add_argument("indel_vcf_out", type=Path)
        args = parser.parse_args()
        logging.info(args)

        main(args.fasta, args.bed, args.bam, args.snp_vcf_out, args.indel_vcf_out)

    if dev:
        # EGFR
        fasta = Path(
            "/data/references/Homo_sapiens/GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna"
        )
        bed = Path("/data/projects/ROD_1125_variant_improvements/EGFR.bed")
        bam = Path(
            "/data/projects/ROD_1125_variant_improvements/ONT_20221121_EGFR/consensus_aligned/284.taged.bam"
        )
        snp_vcf_out = Path("/data/projects/DAM_0111_vc_multi/testsnp_EGFR.vcf")
        indel_vcf_out = Path("/data/projects/DAM_0111_vc_multi/testindel_EGFR.vcf")
        main(fasta, bed, bam, snp_vcf_out, indel_vcf_out)
