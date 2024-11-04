#!/usr/bin/env python
import io
import json
from collections import OrderedDict
from pathlib import Path
from typing import Any, Tuple

import pandas as pd
import requests


class NoColocatedVariantsException(Exception):
    "Raised when an Ensembl-VEP query does not return any colocated variants."

    pass


def obtain_legacy_cosmic_id(
    cosv_id: str,
    base_url: str = "https://clinicaltables.nlm.nih.gov/api/cosmic/v4/search?terms=",
) -> str:
    """extract legacy information from the clinicaltables.nlm.nih.gov api."""

    query = base_url + cosv_id + "&ef=LegacyMutationID"
    response = requests.get(url=query, headers={"Content-Type": "application/json"})
    try:
        id_list = set(json.loads(response.text)[2]["LegacyMutationID"])
    except KeyError:  # LegacyMutationID was not found in
        return "None"

    legacy_ids = [id for id in id_list if "COSM" in id]
    return ",".join(legacy_ids)


class PanelToHg38Translator:
    def __init__(self, vcf):
        self.vcf = vcf

    def find_amplicon_loci(self, contigs):
        loci = contigs.str.split("[_|-]", expand=False)
        chromosome = [str(locus[-3]) for locus in loci]
        start_pos = [int(locus[-2]) for locus in loci]

        return chromosome, start_pos

    def translate(self):
        chromosome, amplicon_start = self.find_amplicon_loci(self.vcf.CHROM)

        self.vcf.CHROM = chromosome
        self.vcf.POS = self.vcf.POS + amplicon_start

        return self.vcf


class VcfFile:
    def __init__(self, vcf_file, translator=None):
        # self.vcf_file = vcf_file
        self.vcf_header = ""
        self.vcf = self.read_vcf(vcf_file)
        self.translator = translator(self.vcf)

        # def annotate(self):
        if self.translator:
            self.vcf = self.translator.translate()

        # TODO: write new VCF with GRCh38 positions
        # self.write()

    def relaxed_float(self, x: Any) -> float:
        """Return a float, with value error catch"""
        try:
            my_float = float(x)
        except ValueError:
            my_float = float(0)
        return my_float

    def read_vcf(self, path: Path) -> pd.DataFrame:
        """Read in a VCF file and return it as Pandas DataFrame"""
        with open(path, "r") as f:
            header = []
            lines = []
            for l in f:
                if l.startswith("##"):
                    header.append(l)
                else:
                    lines.append(l)

        self.vcf_header = header
        df = pd.read_csv(
            io.StringIO("".join(lines)),
            dtype={
                "#CHROM": str,
                "POS": int,
                "ID": str,
                "REF": str,
                "ALT": str,
                "QUAL": str,
                "FILTER": str,
                "INFO": str,
            },
            sep="\t",
        ).rename(columns={"#CHROM": "CHROM"})

        if not df.empty:
            df.columns = df.columns.str.upper()
            formats = df.FORMAT[0].split(":")
            for i, fmt in enumerate(formats):
                df[fmt] = df.SAMPLE1.apply(
                    lambda x: self.relaxed_float(x.split(":")[i])
                    if (x.split(":")[i])
                    else 0
                )

        return df

    def write(self, path: Path):
        """Write output VCF file"""
        with open(path, "w") as new_vcf:
            new_vcf.writelines(self.vcf_header)

            writeable_vcf = self.vcf.rename(columns={"CHROM": "#CHROM"})
            writeable_vcf = writeable_vcf[
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    "SAMPLE1",
                ]
            ]

            new_vcf.writelines(writeable_vcf.to_csv(sep="\t", index=False))

    @staticmethod
    def get_query_allele_positions(
        chromosome: str, start: int, ref_allele: str, alt_allele: str
    ) -> Tuple[int, int]:
        """Determine start and end positions for Ensembl query

        Positions are determined based on reference and alternative
        alleles, which inform if a variant is a SNP, an 1 or a
        deletion. Index adjustments have to be made for indels.
        """

        if len(ref_allele) > len(alt_allele):
            # this is a deletion
            # remove first base of ref which is maintained and not part of the deletion
            end = start + len(ref_allele) - 1
            start += 1
            alt_allele = "-"

        elif len(ref_allele) < len(alt_allele):
            # This is an insertion
            end = start
            # remove first base of ref which is maintained and not part of the deletion
            start += 1
            ref_allele = "-"

        elif len(ref_allele) == len(alt_allele):
            # This is a snp (ignoring insdel events)
            end = start

        # Ensembl-VEP needs REF to be the strand, which is always forward (1)
        ref_allele = "1"

        return (start, end, ref_allele, alt_allele)

    @staticmethod
    def ensembl_vep(query: str) -> dict:
        """Query Ensembl-VEP through API to annotate variants"""
        try:
            response = requests.get(
                url=query, headers={"Content-Type": "application/json"}
            )
            json_data = json.loads(response.text)[0]

        except (KeyError, ConnectionError) as error:
            # No response data found in this query location
            # or connection was aborted
            print(error)
            return

        return json_data

    @staticmethod
    def get_annotation_text(vep_json: dict) -> str:
        """Parse VEP response dict and find relevant annotations"""
        # Initialize variant annotations
        variant_class = None
        consequence = None
        mutation_ids = None
        gene = None
        impact = None
        biotype = None
        amino_acids = None
        canonical = None
        sift = None
        polyphen = None
        refseq_transcripts = None
        transcript_id = None
        hgvsc = None
        hgvsp = None

        if not vep_json:
            # Ensembl-VEP query did not return any response
            # e.g. because the variant was in a backbone sequence
            annot_text = "."
            return annot_text

        # Find relevant annotations in response dict
        variant_class = vep_json.get("variant_class")
        consequence = vep_json.get("most_severe_consequence")

        # Try-Except block: cannot always find colocated variants
        # and thus COSMIC IDs
        try:
            colocated_variants = vep_json.get("colocated_variants")

            if not colocated_variants:
                raise NoColocatedVariantsException(
                    "Ensembl-VEP query does not return any colocated variants."
                )

            else:
                mutation_ids = []
                cosmic_legacy_ids = []
                for xref in colocated_variants:
                    if xref["allele_string"] == "COSMIC_MUTATION":
                        cosv = xref["id"]
                        if cosv:
                            mutation_ids.append(cosv)
                            cosm = obtain_legacy_cosmic_id(cosv)
                            cosmic_legacy_ids.append(cosm)
                        else:
                            continue

                # Join list of found IDs into comma-separated string
                mutation_ids = ",".join(mutation_ids) if mutation_ids else "None"
                cosmic_legacy_ids = (
                    ",".join(cosmic_legacy_ids) if cosmic_legacy_ids else "None"
                )

        except NoColocatedVariantsException:
            mutation_ids = "None"
            cosmic_legacy_ids = "None"

        transcript_cons = vep_json.get("transcript_consequences")
        # Transcript consequences can differ a lot per query
        # If something is not found, will be returned as 'None'
        if transcript_cons:
            gene = transcript_cons[0].get("gene_symbol")
            impact = transcript_cons[0].get("impact")
            biotype = transcript_cons[0].get("biotype")
            amino_acids = transcript_cons[0].get("amino_acids")
            canonical = transcript_cons[0].get("canonical")

            sift_prediction = transcript_cons[0].get("sift_prediction")
            sift_score = transcript_cons[0].get("sift_score")
            if sift_prediction:
                sift = f"{sift_prediction}({sift_score})"

            polyphen_prediction = transcript_cons[0].get("polyphen_prediction")
            polyphen_score = transcript_cons[0].get("polyphen_score")
            if polyphen_prediction:
                polyphen = f"{polyphen_prediction}({polyphen_score})"

            refseq_transcripts = transcript_cons[0].get("refseq_transcript_ids")
            refseq_transcript_ids = []
            if refseq_transcripts:
                for id in refseq_transcripts:
                    refseq_transcript_ids.append(id)

            refseq_transcript_ids = (
                ",".join(refseq_transcript_ids) if refseq_transcript_ids else "None"
            )

            transcript_id = transcript_cons[0].get("transcript_id")
            hgvsc = transcript_cons[0].get("hgvsc")
            hgvsp = transcript_cons[0].get("hgvsp")

        # Merge all annotations into a string to be returned
        annot_dict = OrderedDict(
            {
                "variant_class": variant_class,
                "consequence": consequence,
                "COSMIC": mutation_ids,
                "COSMIC legacy": cosmic_legacy_ids,
                "gene": gene,
                "impact": impact,
                "biotype": biotype,
                "amino_acids": amino_acids,
                "canonical": canonical,
                "SIFT": sift,
                "PolyPhen": polyphen,
                "RefSeq transcripts": refseq_transcript_ids,
                "Transcript": transcript_id,
                "HGVSC": hgvsc,
                "HGVSP": hgvsp,
            }
        )

        annot_text = "ANNOTATION;"
        annot_text += ";".join([f"{k}={v}" for k, v in annot_dict.items()])
        # If annotation is None, don't print it
        # annot_text = ";".join([f"{k}={v}" for k, v in annot_dict.items() if v])

        return annot_text

    def annotate_vep(self, server: str):
        """Annotate a set of variants from a VCF file with Ensembl-VEP

        Input: Server URL (string), e.g. "https://rest.ensembl.org"
        """
        # Set API search options
        params = "canonical=1&variant_class=1&hgvs=1&xref_refseq=1&vcf_string=1&pick=1"

        if self.vcf.empty:
            # There are no variants to annotate
            return

        annotations = []
        # Loop over variants in VCF file, annotate one at a time
        # TODO: Parallelize
        for var in self.vcf.iterrows():
            chromosome = var[1]["CHROM"]
            start = var[1]["POS"]
            ref_allele = var[1]["REF"]
            alt_allele = var[1]["ALT"]

            start, end, ref_allele, alt_allele = self.get_query_allele_positions(
                chromosome, start, ref_allele, alt_allele
            )

            query = (
                f"{server}/vep/human/region/"
                f"{chromosome}:{start}-{end}:"
                f"{ref_allele}/{alt_allele}?"
                f"{params}?"
            )

            # Query Ensembl with the VEP API, returns a JSON dict
            vep_json = self.ensembl_vep(query)
            # Parse JSON response to get annotations in a string
            annotation_text = self.get_annotation_text(vep_json)
            # Add to annotations list
            annotations.append(annotation_text)

        # Write annotations to INFO column in VCF output file
        self.vcf["INFO"] = annotations


if __name__ == "__main__":
    dev = False
    if not dev:
        import argparse

        parser = argparse.ArgumentParser(description="Filter a vcf")

        parser.add_argument("variant_vcf", type=Path)
        parser.add_argument("file_out", type=Path)
        args = parser.parse_args()

        # Can be added to argparse
        server = "https://rest.ensembl.org"

        vcf = VcfFile(args.variant_vcf, PanelToHg38Translator)
        vcf.annotate_vep(server)
        vcf.write(args.file_out)

    if dev:
        variant_vcf = "/scratch/spellbook/ROD/cauldron/cyclomicsseq-oralscreen/testing/FAY73116_filtered.vcf"
        server = "https://rest.ensembl.org"

        vcf = VcfFile(variant_vcf, PanelToHg38Translator)
        vcf.annotate_vep(server)
        vcf.write(
            "/scratch/spellbook/ROD/cauldron/cyclomicsseq-oralscreen/testing/FAY73116_filtered_testannotated.vcf"
        )
