import argparse

import pysam
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# modified version of the RCA processor written by liting to rotate bams to fasta



class RcaProcessor:
    """
    A processor for modifying bam files concerning circular reads generated by Rolling circle amplification.
    """
    def __init__(self, path, write_fasta=None):
        self.path = path
        self.write_fasta = write_fasta
        self.bam = pysam.AlignmentFile(self.path)


    def _find_clipped_regions(self, read):
        left_soft = False
        right_soft = False
        rotate_head_to_tail = 0
        rotate_tail_to_head = 0
         # If a read starts with soft clip or hard clip:
        if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:
            left_soft = True
            rotate_head_to_tail = read.cigartuples[0][1]
        if read.cigartuples[-1][0] == 4 or read.cigartuples[-1][0] == 5:
            right_soft = True
            rotate_tail_to_head = read.cigartuples[-1][1]
        #If both side has clipping, solution now: rotate the longer fragment to the shorter side.
        if left_soft and right_soft:
            # check each side length
            if rotate_tail_to_head - rotate_head_to_tail > 0:
                left_soft = False
            else:
                right_soft = False
        # if the read is originally on the Crick strand (mapped reverse), all information need to swap
        if read.is_reverse:
            right_soft, left_soft = left_soft, right_soft
            rotate_head_to_tail, rotate_tail_to_head = rotate_tail_to_head, rotate_head_to_tail

        return left_soft, right_soft, rotate_head_to_tail, rotate_tail_to_head     

    def rotate_bam_by_cigar(self):
        with open(self.write_fasta, "w+") as fasta:
            for i, read in enumerate(self.bam):
                read_seq = read.get_forward_sequence()
                if read.is_unmapped:
                    fasta.write(f">{read.qname}\n")
                    fasta.write(f"{read_seq}\n")
                else:
                    if read.is_supplementary:
                        continue

                    read_len = len(read_seq)
                    left_soft, right_soft, rotate_head_to_tail, rotate_tail_to_head  = self._find_clipped_regions(read)
                    
                    if right_soft:
                        rotated_seq = read_seq[-rotate_tail_to_head:] + read_seq[:read_len - rotate_tail_to_head]
                    if left_soft:
                        rotated_seq = read_seq[-(read_len - rotate_head_to_tail):] + read_seq[:rotate_head_to_tail]
                
                    # write results:
                    fasta.write(f">{read.qname}\n")
                    fasta.write(f">{rotated_seq}\n")



class BamToFastqRotator:
    """
    A processor for modifying bam files concerning circular reads generated by Rolling circle amplification.
    """
    def __init__(self, path, write_fastq=None):
        self.path = path
        self.write_fastq = write_fastq
        self.bam = pysam.AlignmentFile(self.path)


    def _find_clipped_regions(self, read):
        left_soft = False
        right_soft = False
        rotate_head_to_tail = 0
        rotate_tail_to_head = 0
         # If a read starts with soft clip or hard clip:
        if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:
            left_soft = True
            rotate_head_to_tail = read.cigartuples[0][1]
        if read.cigartuples[-1][0] == 4 or read.cigartuples[-1][0] == 5:
            right_soft = True
            rotate_tail_to_head = read.cigartuples[-1][1]
        #If both side has clipping, solution now: rotate the longer fragment to the shorter side.
        if left_soft and right_soft:
            # check each side length
            if rotate_tail_to_head - rotate_head_to_tail > 0:
                left_soft = False
            else:
                right_soft = False
        # if the read is originally on the Crick strand (mapped reverse), all information need to swap
        if read.is_reverse:
            right_soft, left_soft = left_soft, right_soft
            rotate_head_to_tail, rotate_tail_to_head = rotate_tail_to_head, rotate_head_to_tail

        return left_soft, right_soft, rotate_head_to_tail, rotate_tail_to_head     

    def rotate_bam_by_cigar(self):
        with open(self.write_fastq, "w") as fastq_handle:
            for i, read in enumerate(self.bam):
                # read is a pysam alignedSegment, not a biopython object
                read_seq = read.get_forward_sequence()
                read_qual = read.qual
                if read.is_unmapped:
                    fastq_read = SeqRecord(
                        Seq(read_seq),
                        id = "", 
                        name =  read.qname,
                        description= "",
                    )
                    rotated_qual = [ord(x) - 33 for x in read_qual]
                    fastq_read.letter_annotations["phred_quality"] = rotated_qual
                    SeqIO.write(fastq_read, handle=fastq_handle, format = "fastq")

                else:
                    if read.is_supplementary or read.is_secondary:
                        continue

                    read_len = len(read_seq)
                    left_soft, right_soft, rotate_head_to_tail, rotate_tail_to_head  = self._find_clipped_regions(read)
                    
                    # rotate based on softclipping findings
                    if right_soft:
                        rotated_seq = read_seq[-rotate_tail_to_head:] + read_seq[:read_len - rotate_tail_to_head]
                        rotated_qual = read_qual[-rotate_tail_to_head:] + read_qual[:read_len - rotate_tail_to_head]
                    elif left_soft:
                        rotated_seq = read_seq[-(read_len - rotate_head_to_tail):] + read_seq[:rotate_head_to_tail]
                        rotated_qual = read_qual[-(read_len - rotate_head_to_tail):] + read_qual[:rotate_head_to_tail]
                    else:
                        rotated_seq = read_seq
                        rotated_qual = read_qual
                    
                    if len(read_seq) != len(rotated_seq):
                        raise IOError("mistake in rotating algo")
                    if len(read_qual) != len(rotated_qual):
                        raise IOError("mistake in rotating algo with quals")

                    if read.is_reverse:
                        # reverse the forward sequence to get back to reverse data
                        final_sequence = Seq(rotated_seq[::-1])
                        rotated_qual = rotated_qual[::-1]
                    else:
                        final_sequence = Seq(rotated_seq)
                    
                    fastq_read = SeqRecord(
                        final_sequence,
                        id = read.qname, 
                        name =  "somename",
                        description= "some descr"
                    )
                    rotated_qual = [ord(x) - 33 for x in rotated_qual]
                    fastq_read.letter_annotations["phred_quality"] = rotated_qual

                    SeqIO.write(sequences=fastq_read, handle=fastq_handle, format="fastq")
                    

class BamToFastaRotator:
    """
    A processor for modifying bam files concerning circular reads generated by Rolling circle amplification.
    """
    def __init__(self, path, write_fasta=None):
        self.path = path
        self.write_fasta = write_fasta
        self.bam = pysam.AlignmentFile(self.path)


    def _find_clipped_regions(self, read):
        left_soft = False
        right_soft = False
        rotate_head_to_tail = 0
        rotate_tail_to_head = 0
         # If a read starts with soft clip or hard clip:
        if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:
            left_soft = True
            rotate_head_to_tail = read.cigartuples[0][1]
        if read.cigartuples[-1][0] == 4 or read.cigartuples[-1][0] == 5:
            right_soft = True
            rotate_tail_to_head = read.cigartuples[-1][1]
        #If both side has clipping, solution now: rotate the longer fragment to the shorter side.
        if left_soft and right_soft:
            # check each side length
            if rotate_tail_to_head - rotate_head_to_tail > 0:
                left_soft = False
            else:
                right_soft = False
        # if the read is originally on the Crick strand (mapped reverse), all information need to swap
        if read.is_reverse:
            right_soft, left_soft = left_soft, right_soft
            rotate_head_to_tail, rotate_tail_to_head = rotate_tail_to_head, rotate_head_to_tail

        return left_soft, right_soft, rotate_head_to_tail, rotate_tail_to_head     

    def rotate_bam_by_cigar(self):
        with open(self.write_fasta, "w") as fasta_handle:
            for i, read in enumerate(self.bam):
                # read is a pysam alignedSegment, not a biopython object
                read_seq = read.get_forward_sequence()
                if read.is_unmapped:
                    fasta_read = SeqRecord(
                        Seq(read_seq),
                        id = "", 
                        name =  read.qname,
                        description= "",
                    )
                    SeqIO.write(fasta_read, handle=fasta_handle, format = "fasta")

                else:
                    if read.is_supplementary or read.is_secondary:
                        continue

                    read_len = len(read_seq)
                    left_soft, right_soft, rotate_head_to_tail, rotate_tail_to_head  = self._find_clipped_regions(read)
                    
                    # rotate based on softclipping findings
                    if right_soft:
                        rotated_seq = read_seq[-rotate_tail_to_head:] + read_seq[:read_len - rotate_tail_to_head]
                    elif left_soft:
                        rotated_seq = read_seq[-(read_len - rotate_head_to_tail):] + read_seq[:rotate_head_to_tail]
                    else:
                        rotated_seq = read_seq
                    
                    if len(read_seq) != len(rotated_seq):
                        raise IOError("mistake in rotating algo")

                    if read.is_reverse:
                        # reverse the forward sequence to get back to reverse data
                        final_sequence = Seq(rotated_seq[::-1])
                    else:
                        final_sequence = Seq(rotated_seq)
                    
                    fasta_read = SeqRecord(
                        final_sequence,
                        id = read.qname, 
                        name =  "somename",
                        description= "some descr"
                    )

                    SeqIO.write(sequences=fasta_read, handle=fasta_handle, format="fasta")
       

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Rotate read records based on bam mapping to completely map circular reads. Output in fasta or fastq.')
    parser.add_argument('-i','--input_bam', type=str,
                        help='input bam file path')
    parser.add_argument('-o','--output_file', type=str,
                        help='output file path')
    parser.add_argument('-t','--output_type', type=str,
                        help='output fasta file')
                        
    args = parser.parse_args()
    output_type = str(args.output_type)
    if output_type == "fastq":
        Rotator = BamToFastqRotator
    elif output_type == "fasta":
        Rotator = BamToFastaRotator
    else:
        raise IOError(f"output_type: {output_type} not supported")
    
    bam = Rotator(args.input_bam, args.output_file)
    bam.rotate_bam_by_cigar()
