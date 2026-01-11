from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os
import re
import pandas as pd


def split_fasta_at_ns(input_fasta, output_prefix):
    """
    Splits a FASTA file at any sequence of 'N's into multiple records
    and writes them all into one output FASTA.
    """

    output_file = f"{output_prefix}_split.fasta"
    
    if os.path.exists(output_file):
        print(f"Output file {output_file} already exists. Skipping split.")
        return
    with open(input_fasta, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # Use regex to keep Ns as separators
            segments = re.split(r"(N+)", str(record.seq))
            seq_index = 0
            for segment in segments:
                if len(segment) == 0:
                    continue
                seq_index += 1
                new_id = f"{record.id}_seg{seq_index}"
                new_record = SeqRecord(Seq(segment), id=new_id, description=f"Original: {record.description}")
                SeqIO.write(new_record, outfile, "fasta")
    
    print(f"Written split sequences to {output_file}")
    return record.id

def getContigOrder(input_split_fasta, gapOrder):
    
    splitFa = SeqIO.to_dict(SeqIO.parse(input_split_fasta, "fasta"))
    contiglist = list(splitFa.keys())
    print(f"Original contig list: {contiglist}")

    idx = 2 * gapOrder - 1
    if 0 <= idx < len(contiglist):
        contiglist[idx] = ["newGap", "rDNA_consensus", "rDNA_consensus", "newGap"]
    else:
        print(f"Warning: gapOrder {gapOrder} is out of range for contiglist of length {len(contiglist)}")

    contiglist = [y for x in contiglist for y in (x if isinstance(x, list) else [x])]
    
    print(f"Modified contig list: {contiglist}")
    return contiglist, splitFa
def concatFa(rDNA_morph_fasta,gapSize,contiglist,output_fasta, splitFa, contig):# rDNA_consensus
    if not os.path.exists(rDNA_morph_fasta):
        print(f"Error: rDNA consensus file {rDNA_morph_fasta} does not exist.")
        sys.exit(1)
    if os.path.exists(output_fasta):
        print(f"Output file {output_fasta} already exists. Skipping concatenation.")
        return
    # read the rDNA consensus
    rDNA_consensus = SeqIO.to_dict(SeqIO.parse(rDNA_morph_fasta, "fasta"))
    rDNA_consensus_name = list(rDNA_consensus.keys())[0]
    # gap seq
    newGap = ("N" * gapSize)
    # concat the sequences
    new_seq_str = ""

    for newcontig in contiglist:
        if newcontig == "newGap":
            seq_to_add = newGap
        elif newcontig == "rDNA_consensus":
            seq_to_add = str(rDNA_consensus[rDNA_consensus_name].seq)
        else:
            seq_to_add = str(splitFa[newcontig].seq)
        
        new_seq_str += seq_to_add

    # Wrap string in Seq object
    new_record = SeqRecord(
        Seq(new_seq_str),
        id=f"{contig}_rDNA_patched",
        description=f"{contig} with rDNA and gap inserted(size={gapSize})"
    )
    # Write the fasta
    SeqIO.write(new_record, output_fasta, "fasta")
    
    print(f"Written patched sequence to {output_fasta}")

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python insert_rDNAcns.py <input.fasta> <output.fasta> <rDNA_cns.fasta> <size_of_gap> <gapOrder> <force> <tmp_file>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    rDNA_morph_fasta = sys.argv[3]
    gapSize = int(sys.argv[4])
    gapOrder = int(sys.argv[5])
    force = sys.argv[6].lower() == "true" if len(sys.argv) > 6 else False
    tmp_file = sys.argv[7] if len(sys.argv) > 7 else 1

    # check the inputs 
    print(f"Checking inputs:",
          f"Input fasta: {input_fasta}"
          f"\nOutput fasta: {output_fasta}"
          f"\nrDNA morph fasta: {rDNA_morph_fasta}"
          f"\nGap size: {gapSize}"
          f"\nGap order: {gapOrder}"
          f"\nForce overwrite: {force}"
          f"\nTemporary file prefix: {tmp_file if tmp_file != 1 else input_fasta + '.tmp'}\n")

    # check if output file exists
    if force == True:
        if os.path.exists(f"{tmp_file}_split.fasta"):
            os.remove(f"{tmp_file}_split.fasta")
        if os.path.exists(output_fasta):
            os.remove(output_fasta)

    if tmp_file == 1:
        tmp_file = input_fasta + ".tmp"

    # split contig
    contig = split_fasta_at_ns(input_fasta, tmp_file)
    input_split_fasta = f"{tmp_file}_split.fasta"
    # fix the orders
    contiglist, splitFa = getContigOrder(input_split_fasta, gapOrder)
    # concat the sequences
    concatFa(rDNA_morph_fasta, gapSize, contiglist, output_fasta, splitFa, contig)
    # clean up
    os.remove(input_split_fasta)
    