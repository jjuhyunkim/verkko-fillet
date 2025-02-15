import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

# logfile = open('hg002.log', 'w')

def compress_and_build_maps(fasta_file, density):
    '''Compresses the input fasta file and returns two maps

    Args:
        fasta_file (str): The file location of the input fasta
    
    Returns:
        id_to_info (dict): A dictionary mapping a fasta id to a dictionary containing its run length encoding, as well as two dictionaries mapping indices between compressed and uncompressed space
        id_to_compressed_seq (dict): A dictionary mapping a fasta id to its compressed sequence
    '''
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    id_to_info = dict()
    id_to_compressed_seq = dict()
    for fasta in fasta_sequences:
        # logfile.write(fasta.id)
        compressed_seq = ''
        compressed_to_uncompressed = dict()
        uncompressed_to_compressed = dict()
        rle = []
        id, uncompressed_seq = fasta.id, str(fasta.seq.upper())
        uncompressed_seq_len = len(uncompressed_seq)
        uncompressed_idx = 0
        compressed_idx = 0
        run = 1
        while uncompressed_idx < len(uncompressed_seq):
            if uncompressed_idx % density == 0:
                uncompressed_to_compressed[uncompressed_idx] = compressed_idx
            if compressed_idx % density == 0 and compressed_idx not in compressed_to_uncompressed:
                compressed_to_uncompressed[compressed_idx] = uncompressed_idx
            if uncompressed_idx < len(uncompressed_seq) - 1 and uncompressed_seq[uncompressed_idx] == uncompressed_seq[uncompressed_idx + 1]:
                run += 1
                uncompressed_idx += 1
            else:
                compressed_seq += uncompressed_seq[uncompressed_idx]
                rle.append(run)
                run = 1
                uncompressed_idx += 1
                compressed_idx += 1
        id_to_info[id]= dict()
        id_to_compressed_seq[id] = SeqRecord(
            Seq(compressed_seq),
            id=fasta.id,
            name=fasta.name,
            description=fasta.description
        )
        compressed_seq_len = len(compressed_seq)
        id_to_info[id]['compression_map'] = uncompressed_to_compressed
        id_to_info[id]['uncompression_map'] = compressed_to_uncompressed
        id_to_info[id]['rle'] = rle
        id_to_info[id]['uncompressed_len'] = uncompressed_seq_len
        id_to_info[id]['compressed_len'] = compressed_seq_len
    return id_to_info, id_to_compressed_seq

def write_compressed_fasta(id_to_compressed_seq, compressed_fasta_file) :
    '''Write the compressed fasta to file

    Args:
        id_to_compressed_seq (dict)ionary mapping fasta ids to their compressed sequences
        compressed_fasta_file (str): The file location where the output fasta will be written
    
    Returns:
        None
    '''
    f = open(compressed_fasta_file, 'w')
    SeqIO.write(id_to_compressed_seq.values(), f, 'fasta')
    f.close()

def save_map(id_to_info, map_file) :
    '''Write the map to file

    Args:
        id_to_info (dict)ionary mapping fasta ids to a dictionary containing information to be stored in the map for liftover
        map_file (str): The file location where the output map will be written
    
    Returns:
        None
    '''
    j = json.dumps(id_to_info)
    f = open(map_file, 'w')
    f.write(j)
    f.close()

def open_map(map_file) :
    '''Open the map file

    Args:
        map_file (str): The file location where the map is located
    
    Returns:
        map (dict): A dictionary containing the map
    '''
    f = open(map_file, 'r')
    map = json.load(f)
    f.close()
    return map

def write_uncompressed_fasta(in_fasta_file, out_fasta_file, map) :
    '''Uncompress the input fasta and write it to file

    Args:
        in_fasta_file (str): The file location where the input compressed fasta is located
        out_fasta_file (str): The file location where the output uncompressed fasta will be written
        map (dict): A dictionary containing the map used to uncompress the fasta
    
    Returns:
        None
    '''
    out_file = open(out_fasta_file, 'w')
    fasta_sequences = SeqIO.parse(open(in_fasta_file),'fasta')
    id_to_uncompressed_seq = dict()
    for fasta in fasta_sequences:
        uncompressed_seq = ''
        rle = map[fasta.id]['rle']
        compressed_seq = fasta.seq
        for i in range(len(compressed_seq)):
            uncompressed_seq += compressed_seq[i] * rle[i]
        id_to_uncompressed_seq[fasta.id] = SeqRecord(
            Seq(uncompressed_seq),
            id=fasta.id,
            name=fasta.name,
            description=fasta.description
        )
    SeqIO.write(id_to_uncompressed_seq.values(), out_file, 'fasta')
    out_file.close()

def compress_writefasta_savemap(in_fasta_file, out_fasta_file, map_file, density) :
    '''Wrapper function to execute process of compressing the input fasta, writing the compressed fasta, and saving the mapping necessary for liftover

    Args:
        in_fasta_file (str): The file location where the input uncompressed fasta is located
        out_fasta_file (str): The file location where the output compressed fasta will be written
        map_file (str): The file location where the output map will be written
        density (int): The sparseness of the mapping
    
    Returns:
        None
    '''
    id_to_info, id_to_compressed_seq = compress_and_build_maps(in_fasta_file, density)
    write_compressed_fasta(id_to_compressed_seq, out_fasta_file)
    save_map(id_to_info, map_file) 

def uncompress_writefasta(in_fasta_file, out_fasta_file, map_file) :
    '''Wrapper function to execute process of uncompressing the input fasta, and writing the uncompressed fasta

    Args:
        in_fasta_file (str): The file location where the input compressed fasta is located
        out_fasta_file (str): The file location where the output uncompressed fasta will be written
        map_file (str): The file location where the output map will be written
    
    Returns:
        None
    '''
    map = open_map(map_file)
    write_uncompressed_fasta(in_fasta_file, out_fasta_file, map)

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-c', '--compress', help='Take uncompressed reference and compress it', action='store_true')
group.add_argument('-u', '--uncompress', help='Take compressed reference and uncompress it', action='store_true')
parser.add_argument('in_fasta_file', type=str, help='Path to the input fasta file')
parser.add_argument('out_fasta_file', type=str, help='Path to the file where the program generated fasta will be written')
parser.add_argument('map_file', type=str, help='Path to the file where the program generated map converting between compressed and uncompressed space will be written')
parser.add_argument('-d', '--density', nargs='?', type=int, default=1000, help='Density of generated maps, tradeoff between space and eventual lookup time for lifting')

args = parser.parse_args()

if args.compress:
    compress_writefasta_savemap(in_fasta_file=args.in_fasta_file, out_fasta_file=args.out_fasta_file, map_file=args.map_file, density=args.density)
else:
    uncompress_writefasta(in_fasta_file=args.in_fasta_file, out_fasta_file=args.out_fasta_file, map_file=args.map_file)
