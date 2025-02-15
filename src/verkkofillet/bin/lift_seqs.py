import argparse
from Bio import SeqIO
import json


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

def get_closest_value(arr, val):
    '''Get closest value to val in an unsorted list

    Args:
        arr (list): The list containing all values to check
        val (int): The value of interest
    
    Returns:
        (int): The element of arr closest to val
    '''
    return arr[min(range(len(arr)), key = lambda i: abs(arr[i]-val))]

def get_closest_value_sorted(arr, val):
    '''Get closest value to val in a sorted list

    Args:
        arr (list): The list containing all values to check
        val (int): The value of interest
    
    Returns:
        (int): The element of arr closest to val
    '''
    n = len(arr)
    if (val <= arr[0]):
        return arr[0]
    if (val >= arr[n - 1]):
        return arr[n - 1]

    left = 0
    right = n
    mid = 0
    while (left < right):
        mid = (left + right) // 2
        if (arr[mid] == val):
            return arr[mid]
        if (val < arr[mid]) :
            if (mid > 0 and val > arr[mid - 1]):
                return pick_closer_of_two(arr[mid - 1], arr[mid], val)
            right = mid
        else :
            if (mid < n - 1 and val < arr[mid + 1]):
                return pick_closer_of_two(arr[mid], arr[mid + 1], val) 
            left = mid + 1
    return arr[mid]
 
def pick_closer_of_two(a, b, val):
    if abs(a - val) <= abs(b - val):
        return a
    else:
        return b

def get_end_of_hp(seq, idx):
    char = seq[idx]
    i = idx
    while seq[i] == char:
        i += 1
    return i - 1

def get_start_of_hp(seq, idx):
    char = seq[idx]
    i = idx
    while seq[i] == char:
        i -= 1
    return i + 1

def hps_between(seq, idx1, idx2):
    start = get_end_of_hp(seq, min(idx1, idx2))
    end = get_start_of_hp(seq, max(idx1, idx2))
    curr_char = seq[start]
    hps = 0
    for i in range(start, end):
        if seq[i] != curr_char:
            hps += 1
        curr_char = seq[i]
    return hps

def same_hp(seq, idx1, idx2):
    return seq[idx1] == seq[idx2] and hps_between(seq, idx1, idx2) == 0

def lift_c_to_u_from_smaller_idx(compressed_coord, curr_compressed_coord, map, rle, is_end):
    '''A helper function for lifting a coordinate from compressed to uncompressed space in the case when the nearest map key pertains to a smaller compressed index than the one being lifted

    Args:
        compressed_coord (int): The compressed coordinate being lifted
        curr_compressed_coord (int): The nearest compressed coordinate present in the map
        map (dict): A dictionary mapping compresed indices to their corresponding uncompressed indices
        rle (list): A run length encoding of the uncompressed sequence
        is_end (bool): True if this corresponds to the end index of the sequence being lifted
    
    Returns:
        curr_uncompressed_coord (int): The corresponding coordinate in uncompressed space
    '''
    curr_uncompressed_coord = map[curr_compressed_coord]
    while curr_compressed_coord < compressed_coord:
        curr_uncompressed_coord += rle[curr_compressed_coord]
        curr_compressed_coord += 1
    if is_end:
        curr_uncompressed_coord += rle[curr_compressed_coord]
    return curr_uncompressed_coord

def lift_c_to_u_from_larger_idx(compressed_coord, curr_compressed_coord, map, rle, is_end):
    '''A helper function for lifting a coordinate from compressed to uncompressed space in the case when the nearest map key pertains to a larger compressed index than the one being lifted

    Args:
        compressed_coord (int): The compressed coordinate being lifted
        curr_compressed_coord (int): The nearest compressed coordinate present in the map
        map (dict): A dictionary mapping compresed indices to their corresponding uncompressed indices
        rle (list): A run length encoding of the uncompressed sequence
        is_end (bool): True if this corresponds to the end index of the sequence being lifted
    
    Returns:
        curr_uncompressed_coord (int): The corresponding coordinate in uncompressed space
    '''
    curr_uncompressed_coord = map[curr_compressed_coord]
    while curr_compressed_coord > compressed_coord:
        curr_uncompressed_coord -= rle[curr_compressed_coord - 1]
        curr_compressed_coord -= 1
    if is_end:
        curr_uncompressed_coord += rle[curr_compressed_coord]
    return curr_uncompressed_coord

def lift_c_to_u_helper(compressed_coord, curr_compressed_coord, map, rle, is_end):
    '''A helper function for lifting a coordinate from compressed to uncompressed space

    Args:
        compressed_coord (int): The compressed coordinate being lifted
        curr_compressed_coord (int): The nearest compressed coordinate present in the map
        map (dict): A dictionary mapping compresed indices to their corresponding uncompressed indices
        rle (list): A run length encoding of the uncompressed sequence
        is_end (bool): True if this corresponds to the end index of the sequence being lifted
    
    Returns:
        (int): The corresponding coordinate in uncompressed space
    '''
    if compressed_coord == curr_compressed_coord:
        if is_end:
            return map[curr_compressed_coord] + rle[curr_compressed_coord]
        else:
            return map[curr_compressed_coord]
    elif curr_compressed_coord < compressed_coord:
        return lift_c_to_u_from_smaller_idx(compressed_coord, curr_compressed_coord, map, rle, is_end)
    else:
        return lift_c_to_u_from_larger_idx(compressed_coord, curr_compressed_coord, map, rle, is_end)

def lift_c_to_u(compressed_start, compressed_end, id, map_file):
    '''The main function for lifting 0-based [start, end) coordinates from compressed to uncompressed space

    Args:
        compressed_start (int): The 0-based, inclusive compressed start coordinate being lifted
        compressed_end (int): The 0-based, exclsuive compressed end coordinate being lifted
        id (str): The fasta id of the sequence being lifted
        map_file (str): A path to the file containing the map needed for liftover
    
    Returns:
        uncompressed_start, uncompressed_end (int, int): The corresponding coordinates in uncompressed space, in the form [uncompressed_start, uncompressed_end)
    '''
    map = open_map(map_file)
    rle = map[id]['rle']
    uncompression_map = map[id]['uncompression_map']
    print(rle)
    print(uncompression_map)
    d = {int(k):int(v) for k,v in uncompression_map.items()}
    keys = list(d.keys())
    keys.sort()
    curr_start = get_closest_value_sorted(keys, compressed_start)
    curr_end = get_closest_value_sorted(keys, compressed_end-1)
    uncompressed_start = lift_c_to_u_helper(compressed_start, curr_start, d, rle, False)
    uncompressed_end = lift_c_to_u_helper(compressed_end-1, curr_end, d, rle, True)
    return uncompressed_start, uncompressed_end

def lift_u_to_c_from_smaller_idx(uncompressed_coord, curr_uncompressed_coord, uncompressed_seq, map, rle, is_end):
    '''A helper function for lifting a coordinate from uncompressed to compressed space in the case when the nearest map key pertains to a smaller uncompressed index than the one being lifted

    Args:
        uncompressed_coord (int): The uncompressed coordinate being lifted
        curr_uncompressed_coord (int): The nearest uncompressed coordinate present in the map
        uncompressed_seq (str): The uncompressed sequence being lifted
        map (dict): A dictionary mapping uncompresed indices to their corresponding compressed indices
        rle (list): A run length encoding of the uncompressed sequence
        is_end (bool): True if this corresponds to the end index of the sequence being lifted
    
    Returns:
        curr_compressed_coord (int): The corresponding coordinate in compressed space
    '''
    curr_compressed_coord = map[curr_uncompressed_coord]
    in_same_hp = same_hp(uncompressed_seq, uncompressed_coord, curr_uncompressed_coord)
    if in_same_hp:
        if is_end:
            return curr_compressed_coord + 1
        else:
            return curr_compressed_coord
    else:
        hps = hps_between(uncompressed_seq, uncompressed_coord, curr_uncompressed_coord) + 1
        curr_compressed_coord += hps
    
    if is_end and not same_hp(uncompressed_seq, uncompressed_coord, curr_uncompressed_coord):
        return curr_compressed_coord + 1
    else:
        return curr_compressed_coord

def lift_u_to_c_from_larger_idx(uncompressed_coord, curr_uncompressed_coord, uncompressed_seq, map, rle, is_end):
    '''A helper function for lifting a coordinate from uncompressed to compressed space in the case when the nearest map key pertains to a larger uncompressed index than the one being lifted

    Args:
        uncompressed_coord (int): The uncompressed coordinate being lifted
        curr_uncompressed_coord (int): The nearest uncompressed coordinate present in the map
        uncompressed_seq (str): The uncompressed sequence being lifted
        map (dict): A dictionary mapping uncompresed indices to their corresponding compressed indices
        rle (list): A run length encoding of the uncompressed sequence
        is_end (bool): True if this corresponds to the end index of the sequence being lifted
    
    Returns:
        curr_compressed_coord (int): The corresponding coordinate in compressed space
    '''
    curr_compressed_coord = map[curr_uncompressed_coord]
    in_same_hp = same_hp(uncompressed_seq, uncompressed_coord, curr_uncompressed_coord)
    if in_same_hp:
        if is_end:
            return curr_compressed_coord + 1
        else:
            return curr_compressed_coord
    else:
        hps = hps_between(uncompressed_seq, uncompressed_coord, curr_uncompressed_coord) + 1
        curr_compressed_coord -= hps

    if is_end and not same_hp(uncompressed_seq, uncompressed_coord, curr_uncompressed_coord):
        return curr_compressed_coord + 1
    else:
        return curr_compressed_coord

def lift_u_to_c_helper(uncompressed_coord, curr_uncompressed_coord, uncompressed_seq, map, rle, is_end):
    '''A helper function for lifting a coordinate from uncompressed to compressed space

    Args:
        uncompressed_coord (int): The uncompressed coordinate being lifted
        curr_uncompressed_coord (int): The nearest uncompressed coordinate present in the map
        uncompressed_seq (str): The uncompressed sequence being lifted
        map (dict): A dictionary mapping uncompresed indices to their corresponding compressed indices
        rle (list): A run length encoding of the uncompressed sequence
        is_end (bool): True if this corresponds to the end index of the sequence being lifted
    
    Returns:
        (int): The corresponding coordinate in compressed space
    '''
    if uncompressed_coord == curr_uncompressed_coord:
        if is_end:
            return map[curr_uncompressed_coord] + 1
        else:
            return map[curr_uncompressed_coord]
    elif curr_uncompressed_coord < uncompressed_coord:
        return lift_u_to_c_from_smaller_idx(uncompressed_coord, curr_uncompressed_coord, uncompressed_seq, map, rle, is_end)
    else:
        return lift_u_to_c_from_larger_idx(uncompressed_coord, curr_uncompressed_coord, uncompressed_seq, map, rle, is_end)

def lift_u_to_c(uncompressed_start, uncompressed_end, uncompressed_seq, id, map_file):
    '''A helper function for lifting a coordinate from uncompressed to compressed space

    Args:
        uncompressed_start (int): The 0-based, inclusive uncompressed start coordinate being lifted
        uncompressed_end (int): The 0-based, exclusive uncompressed end coordinate being lifted
        uncompressed_seq (str): The uncompressed sequence being lifted
        id (int): The fasta id of the uncompressed sequence being lifted
        map_file (str): The path to the file containing the map needed for liftover
    
    Returns:
        compressed_start, compressed_end (int, int): The corresponding coordinates in compressed space, in the form [compressed_start, compressed_end)
    '''
    map = open_map(map_file)
    rle = map[id]['rle']
    compression_map = map[id]['compression_map']
    d = {int(k):int(v) for k,v in compression_map.items()}
    keys = list(d.keys())
    keys.sort()
    curr_start = get_closest_value(keys, uncompressed_start)
    curr_end = get_closest_value(keys, uncompressed_end-1)
    compressed_start = lift_u_to_c_helper(uncompressed_start, curr_start, uncompressed_seq, d, rle, False)
    compressed_end = lift_u_to_c_helper(uncompressed_end-1, curr_end, uncompressed_seq, d, rle, True)
    return compressed_start, compressed_end

def lift_seqs(uncompressed_ref_fasta_file, seqs_bed_file, map_file, lift_compressed_to_uncompressed, out_bed_file):
    # print(uncompressed_ref_fasta_file)
    ref_fasta_seqs = SeqIO.to_dict(SeqIO.parse(open(uncompressed_ref_fasta_file),'fasta'))
    in_bed = open(seqs_bed_file, 'r')
    out_bed = open(out_bed_file, 'w')
    for line in in_bed:
        # print(line)
        line = line.strip()
        line_split = line.split('\t')
        id = line_split[0]
        fasta = ref_fasta_seqs[id]
        if lift_compressed_to_uncompressed:
            # print(int(line_split[1]), int(line_split[2]))
            start, end = lift_c_to_u(int(line_split[1]), int(line_split[2]), id, map_file)
        else:
            start, end = lift_u_to_c(int(line_split[1]), int(line_split[2]), str(fasta.seq), id, map_file)
        out_bed.write(id + '\t' + str(start) + '\t' + str(end) + '\n')
    out_bed.close()

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-cu', '--compressed_to_uncompressed', help='Take compressed coordinates and lift them to compressed space', action='store_true')
group.add_argument('-uc', '--uncompressed_to_compressed', help='Take uncompressed coordinates and lift them to uncompressed space', action='store_true')
parser.add_argument('uncompressed_ref_fasta_file', type=str, help='Path to the file containing the uncompressed reference')
parser.add_argument('map_file', type=str, help='Path to the file containing the map between compressed and uncompressed space')
parser.add_argument('in_bed_file', type=str, help='Path to the bed file containing the coordinates to be lifted over')
parser.add_argument('out_bed_file', type=str, help='Path to the file where the program generated bed file containing the new lifted coordinates will be written')
args = parser.parse_args()

lift_seqs(uncompressed_ref_fasta_file=args.uncompressed_ref_fasta_file, seqs_bed_file=args.in_bed_file, map_file=args.map_file, lift_compressed_to_uncompressed=args.compressed_to_uncompressed, out_bed_file=args.out_bed_file)