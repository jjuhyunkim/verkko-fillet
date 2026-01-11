import pandas as pd
from dataclasses import dataclass
from typing import List, Tuple, Optional
import re 
import os 
import subprocess
import copy
from ._findNode_from_region import readGraph, getNodeSpace_from_onePath
from ..tools._insert_gap import insertGap

@dataclass
class NodeCoordinates:
    start: int
    end: int

@dataclass
class ClippedPAF:
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alnlen: int
    mapq: int

# 


    



# ===== Expected core columns in your DataFrame =====
PAF_COLS = [
    "qname","qlen","qstart","qend","strand",
    "tname","tlen","tstart","tend","nmatch","alnlen","mapq"
]

# ===== CIGAR helpers (supports = X M I D) =====

def parse_len_op_cigar(cigar: str, assume_final_equals: bool = True) -> List[Tuple[str, int]]:
    """Parse len-op CIGAR like '100=2D50X3I20=' -> [('=',100),('D',2),('X',50),('I',3),('=',20)]."""
    out, num = [], []
    for ch in cigar:
        if ch.isdigit():
            num.append(ch)
        else:
            if not num:
                raise ValueError(f"Bad CIGAR near '{ch}' in '{cigar[:80]}...'")
            out.append((ch, int("".join(num))))
            num = []
    if num:
        if assume_final_equals:
            out.append(("=", int("".join(num))))
        else:
            raise ValueError("Trailing number without op in CIGAR.")
    return out

def compress_cigar(ops: List[Tuple[str,int]]) -> str:
    """Coalesce adjacent same-op runs and drop zero-length ops."""
    out: List[Tuple[str,int]] = []
    for op, n in ops:
        if n <= 0:
            continue
        if out and out[-1][0] == op:
            out[-1] = (op, out[-1][1] + n)
        else:
            out.append((op, n))
    return "".join(f"{n}{op}" for op, n in out)

def cigar_counts_and_ends(cigar: str, qstart: int, tstart: int):
    """
    Return:
      qend, tend, nmatch, alnlen, XI, XD, NM
    where (PAF spec):
      nmatch = exact matches only ( '=' )
      alnlen = matches + mismatches + gaps  (= + X + M + I + D)
      NM     = mismatches + gaps            ( X + I + D )
    """
    ops = parse_len_op_cigar(cigar, assume_final_equals=True)
    q = qstart; t = tstart
    n_eq = n_x = n_m = n_I = n_D = 0
    for op, n in ops:
        if op in ("=", "X", "M"):
            q += n; t += n
            if op == "=": n_eq += n
            elif op == "X": n_x += n
            else: n_m += n
        elif op == "I":
            q += n; n_I += n
        elif op == "D":
            t += n; n_D += n
        else:
            # ignore S/H/P etc. for PAF core stats
            pass
    qend = q; tend = t
    nmatch = n_eq
    alnlen = n_eq + n_x + n_m + n_I + n_D  # **includes gaps**
    NM = n_x + n_I + n_D
    return qend, tend, nmatch, alnlen, n_I, n_D, NM

# ===== Clip one alignment row to [keep_t0, keep_t1) on target =====

@dataclass
class ClipResult:
    qstart_new: int
    qend_new: int
    tstart_new: int
    tend_new: int
    cigar_new: str
    nmatch: int
    alnlen: int
    XI: int
    XD: int
    NM: int

def clip_cigar_to_target_interval(
    cigar: str,
    qstart: int,
    tstart: int,
    keep_t0: int,
    keep_t1: int,
) -> Optional[ClipResult]:
    """
    Keep subsegments that occur while target t is in [keep_t0, keep_t1).
    Rules:
      =/X/M -> keep overlapping length (consumes both)
      D     -> keep overlapping length on target (consumes t only)
      I     -> keep if it occurs while current t is inside the interval (consumes q only)
    """
    ops = parse_len_op_cigar(cigar, assume_final_equals=True)
    q = qstart; t = tstart

    started = False
    out_ops: List[Tuple[str,int]] = []
    qstart_new: Optional[int] = None

    def begin_if_needed(qpos):
        nonlocal started, qstart_new
        if not started:
            started = True
            qstart_new = qpos

    for op, n in ops:
        if op in ("=", "X", "M"):
            seg_t0, seg_t1 = t, t + n
            inter0 = max(seg_t0, keep_t0)
            inter1 = min(seg_t1, keep_t1)
            inter_len = inter1 - inter0
            if inter_len > 0:
                offset = inter0 - seg_t0
                # advance into this op up to first kept base
                q += offset; t += offset
                begin_if_needed(q)
                # keep the overlapping paired chunk (preserve op; don't coerce 'M' to '=')
                out_ops.append((op, inter_len))
                # advance past kept part
                q += inter_len; t += inter_len
                # advance any trailing outside part
                tail = seg_t1 - inter1
                q += tail; t += tail
            else:
                q += n; t += n

        elif op == "D":
            seg_t0, seg_t1 = t, t + n
            inter0 = max(seg_t0, keep_t0)
            inter1 = min(seg_t1, keep_t1)
            inter_len = inter1 - inter0
            if inter_len > 0:
                begin_if_needed(q)
                out_ops.append(("D", inter_len))
            t += n

        elif op == "I":
            if keep_t0 <= t < keep_t1:
                begin_if_needed(q)
                out_ops.append(("I", n))
            q += n

        else:
            # ignore S/H/P/etc (no target consumption)
            pass

    if not started or not out_ops:
        return None

    cigar_new = compress_cigar(out_ops)
    tstart_new = max(tstart, keep_t0)
    qend_new, tend_new, nmatch, alnlen, XI, XD, NM = cigar_counts_and_ends(cigar_new, qstart_new, tstart_new)

    return ClipResult(
        qstart_new=qstart_new, qend_new=qend_new,
        tstart_new=tstart_new, tend_new=tend_new,
        cigar_new=cigar_new,
        nmatch=nmatch, alnlen=alnlen,
        XI=XI, XD=XD, NM=NM
    )

# ===== Utilities to extract/merge tags in a flexible DF =====

def _extract_cg_from_row(row: pd.Series) -> Optional[str]:
    """
    Try to find cg (CIGAR) in common layouts:
      1) 'tags' dict with 'cg'
      2) a 'cg' column (value may be 'cg:Z:...' or just the cigar)
      3) any extra column that contains 'cg:Z:...'
    """
    # 1) tags dict
    if "tags" in row and isinstance(row["tags"], dict) and "cg" in row["tags"]:
        return str(row["tags"]["cg"])

    # 2) 'cg' column
    if "cg" in row.index and pd.notna(row["cg"]):
        s = str(row["cg"])
        return s.split("cg:Z:")[-1] if "cg:Z:" in s else s

    # 3) scan other columns for a raw tag string
    for col in row.index:
        if col in PAF_COLS or col == "tags":
            continue
        v = row[col]
        if isinstance(v, str) and v.startswith("cg:Z:"):
            return v[5:]
    return None

def _merge_tags(row: pd.Series, cigar_new: str, NM: int, XI: int, XD: int) -> dict:
    """
    Preserve existing tags as much as possible; overwrite cg/NM/XI/XD.
    If row['tags'] exists and is a dict, start from it. Otherwise, try to
    pick up any raw tag-like strings from extra columns.
    """
    tags = {}
    if "tags" in row and isinstance(row["tags"], dict):
        tags = dict(row["tags"])
    else:
        # Slurp any 'k:TYPE:val' strings found in non-core columns
        for col in row.index:
            if col in PAF_COLS:
                continue
            v = row[col]
            if isinstance(v, str) and v.count(":") >= 2:
                k = v.split(":", 1)[0]
                if k in ("cg","NM","XI","XD"):
                    continue
                tags[k] = v  # keep raw

    # overwrite core edit tags
    tags["cg"] = cigar_new
    tags["NM"] = NM
    tags["XI"] = XI
    tags["XD"] = XD
    return tags

# ===== Main DF function (no file I/O) =====

def clip_paf_df(df: pd.DataFrame, tname: str, tstart: int, tend: int) -> pd.DataFrame:
    """
    Clip an already-loaded PAF DataFrame to target interval [tstart, tend) on chromosome tname.
    Returns a NEW DataFrame with columns PAF_COLS + 'tags' (dict). Rows without overlap or without cg are dropped.
    """
    kept = []
    for _, row in df.iterrows():
        if row["tname"] != tname:
            continue
        if row["tend"] <= tstart or row["tstart"] >= tend:
            continue

        cg = _extract_cg_from_row(row)
        if not cg:
            continue  # require CIGAR to clip

        res = clip_cigar_to_target_interval(
            cigar=cg,
            qstart=int(row["qstart"]),
            tstart=int(row["tstart"]),
            keep_t0=int(tstart),
            keep_t1=int(tend),
        )
        if res is None:
            continue

        tags_new = _merge_tags(row, res.cigar_new, res.NM, res.XI, res.XD)

        kept.append({
            "qname":  row["qname"],
            "qlen":   int(row["qlen"]),
            "qstart": int(res.qstart_new),
            "qend":   int(res.qend_new),
            "strand": row["strand"],
            "tname":  row["tname"],
            "tlen":   int(row["tlen"]),
            "tstart": int(res.tstart_new),
            "tend":   int(res.tend_new),
            "nmatch": int(res.nmatch),     # '=' only
            "alnlen": int(res.alnlen),     # matches + mismatches + gaps (I/D/M/X)
            "mapq":   int(row["mapq"]),
            "tags":   tags_new
        })

    return pd.DataFrame(kept, columns=PAF_COLS + ["tags"])

# ===== Optional: serialize clipped DF back to PAF text =====

def paf_line_from_row(row: pd.Series) -> str:
    """Serialize one clipped row (with 'tags' dict) back into a PAF line."""
    core = [str(row[c]) for c in PAF_COLS]
    tag_strings = []
    for k, v in row["tags"].items():
        # Pass through raw tags like 'tp:A:P' if present
        if isinstance(v, str) and v.count(":") >= 2 and v.split(":", 1)[0] == k:
            tag_strings.append(v)
        elif isinstance(v, int):
            tag_strings.append(f"{k}:i:{v}")
        else:
            tag_strings.append(f"{k}:Z:{v}")
    return "\t".join(core + tag_strings)


## Update gap 
#

def AddNewGap(obj, contig, start_end_node):
    gapdb = obj.gaps.copy()
    pathdb = obj.paths.copy()

    pathdb = pathdb[pathdb['name'] == contig].reset_index(drop=True)
    path = pathdb['path'].str.split(",").tolist()[0]

    # index searching to check if the order is correct
    with_orientation = []
    for node in start_end_node:
        for orientation in ["-", "+"]:
            if node + orientation not in path:
                continue
            else:
                print(f"Found {node + orientation} in path")
                break
        node = node + orientation
        with_orientation.append(node)

    idx = [path.index(node) for node in with_orientation]

    ## Get Original list between start and end
    min_idx = min(idx)
    max_idx = max(idx)
    original_node_list = path[min_idx:max_idx+1]
    # print(original_node_list)
    gap_id = "gapid_" + str(len(gapdb))
    print(f"New gap_id: {gap_id}")

    newGapdb = pd.DataFrame({
        'gapId': gap_id,
        'name': contig,
        'gaps' : [original_node_list],
    })

    # Ensure newGapdb has all columns as in gapdb, fill with "" if missing
    for col in gapdb.columns:
        if col not in newGapdb.columns:
            newGapdb[col] = ""
    newGapdb = newGapdb[gapdb.columns]  # Ensure column order matches

    gapdb = pd.concat([gapdb, newGapdb], ignore_index=True)
    obj.gaps = gapdb
    # print(f"Updated gapdb:\n{obj.gaps.tail()}")
    return obj, gap_id



def patchFromOtherAsm(graph, start_end_node, aln_edge, gap_id, gaf_sub, patching_sequence):
    # get node location on the path
    # get graph topology 
    segment, link = readGraph(graph)
    # segment.head(2)
    # link.head(2)
    node_list_in_path = [t for t in re.split(r'(?=<)|(?=>)', gaf_sub['tname'].values[0]) if t]
    node_list_in_path_plain = [t for t in re.split(r'<|>', gaf_sub['tname'].values[0]) if t]
    node_space = getNodeSpace_from_onePath(node_list_in_path_plain, segment, link)
    start_node_coor_start, start_node_coor_end = node_space.loc[node_space['node'] == start_end_node[0], ['start_coor', 'end_coor']].values[0]
    end_node_coor_start, end_node_coor_end = node_space.loc[node_space['node'] == start_end_node[1], ['start_coor', 'end_coor']].values[0]

    ## Check the size of the nodes 
    for node in start_end_node:
        node_len = segment.loc[segment['node']==node, 'len'].values[0]
        if node_len < aln_edge:
            raise ValueError(f"Node {node} length {node_len} is less than alignment block size {aln_edge}. Please choose different nodes.")
        # continue
    # Extract the elements that are ">" or "<" + startNode and endNode
    with_orientation = []
    for node in start_end_node:
        for orientation in ["<", ">"]:
            if orientation + node not in node_list_in_path:
                continue
            with_orientation.append(orientation + node)
    if len(with_orientation) != 2:
        raise ValueError(f"Error: Could not find both start and end nodes with orientation in the path: {node_list_in_path}")
    print(f"Found start and end nodes with orientation: {with_orientation}")

    ## 02 get numbers
    # get position
    tstart = gaf_sub.iloc[0]['tstart']
    tend = gaf_sub.iloc[0]['tend']
    # get list of nodes 

    clipped_start_node_paf = clip_paf_df(gaf_sub, tname=gaf_sub['tname'].values[0], tstart=max(start_node_coor_start, tstart), tend=start_node_coor_end-aln_edge)
    # clipped_start_node_paf['tname'] = with_ori_in_path[0]
    clipped_start_node_paf['tname'] = node_list_in_path[0]
    clipped_start_node_paf['tlen'] = segment.loc[segment['node']==start_end_node[0], 'len'].values[0]

    clipped_end_node_paf = clip_paf_df(gaf_sub, tname=gaf_sub['tname'].values[0], tstart=end_node_coor_start+aln_edge, tend=min(end_node_coor_end, tend))
    # clipped_end_node_paf['tname'] = with_ori_in_path[1]
    clipped_end_node_paf['tname'] = node_list_in_path[-1]
    clipped_end_node_paf['tlen'] = segment.loc[segment['node']==start_end_node[1], 'len'].values[0]
    clipped_end_node_paf['tstart'] = clipped_end_node_paf['tstart'] - end_node_coor_start
    clipped_end_node_paf['tend'] = clipped_end_node_paf['tend'] - end_node_coor_start

    qend = clipped_end_node_paf['qend'].values[0]
    
    qname = clipped_end_node_paf['qname'].values[0]
    clipped_paf = pd.concat([clipped_start_node_paf, clipped_end_node_paf], ignore_index=True)
    clipped_paf = clipped_paf[['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 'tstart', 'tend', 
                    'nmatch', 'alnlen', 'mapq']]
    clipped_paf['qlen']= qend
    # clipped_paf['qname']= f"{qname}_{gap_id}"
    clipped_paf.to_csv(f"missing_edge/patch.nogap.{gap_id}.filtered.gaf", sep="\t", header=False, index=False)
    print(f"Written clipped paf to missing_edge/patch.nogap.{gap_id}.filtered.gaf")

    ## Final path 
    print(f"Making new node: {gap_id} for patching")
    insertGap(
        gap_id,
        gaf_sub,
        outputDir="missing_edge",
        alignGAF=f"missing_edge/patch.nogap.{gap_id}.filtered.gaf",
        graph="assembly.homopolymer-compressed.gfa"
    )

    os.rename(f"missing_edge/{gap_id}.missing_edge.ont_list.txt", f"missing_edge/{gap_id}.missing_edge.patching_list.txt")

    qstart = min(clipped_paf['qstart'])
    qend = max(clipped_paf['qend'])

    cmd = f"samtools faidx {patching_sequence} {qname}:{qstart+1}-{qend} | sed -e '1d' | sed -e '1i >{qname}_{gap_id} ' > missing_edge/{gap_id}.missing_edge.patching.fasta"
    subprocess.run(cmd, shell=True, check=True)

    for roi in [qname] :
        cmd = f"sed -i 's/{roi}/{roi}_{gap_id}/g' missing_edge/patch.nogap.{gap_id}.gaf"
        subprocess.run(cmd, shell=True, check=True)     

        cmd = f"sed -i 's/{roi}/{roi}_{gap_id}/g' missing_edge/patch.{gap_id}.gaf"
        subprocess.run(cmd, shell=True, check=True)

        cmd = f"sed -i 's/{roi}/{roi}_{gap_id}/g' missing_edge/{gap_id}.missing_edge.patching.gfa"
        subprocess.run(cmd, shell=True, check=True)

    return clipped_paf