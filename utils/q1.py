

from typing import Any
from Bio import Align
import hashlib 


def read_fastq(fastq_path: str) -> list[str]:
    out = []
    with open(fastq_path, 'r') as fp:               
        line = fp.readline()                        
        while line:
            sequence = fp.readline().strip()        
            fp.readline()                           
            fp.readline()                           
            out.append(sequence)
            line = fp.readline()  
    return out   

def calc_overlap_sg(read1: str, read2: str, should_print: bool=False) -> Any:
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0 
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -1
    aligner.target_end_open_gap_score = 0
    aligner.target_end_extend_gap_score = 0
    aligner.query_end_open_gap_score = 0
    aligner.query_end_extend_gap_score = 0
    alignments = aligner.align(read1, read2)
    aln = alignments[0]
    if should_print:
        print(aln)
    score = aln.score
    if score < 5:
        return None
    read1_span = (aln.aligned[0][0][0], aln.aligned[0][-1][-1])
    read2_span = (aln.aligned[1][0][0], aln.aligned[1][-1][-1])
    return (score, read1_span, read2_span)

def hashfunc(text: str) -> int:
    return int(hashlib.md5(text.encode('utf-8')).hexdigest(), 16)