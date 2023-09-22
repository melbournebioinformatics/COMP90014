

from typing import Optional, Any
from Bio import Align
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.transforms import Affine2D
from matplotlib.markers import MarkerStyle
import numpy as np


def segment_plot(segments: list[Any], title: Optional[str]=None) -> None:
    colors = cm.rainbow(np.linspace(0, 1, len(segments)))
    t = Affine2D().rotate_deg(45)
    marker = MarkerStyle('|', 'left', t)
    for pair, c in zip(segments, colors):
        plt.plot(pair[0], pair[1], color=c, marker=marker)
    if title:
        plt.title(title)
    plt.show()
    plt.savefig('segment_plot.png')


def get_alignment_score(read1: str, read2: str, should_print: bool=False) -> int:
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0 
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -1
    alignments = aligner.align(read1, read2)
    aln = alignments[0]
    if should_print:
        print(aln)
    score = aln.score
    return score
