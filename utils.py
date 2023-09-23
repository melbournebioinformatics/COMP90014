

from Bio import Align
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.transforms import Affine2D
from matplotlib.markers import MarkerStyle

###############
### FILE IO ###
###############


def read_fastq(fastq_path):
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


#####################
### VISUALISATION ###
#####################

def draw_overlap_graph(graph):
    fig = plt.figure(1, figsize=(5, 5))
    pos = nx.circular_layout(graph)
    nx.draw_networkx_edges(graph, pos, connectionstyle='arc3', width=0.7, arrows=True, arrowstyle='-|>', arrowsize=10, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_labels(graph, pos, font_size=9, font_family="sans-serif")
    nx.draw_networkx_edge_labels(
        graph, pos, font_color='red', font_size=8, 
        edge_labels={e: graph.edges[e]['label'] for e in graph.edges}
    )
    plt.tight_layout()
    plt.show()

def draw_interval_tree(graph, title):
    from networkx.drawing.nx_agraph import graphviz_layout  # TODO REMOVE
    fig = plt.figure(1, figsize=(20, 30), dpi=60)
    if title is not None:
        plt.title(title, fontsize=20)
    pos = graphviz_layout(graph, prog="dot")
    nx.draw_networkx_nodes(graph, pos, node_color='white', node_size=5000, edgecolors='black', linewidths=2)
    nx.draw_networkx_edges(graph, pos, width=2)
    nx.draw_networkx_labels(graph, pos, font_size=25, font_family="sans-serif")
    nx.draw_networkx_edge_labels(
        graph, pos, font_color='red', font_size=40, 
        edge_labels={e: graph.edges[e]['label'] for e in graph.edges}
    )
    plt.tight_layout()
    plt.show()
    plt.savefig('IST.png') # TODO REMOVE

def segment_plot(segments, title=None):
    colors = cm.rainbow(np.linspace(0, 1, len(segments)))
    t = Affine2D().rotate_deg(45)
    marker = MarkerStyle('|', 'left', t)
    for pair, c in zip(segments, colors):
        plt.plot(pair[0], pair[1], color=c, marker=marker)
    if title:
        plt.title(title)
    plt.show()


#################
### ALIGNMENT ###
#################

def calc_overlap_sg(read1, read2, should_print=False):
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
    if score == 0:
        return None
    read1_span = (aln.aligned[0][0][0], aln.aligned[0][-1][-1])
    read2_span = (aln.aligned[1][0][0], aln.aligned[1][-1][-1])
    return (score, read1_span, read2_span)

def get_alignment_score(read1, read2, should_print=False):
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