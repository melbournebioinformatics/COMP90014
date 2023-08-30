

from typing import Optional
import networkx as nx
from collections import Counter

def kmer_dist(seq1: str, seq2: str, k: int) -> Optional[int]:
    if len(seq1) < k or len(seq2) < k:
        return None
    
    kmers1 = []
    for i in range(len(seq1) - k + 1):
        kmers1.append(seq1[i:i+k])
    
    kmers2 = []
    for i in range(len(seq2) - k + 1):
        kmers2.append(seq2[i:i+k])
    
    all_kmers = set(kmers1 + kmers2)
    counts1 = Counter(kmers1)
    counts2 = Counter(kmers2)

    dist = 0
    for kmer in all_kmers:
        dist += abs(counts1[kmer] - counts2[kmer])
    return dist 

def gen_cluster_name(c1: str, c2: str) -> str:
    inode_charlist = [*c1] + [*c2]
    inode_charlist.sort()
    cluster_name = ''.join(inode_charlist)
    return cluster_name

def split_cluster_name(name: str) -> list[str]:
    cname = [*name]
    cname.sort()
    return cname

def get_branch_length(node: str, tree: nx.DiGraph) -> float:
    length = 0
    cnode = node
    while len(list(tree.neighbors(cnode))):
        neighbors = list(tree.neighbors(cnode))
        max_label = max([tree.edges[(cnode, n)]['label'] for n in neighbors])
        edge = [e for e in tree.edges if tree.edges[e]['label'] == max_label][0]
        length += max_label # type: ignore
        nnode = edge[1]
        cnode = nnode
    return length