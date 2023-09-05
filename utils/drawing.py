

import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import graphviz_layout


def draw_suffix_trie(graph, title=None):
    fig = plt.figure(1, figsize=(7, 10), dpi=60)
    if title is not None:
        plt.title(title, fontsize=20)
    pos = graphviz_layout(graph, prog="dot")
    nx.draw(
        graph, pos, edge_color='black', width=1, linewidths=1,
        node_size=100, node_color='pink', alpha=0.9
    )
    nx.draw_networkx_edge_labels(
        graph, pos, font_color='red', font_size=20, 
        edge_labels={e: graph.edges[e]['label'] for e in graph.edges}
    )
    plt.tight_layout()
    plt.show()

def draw_travel_undirected(graph, pos=None):
    fig = plt.figure(1, figsize=(30, 30), dpi=60)
    if not pos:
        pos = nx.spring_layout(graph, seed=9, k=0.8, iterations=13)  
    nx.draw_networkx_nodes(graph, pos, node_color='white', node_size=25000, edgecolors='black', linewidths=2)
    nx.draw_networkx_edges(graph, pos, width=2)
    nx.draw_networkx_labels(graph, pos, font_size=25, font_family="sans-serif")
    nx.draw_networkx_edge_labels(
        graph, pos, font_color='red', font_size=40, 
        edge_labels={e: graph.edges[e]['label'] for e in graph.edges}
    )
    plt.tight_layout()
    plt.show()
    return pos

def draw_travel_directed(graph, pos):
    fig = plt.figure(1, figsize=(30, 30), dpi=60)
    nx.draw_networkx_nodes(graph, pos, node_color='white', node_size=25000, edgecolors='black', linewidths=2)
    nx.draw_networkx_edges(graph, pos, width=3, arrows=True, arrowstyle='-|>', arrowsize=50, min_target_margin=80)
    nx.draw_networkx_labels(graph, pos, font_size=25, font_family="sans-serif")
    nx.draw_networkx_edge_labels(
        graph, pos, font_color='red', font_size=40, 
        edge_labels={e: graph.edges[e]['label'] for e in graph.edges}
    )
    plt.tight_layout()
    plt.show()

def draw_evo_tree(graph):
    fig = plt.figure(1, figsize=(7, 7), dpi=60)
    pos = graphviz_layout(graph, prog="dot")
    nx.draw(
        graph, pos, edge_color='black', width=1, linewidths=1,
        node_size=1500, node_color='pink', alpha=0.9
    )
    nx.draw_networkx_labels(graph, pos, font_size=16, font_family="sans-serif")
    nx.draw_networkx_edge_labels(
        graph, pos, font_color='red', font_size=16, 
        edge_labels={e: f"{graph.edges[e]['label']:0.2f}" for e in graph.edges}
    )
    plt.show()

