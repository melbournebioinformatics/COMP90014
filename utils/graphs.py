

def check_graph_nodes(edges, graph):
    for edge in edges:
        if edge[0] not in graph.nodes:
            return False
        if edge[1] not in graph.nodes:
            return False
    return True

def check_graph_edges(edges, graph):
    for e in edges:
        try:
            graph.edges[e[0], e[1]]
        except KeyError:
            return False
    return True

def check_graph_labels(edges, graph):
    all_edges = list(graph.edges(data=True))
    for e in edges:
        for edge in all_edges:
            if e[0] == edge[0] and e[1] == edge[1]:
                if e[2] != edge[2]['label']:
                    return False
    return True