# Bubble removal using Tourbus algorithm

# current state: figuring out what to do about twin nodes (reverse direction)
# should they still be separated or can they count as one node now?

import graph
from graph import Node

import heapq

def tourbus(graph, source_id):
    # dist: Dict[Tuple[Node, Node], int]
    # parent: Dict[Tuple[Node, Node], Tuple[Node, Node]]
    # discovered: Dict[Tuple[Node, Node], bool]
    # pq: List[Tuple[float, int, Tuple[Node, Node]]]
    dist: Dict[int, int]
    parent: Dict[int, int]
    discovered: Dict[int, bool]
    pq: List[Tuple[float, int]]
    
    dist = {node_id: float('inf') for node_id in graph.nodes}
    parent = {node_id: None for node_id in graph.nodes}
    discovered = {node_id: False for node_id in graph.nodes}

    dist[source_id] = 0.0
    discovered[source_id] = True
    pq = [(0.0, source_id)]

    while pq:
        cur_dist, u_id = heapq.heappop(pq)
        node_obj = graph.nodes[u_id]

        if cur_dist > dist[u_id]: continue
        # # if graph.nodes[node_obj].deleted: continue

        # explore all outgoing neighbors of current vertex
        for edge in node_obj.in_edges:
            dest = graph.nodes[edge.dest]

            # edge cost is the length of s(B) divided by the multiplicity of the arc leading from A to B
            edge_weight = len(dest.seq) / edge.multiplicity
            new_dist = cur_dist + edge_weight

            # if not discovered[neighbor]:
            #     discovered[neighbor] = True
            #     dist[neighbor] = new_dist
            #     parent[neighbor] = u
                # push onto heap (remember the counter tiebreaker)
            # else:
                # neighbor has already been discovered, so backtrack!

                # decide which path to keep


    return dist