# Bubble removal using Tourbus algorithm

# current state: figuring out what to do about twin nodes (reverse direction)
# should they still be separated or can they count as one node now?

import graph
from graph import Node

import heapq
import itertools
from collections import Counter

def tourbus(graph, source):
    dist: Dict[Tuple[Node, Node], int]
    parent: Dict[Tuple[Node, Node], Tuple[Node, Node]]
    discovered: Dict[Tuple[Node, Node], bool]
    pq: List[Tuple[float, int, Tuple[Node, Node]]]
    multiplicities: Dict[Tuple[Node, Node], Tuple[Tuple[Node, Node], int]]
    
    dist = {node: float('inf') for node in graph.nodes}
    parent = {node: None for node in graph.nodes}
    discovered = {node: False for node in graph.nodes}

    dist[source] = 0.0
    discovered[source] = True

    # for use as tiebreaker in priority queue
    tie_counter = itertools.count()
    pq = [(0.0, next(tie_counter), source)]

    while pq:
        cur_dist, _, u = heapq.heappop(pq)

        if cur_dist > dist[u]: continue
        # if graph.nodes[u].deleted: continue

        # precompute dictionary of multiplicities
        multiplicities = dict(Counter(u[0].out_edges))

        # explore all neighbors (outgoing) of current vertex
        for neighbor in multiplicities:
            # edge cost is the length of s(B) divided by the multiplicity of the arc leading from A to B
            edge_weight = len(neighbor.descriptor) / multiplicities[neighbor]
            new_dist = cur_dist + edge_weight

            # print("hello node", u)
            # print("hey neighbor", neighbor)
            # if not discovered[neighbor]:
            #     discovered[neighbor] = True
            #     dist[neighbor] = new_dist
            #     parent[neighbor] = u
                # push onto heap (remember the counter tiebreaker)
            # else:
                # neighbor has already been discovered, so backtrack!

                # decide which path to keep


    return dist