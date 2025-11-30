# Bubble removal using Tourbus algorithm

import graph
from graph import Node

import heapq
import itertools

def tourbus(graph, source):
    dist: Dict[Tuple[Node, Node], int]
    parent: Dict[Tuple[Node, Node], Tuple[Node, Node]]
    discovered: Dict[Tuple[Node, Node], bool]
    pq: List[Tuple[int, int, Tuple[Node, Node]]]
    
    dist = {node: float('inf') for node in graph.nodes}
    parent = {node: None for node in graph.nodes}
    discovered = {node: False for node in graph.nodes}

    dist[source] = 0
    discovered[source] = True

    # for use as tiebreaker during dijkstra
    tie_counter = itertools.count()
    pq = [(0, next(tie_counter), source)]

    while pq:
        cur_dist, _, u = heapq.heappop(pq)

    return dist