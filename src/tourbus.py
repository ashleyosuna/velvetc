# Bubble removal using Tourbus algorithm

import graph

import heapq

def tourbus(graph, source):
    dist = {node: float('inf') for node in graph.nodes}
    parent = {node: None for node in graph.nodes}
    discovered = {node: False for node in graph.nodes}

    dist[source] = 0
    discovered[source] = True

    pq = [(0, source)]
    return dist