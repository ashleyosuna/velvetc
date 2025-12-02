# Bubble removal using Tourbus algorithm

# current state: figuring out what to do about twin nodes (reverse direction)
# should they still be separated or can they count as one node now?

import graph
from graph import Node

import heapq

# i think this is how they calculate edge weights?
def calculate_edge_weight(dest_node, multiplicity):
    '''edge cost is the length of s(B) divided by the multiplicity of the arc leading from A to B'''
    return len(dest_node.seq) / multiplicity

def backtrack_lca(parent: Dict[int, Optional[int]], id_node_a: int, id_node_b: int) -> Tuple[List[int], List[int]]:
    ''' 
    Backtracks and finds the lowest common ancestor of node_a and node_b
    Returns both paths from lca onwards 
    '''

    path_a: List[int] = []
    path_b: List[int] = []

    cur_a: int = id_node_a
    cur_b: int = id_node_b

    ancestor: int
    # just return the path with ancestor for now, then figure out splicing later
    # issue: codon doesn't keep order of list?

    # this seems highly inefficient, searching through list every time
    while cur_a and cur_b:
        if cur_a == cur_b:
            path_a.append(cur_a)
            path_b.append(cur_b)
            ancestor = cur_a
            break
        elif cur_a in path_b:
            path_a.append(cur_a)
            ancestor = cur_a
            break
        elif cur_b in path_a:
            path_b.append(cur_b)
            ancestor = cur_b
            break
        else:
            path_a.append(cur_a)
            path_b.append(cur_b)
            try:
                cur_a = parent[cur_a]
                cur_b = parent[cur_b]
            except ValueError: # this means there is no lowest common ancestor -- shouldn't happen, right?
                print("no lowest common ancestor")
                path_a.reverse, path_b.reverse
                return (path_a, path_b)

    path_a.reverse, path_b.reverse
    return (path_a, path_b)

def tourbus(graph):
    dist: Dict[int, int] = {}
    parent: Dict[int, Optional[int]] = {}
    discovered: Dict[int, bool] = {}
    pq: List[Tuple[float, int]] = []
    
    dist = {node_id: float('inf') for node_id in graph.nodes}
    parent = {node_id: None for node_id in graph.nodes}
    discovered = {node_id: False for node_id in graph.nodes}

    start_nodes = list(graph.starts)
    for s in start_nodes:
        dist[s] = 0.0
        discovered[s] = True
        heapq.heappush(pq, (0.0, s))

    while pq:
        cur_dist, cur_id = heapq.heappop(pq)
        node_obj = graph.nodes[cur_id]
        # print(node_obj.in_edges, node_obj.out_edges)

        if cur_dist > dist[cur_id]: continue
        # # if graph.nodes[node_obj].deleted: continue

        # explore all outgoing neighbors of current vertex
        for edge in node_obj.out_edges:
            dest_id = edge.dest
            dest_node = graph.nodes[dest_id]

            edge_weight = calculate_edge_weight(dest_node, edge.multiplicity)
            new_dist = cur_dist + edge_weight

            if not discovered[dest_id]:
                discovered[dest_id] = True
                dist[dest_id] = new_dist
                parent[dest_id] = cur_id
                heapq.heappush(pq, (new_dist, dest_id))
            else:
                # neighbor has already been discovered => potential bubble
                # backtrack time!
                print(f"cur id ", cur_id)
                print(f"found neighbor ", dest_node)
                path_a, path_b = backtrack_lca(parent, cur_id, dest_id)
                print(path_a, path_b)

                # decide which path to keep
                # scoring using global alignment


    return dist