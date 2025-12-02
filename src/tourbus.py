# Bubble removal using Tourbus algorithm

# current state: figuring out what to do about twin nodes (reverse direction)

import graph
from graph import Node

import heapq

# i believe this is how they calculate edge weights
def calculate_edge_weight(dest_node, multiplicity):
    '''edge cost is the length of s(B) divided by the multiplicity of the arc leading from A to B'''
    return len(dest_node.seq) / multiplicity

# probably could be made more efficient
def max_overlap(a: str, b: str) -> int:
    """
    Returns length of longest suffix of 'a' that matches a prefix of 'b'
    """
    max_overlap = min(len(a), len(b))
    for l in range(max_overlap, 0, -1):
        if a[-l:] == b[:l]:
            return l
    return 0

def concatenate_path(path: List[str]):
    result: str = path[0]

    for seq in path[1:]:
        overlap_len = max_overlap(result, seq)
        result += seq[overlap_len:]
    
    return result

def backtrack_lca(parent: Dict[int, Optional[int]], 
                  id_node_a: int, id_node_b: int) -> Tuple[List[int], List[int]]:
    ''' 
    Backtracks and finds the lowest common ancestor of node_a and node_b
    Returns both paths from lca onwards 
    '''

    ancestor: int = 0 # there is no node 0

    # collect all ancestors of A
    ancestors_a = set()
    cur = id_node_a
    while cur is not None:
        ancestors_a.add(cur)
        try:
            cur = parent[cur]
        except ValueError:
            # no more parents, cur = None
            break

    # walk B upward until LCA is found
    cur = id_node_b
    while cur is not None:
        if cur in ancestors_a:
            ancestor = cur
            break
        try:
            cur = parent[cur]
        except ValueError:
            break

    if ancestor == 0: # not a bubble
        # print("no ancestor")
        return ([], [])

    # build path A: from ancestor to node_a
    path_a: List[int] = []
    cur = id_node_a
    while cur != ancestor:
        path_a.append(cur)
        cur = parent[cur]
    path_a.append(ancestor)
    path_a.reverse()

    # build path B: from ancestor to node_b
    path_b: List[int] = []
    cur = id_node_b
    while cur != ancestor:
        path_b.append(cur)
        cur = parent[cur]
    path_b.append(ancestor)
    path_b.reverse()

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
                path_a, path_b = backtrack_lca(parent, cur_id, dest_id)
                
                if path_a == []: # not a bubble
                    continue
                
                # bubble => concatenate sequences into strings
                path_a_seqs = [graph.nodes[id_a].seq for id_a in path_a]
                path_b_seqs = [graph.nodes[id_b].seq for id_b in path_b]
                path_a_str = concatenate_path(path_a_seqs)
                path_b_str = concatenate_path(path_b_seqs)
               
                print(path_a_str, path_b_str)

                # decide which path to keep
                # scoring using global alignment


    return dist