# Bubble removal using Tourbus algorithm

# current state: figuring out what to do about twin nodes (reverse direction)
# should they still be separated or can they count as one node now?

import graph
from graph import Node

import heapq

def backtrack_path(parent: Dict[int, Optional[int]], node_id: int) -> List[int]:
    path = []
    cur = node_id
    while cur:
        path.append(cur)
        cur = parent[cur]
    path.reverse
    return path

def tourbus(graph):
    dist: Dict[int, int]
    parent: Dict[int, Optional[int]]
    discovered: Dict[int, bool]
    pq: List[Tuple[float, int]]
    
    dist = {node_id: float('inf') for node_id in graph.nodes}
    parent = {node_id: None for node_id in graph.nodes}
    discovered = {node_id: False for node_id in graph.nodes}

    for s in graph.starts:
        dist[s] = 0.0
        discovered[s] = True
        heapq.heappush(pq, (0.0, s))

    while pq:
        print(pq)
        cur_dist, cur_id = heapq.heappop(pq)
        node_obj = graph.nodes[cur_id]
        # print(node_obj.in_edges, node_obj.out_edges)

        if cur_dist > dist[cur_id]: continue
        # # if graph.nodes[node_obj].deleted: continue

        # explore all outgoing neighbors of current vertex
        for edge in node_obj.out_edges:
            dest_id = edge.dest
            dest_node = graph.nodes[dest_id]

            # edge cost is the length of s(B) divided by the multiplicity of the arc leading from A to B
            edge_weight = len(dest_node.seq) / edge.multiplicity
            new_dist = cur_dist + edge_weight

            if not discovered[dest_id]:
                discovered[dest_id] = True
                dist[dest_id] = new_dist
                parent[dest_id] = cur_id
                heapq.heappush(pq, (new_dist, dest_id))
            else:
                # neighbor has already been discovered => potential bubble
                # backtrack time!
                print(f"found neighbor ", dest_node)
                # path_cur = backtrack_path(parent, cur_id)
                # path_dest = backtrack_path(parent, dest_id)

                # lca = lca(path_cur, path_dest)


                # decide which path to keep
                # scoring using global alignment


    return dist