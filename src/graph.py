from utils import canonical_form, reverse_complement
from typing import List, Tuple, Dict, Set
from collections import defaultdict

class Edge:
    dest: int
    multiplicity: int

    def __init__(self, dest):
        self.dest = dest
        self.multiplicity = 1

    def _increase_multiplicity(self):
        self.multiplicity += 1
    
    def __repr__(self):
        return f"-> ({self.dest}, {self.multiplicity})"

class Node:
    seq: str
    id: int
    in_edges: List[Edge]
    out_edges: List[Edge]

    def __init__(self, seq, id):
        self.seq = seq
        self.id = id
        self.in_edges = List()
        self.out_edges = List()

    def _has_out_edge(self, dest):
        if not len(self.out_edges): return None
        
        for edge in self.out_edges:
            if edge.dest == dest: return edge
        return None
    
    def _has_in_edge(self, dest):
        if not len(self.in_edges): return None
        
        for edge in self.in_edges:
            if edge.dest == dest: return edge
        return None
    
    def __repr__(self):
        return self.seq
    
    def self_loop(self):
        return self._has_out_edge(self.id)
    
    def replace_edge(self, to_replace, new_node, incoming = True):
        if (incoming):
            for e in self.in_edges:
                if e.dest == to_replace: e.dest = new_node 
        else:
            for e in self.out_edges: 
                if e.dest == to_replace: e.dest = new_node 


class Graph:
    def __init__(self):
        self.next_id = 1
        self.nodes: dict[int, Node] = {}
        self.map_to_id: dict[str, int] = {}
        self.starts: List[int] = []

    def _insert_node(self, seq):
        node = Node(seq, self.next_id)
        rc = reverse_complement(node.seq)
        twin = Node(rc, -self.next_id)

        self.map_to_id[node.seq] = self.next_id
        self.map_to_id[twin.seq] = -self.next_id

        self.nodes[self.next_id] = node
        self.nodes[-self.next_id] = twin

        self.next_id += 1
    
    def _get_twin(self, node_id):
        return self.nodes[-node_id]

    def _add_edge(self, prev, next):
        if prev == 0:
            self.starts.append(next)
            return
        
        if abs(prev) == abs(next):
            # belong to the same block
            # TODO: better way to store self-loops?
            # only store self-loops as outgoing edges, so we can detect them easier later?
            # TODO: should they also be stored in the twin node?
            idx = abs(prev)
            node = self.nodes[idx]
            to_next = node._has_out_edge(idx)
            if to_next is not None: to_next._increase_multiplicity()
            else: node.out_edges.append(Edge(idx, 1))

            # to_prev = node._has_in_edge(idx)
            # if to_prev is not None: to_prev._increase_multiplicity()
            # else: node.in_edges.append(Edge(idx, 1))

            return
        
        prev_node = self.nodes[prev]
        next_node = self.nodes[next]

        prev_twin = self._get_twin(prev)
        next_twin = self._get_twin(next)

        to_next = prev_node._has_out_edge(next)
        if to_next is not None: to_next._increase_multiplicity()
        else: prev_node.out_edges.append(Edge(next, 1))

        to_prev = next_node._has_in_edge(prev)
        if to_prev is not None: to_prev._increase_multiplicity()
        else: next_node.in_edges.append(Edge(prev, 1))
        
        # connect twin nodes
        to_next = prev_twin._has_in_edge(-next)
        if to_next is not None: to_next._increase_multiplicity()
        else: prev_twin.in_edges.append(Edge(-next, 1))

        to_prev = next_twin._has_out_edge(-prev)
        if to_prev is not None: to_prev._increase_multiplicity()
        else: next_twin.out_edges.append(Edge(-prev, 1))

    def get_node_count(self):
        return self.next_id - 1

    def can_merge(self, node_a, node_b):
        # TODO: better way to do this?
        # if we store self-loops in both directions (i.e, as incoming and outgoing edges)
        # return ((len(node_a.out_edges) == 1 and not node_a.self_loop()) or (len(node_a.out_edges) == 2 and node_a.self_loop())) \
        # and ((len(node_b.in_edges) == 1 and not node_b.self_loop()) or (len(node_b.in_edges) == 2 and node_b.self_loop()))

        # if we store them only as outgoing edges
        return (len(node_a.out_edges) == 1 and node_a.out_edges[0].dest == node_b.id) and \
                (len(node_b.in_edges) == 1)
    
    def merge_nodes(self, node_a, node_b):
        seq = node_a.seq + node_b.seq[-1]
        node_a.seq = seq

        # node_a only had outgoing edge to node_b, so it simply inherits all node_b's outgoing edges
        node_a.out_edges = node_b.out_edges
        for edge in node_b.out_edges:
            # for each node_b -> node x edge, replace incoming edge from node b with incoming edge from node a
            dest = self.nodes[edge.dest]
            dest.replace_edge(node_b.id, node_a.id)
        
        # remove pointers to node b from node a
        for edge in node_b.in_edges:
            # node_a.in_edges.append(edge)
            e = node_a._has_in_edge(edge.dest)
            if e: e.multiplicity += edge.multiplicity
            else: node_a.in_edges.append(edge)
            dest = self.nodes[edge.dest]
            dest.replace_edge(node_b.id, node_a.id, False)

        # repeat for twin nodes
        twin_a = self._get_twin(node_a.id)
        twin_a.seq = reverse_complement(seq)
        twin_b = self._get_twin(node_b.id)

        twin_a.in_edges = twin_b.in_edges
        for edge in twin_b.in_edges:
            dest = self.nodes[edge.dest]
            dest.replace_edge(twin_b.id, twin_a.id, False)
        
        for edge in twin_b.out_edges:
            # twin_a.out_edges.append(edge)
            e = twin_a._has_out_edge(edge.dest)
            if e: e.multiplicity += edge.multiplicity
            else: twin_a.out_edges.append(edge)
            dest = self.nodes[edge.dest]
            dest.replace_edge(twin_b.id, twin_a.id)

        # delete entries for node b and its twin
        del self.nodes[node_b.id]
        del self.nodes[twin_b.id]


def create_pre_nodes(reads, kmer_table, hash_length, graph):
    for i in range(len(reads)):
        seq = reads[i]
        start = 0

        prev_node = 0
        
        # get initial kmer
        new_kmer = seq[start:hash_length]

        # initialize sequence of uninterrupted kmers
        consecutive_seq = seq[start:hash_length - 1]

        for end in range(hash_length - 1, len(seq)):
            # if not in initial kmer, slide kmer window
            if end >= hash_length: new_kmer = new_kmer[1:] + seq[end]
            
            can_kmer, dir = canonical_form(new_kmer)
            first_occurrence = kmer_table[can_kmer][0]

            # if newly added kmer to the window overlaps with other reads
            if len(kmer_table[can_kmer]) > 1:
                # create a new node for the previously uninterrupted sequence of kmers
                if (end - 1) - start + 1 >= hash_length:
                    curr_node = graph.next_id
                    graph._insert_node(consecutive_seq)

                    graph._add_edge(prev_node, curr_node)
                    prev_node = curr_node
                
                # create a node for the overlapping kmer
                if (i, dir, end - hash_length + 1) == first_occurrence:
                    curr_node = graph.next_id
                    graph._insert_node(new_kmer)
                    graph._add_edge(prev_node, curr_node)
                    prev_node = curr_node
                
                else:
                    curr_node = graph.map_to_id[new_kmer]
                    graph._add_edge(prev_node, curr_node)
                    prev_node = curr_node
                
                # slide window
                start = end - hash_length + 2
                consecutive_seq = consecutive_seq[1:]
            
            consecutive_seq += seq[end]
            
            # if we have reached the end of the read and there is no overlap create a node
            # for this rightmost run of uninterrupted kmers
            if end == len(seq) - 1 and len(kmer_table[can_kmer]) == 1:
                curr_node = graph.next_id
                graph._insert_node(consecutive_seq)
                graph._add_edge(prev_node, curr_node)
                prev_node = curr_node

def concatenate_nodes(graph):
    modified = True

    # keep simplifying until graph stops changing
    while modified:
        modified = False
        
        for node in graph.nodes.values():
            if node.id < 0: continue
            
            # node only has one outgoing edge
            if len(node.out_edges) == 1:
                other_end = node.out_edges[0].dest
                end_node = graph.nodes[other_end]

                if (graph.can_merge(node, end_node)):
                    graph.merge_nodes(node, end_node)
                    modified = True
            
            if modified:
                break
                

# def tour_bus(graph):
#     """
#     Identifies bubbles in the graph, i.e., paths starting from the same node, and ending in the same node.
#     These are identified through BFS.
#     """

#     # BFS

#     # when a node has already been visited -> backtrack until first common ancestor

#     # if paths are very similar -> keep consensus

# def clip_tips(graph):
#     """
#     Removes tips from the graph that are less than 2 * k in length.
#     """
#     return