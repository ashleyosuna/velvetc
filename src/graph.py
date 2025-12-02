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
    
    def __ne__(self, other):
        return self.descriptor != other.descriptor

    def __hash__(self):
        return hash(self.descriptor)

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
    
    def remove_edge(self, dest, incoming = True):
        if (incoming):
            edges = []
            for e in self.in_edges:
                if e.dest != dest: edges.append(e)
            self.in_edges = edges
        else:
            edges = []
            for e in self.out_edges:
                if e.dest != dest: edges.append(e)
            self.out_edges = edges

    def length(self, k):
        return len(self.seq) - k + 1 


class Graph:
    def __init__(self, hash_length):
        self.next_id = 1
        self.nodes: dict[int, Node] = {}
        self.map_to_id: dict[str, int] = {}
        self.starts: List[int] = []
        self.hash_length = hash_length
        self.enable_twin = True # for debugging

    def _insert_node(self, seq):
        node = Node(seq, self.next_id)
        rc = reverse_complement(node.seq)

        self.map_to_id[node.seq] = self.next_id

        self.nodes[self.next_id] = node

        if self.enable_twin:
            twin = Node(rc, -self.next_id)
            self.map_to_id[twin.seq] = -self.next_id
            self.nodes[-self.next_id] = twin

        self.next_id += 1
    
    def length(self, node):
        return node.length(self.hash_length)
    
    def _get_twin(self, node_id):
        return self.nodes[-node_id]

    def _add_edge(self, prev, next):
        if prev == 0:
            self.starts.append(next)
            return
        
        # if abs(prev) == abs(next):
            # belong to the same block
            # TODO: better way to store self-loops?
            # only store self-loops as outgoing edges, so we can detect them easier later?
            # TODO: should they also be stored in the twin node?
            # idx = abs(prev)
            # node = self.nodes[idx]
            # to_next = node._has_out_edge(idx)
            # if to_next is not None: to_next._increase_multiplicity()
            # else: node.out_edges.append(Edge(idx, 1))

            # to_prev = node._has_in_edge(idx)
            # if to_prev is not None: to_prev._increase_multiplicity()
            # else: node.in_edges.append(Edge(idx, 1))

            # return
        
        prev_node = self.nodes[prev]
        next_node = self.nodes[next]

        to_next = prev_node._has_out_edge(next)
        if to_next is not None: to_next._increase_multiplicity()
        else: prev_node.out_edges.append(Edge(next, 1))

        to_prev = next_node._has_in_edge(prev)
        if to_prev is not None: to_prev._increase_multiplicity()
        else: next_node.in_edges.append(Edge(prev, 1))

        if self.enable_twin:
            prev_twin = self._get_twin(prev)
            next_twin = self._get_twin(next)

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
        # return (len(node_a.out_edges) == 1 and node_a.out_edges[0].dest == node_b.id) and \
        #         (len(node_b.in_edges) == 1)
        if abs(node_a.id) == abs(node_b.id): return False

        return len(node_a.out_edges) == 1 and node_a.out_edges[0].dest == node_b.id and \
            len(node_b.in_edges) == 1 and node_b.in_edges[0].dest == node_a.id
    
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


        if self.enable_twin:
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

        new_starts = []
        for s in self.starts:
            if s != node_b.id: new_starts.append(s)
        
        self.starts = new_starts

        # delete entries for node b and its twin
        del self.nodes[node_b.id]
        if self.enable_twin:
            del self.nodes[twin_b.id]

    def destroy_tip(self, node_id):
        curr = self.nodes[node_id]
        twin = self.nodes[-node_id]
        prev = None

        to_delete = []

        while (len(curr.in_edges) == 1 and len(curr.out_edges) == 0 and len(twin.out_edges) == 1 and len(twin.in_edges) == 0) or \
                (len(curr.out_edges) == 1 and len(curr.in_edges) == 0 and len(twin.in_edges) == 1 and len(twin.out_edges) == 0):
            
            twin = self.nodes[-curr.id]
            del self.map_to_id[curr.seq]
            del self.map_to_id[twin.seq]

            to_delete.append(curr.id)
            to_delete.append(-curr.id)

            prev = curr
            curr = self.nodes[curr.out_edges[0].dest]
            twin = self.nodes[-curr.id]
        
        if prev:
            other_end = self.nodes[prev.out_edges[0].dest]

            other_end.remove_edge(prev.id)

            prev_twin = self.nodes[-prev.id]
            other_end = self.nodes[prev_twin.in_edges[0].dest]

            other_end.remove_edge(prev_twin.id, False)
        
        filtered_start = []
        for s in self.starts:
            if s != node_id: filtered_start.append(s)
        self.starts = filtered_start
        
        for n in to_delete: del self.nodes[n]


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
                    if new_kmer not in graph.map_to_id:
                        print("Missing kmer:", new_kmer, " (twin disabled?)")
                        continue

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

def clip_tips(graph):
    """
    Removes tips from the graph that are less than 2 * k in length.
    """    
    modified = True

    while modified:
        modified = False

        for s in graph.starts:
            node = graph.nodes[s]
            twin = graph.nodes[-s]
            length = 0
            simple = False

            while (len(node.in_edges) == 1 and len(node.out_edges) == 0 and len(twin.out_edges) == 1 and len(twin.in_edges) == 0) or \
                (len(node.out_edges) == 1 and len(node.in_edges) == 0 and len(twin.in_edges) == 1 and len(twin.out_edges) == 0):
                if length == 0: length += graph.length(node)
                else: length += graph.length(node) - (graph.hash_length - 1)

                node = graph.nodes[node.out_edges[0].dest]
                twin = graph.nodes[-node.id]

                simple = True

            if simple and length < 2 * graph.hash_length:

                graph.destroy_tip(s)

                modified = True

            if modified: break
    return