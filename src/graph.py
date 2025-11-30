from utils import canonical_form, reverse_complement
from typing import List, Tuple, Dict
from collections import defaultdict

# TODO: convert this into a graph class

class Node:
    descriptor: str
    in_edges: List[Node]
    out_edges: List[Node]
    idx: int

    def __init__(self, descriptor):
        self.descriptor: str = descriptor

        self.out_edges: List[Node] = []   # adjacency list
        self.in_edges:  List[Node] = []
    
    def __repr__(self):
        return self.descriptor
 
    def _has_single_out_edge(self):
        return len(self.out_edges) == 1
    
    def _has_single_in_edge(self):
        return len(self.in_edges) == 1
    
    def _set_index(self, idx): self.idx = idx

    def __eq__(self, other):
        return self.descriptor == other.descriptor

    def __hash__(self):
        return hash(self.descriptor)
    
    def __ne__(self, other):
        return self.descriptor != other.descriptor
    
    def replace_edge(self, to_replace, new_node, incoming = True):
        if (incoming):
            new_in_edges = []
            for e in self.in_edges: 
                if e != to_replace: new_in_edges.append(e)
            new_in_edges.append(new_node)
            self.in_edges = new_in_edges
        else:
            new_out_edges = []
            for e in self.out_edges: 
                if e != to_replace: new_out_edges.append(e)
            new_out_edges.append(new_node)
            self.out_edges = new_out_edges
    
class Graph:
    def __init__(self):
        self.starts: List[Node] =  []
        self.map_to_nodes : Dict[str, int] = {}
        self.nodes : List[Tuple[Node, Node]] = [(Node(""), Node(""))]

    def _add_node(self, node: Node):
        rc = reverse_complement(node.descriptor)
        twin = Node(rc)

        self.map_to_nodes[node.descriptor] = len(self.nodes)
        self.map_to_nodes[rc] = -len(self.nodes)

        self.nodes.append((node, twin))
        return
    
    def _get_node_by_index(self, idx):
        if idx > 0: return self.nodes[idx][0]
        return self.nodes[-idx][1]

    def _get_twin(self, node: Node):
        return self._get_node_by_index(
            -self.map_to_nodes[node.descriptor]
        )
    
    def _add_edge(self, node_a, node_b):
        if not node_a:
            self.starts.append(node_b)
            return
        
        # nodes belong to the same block
        if abs(self.map_to_nodes[node_a.descriptor]) == \
            abs(self.map_to_nodes[node_b.descriptor]): return
        
        a_twin, b_twin = self._get_twin(node_a), self._get_twin(node_b)

        node_a.out_edges.append(node_b)
        node_b.in_edges.append(node_a)
            
        a_twin.in_edges.append(b_twin)
        b_twin.out_edges.append(a_twin)
    
    
    def can_merge(self, node_a, node_b):
        return len(node_a[0].out_edges) == 1 and len(node_b[0].in_edges) == 1
    
    def merge_nodes(self, start, end):
        seq = self.nodes[start][0].descriptor
        for i in range(start + 1, end + 1): seq += self.nodes[i][0].descriptor[-1]
        first_node = self.nodes[start]

        first_node[0].descriptor = seq
        first_node[1].descriptor = reverse_complement(seq)

        last_node = self.nodes[end]

        first_node[0].out_edges = last_node[0].out_edges

        for other_end in last_node[0].out_edges:
            # other end in edges should point to this "new" node and remove old pointer
            other_end.replace_edge(last_node[0], first_node[0])

        first_node[1].in_edges = last_node[1].in_edges
        
        for other_end in last_node[1].in_edges:
            # other end in incoming edges should point to this "new" node and remove old pointer
            other_end.replace_edge(last_node[1], first_node[1], False)

        self.nodes = self.nodes[:start + 1] + self.nodes[end + 1:]
    
    def get_node_by_seq(self, seq):
        return self._get_node_by_index(self.map_to_nodes[seq])

def create_pre_nodes(reads, kmer_table, hash_length, graph):
    for i in range(len(reads)):
        seq = reads[i]
        start = 0

        prev_node = None
        
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
                    node = Node(consecutive_seq)
                    graph._add_node(node)

                    # twin = graph._get_twin(node)

                    graph._add_edge(prev_node, node)
                    prev_node = node
                    # pre_nodes.append(node)
                
                # create a node for the overlapping kmer
                if (i, dir, end - hash_length + 1) == first_occurrence:
                    node = Node(new_kmer)
                    # pre_nodes.append(node)
                    graph._add_node(node)
                    # twin = graph._get_twin(node)

                    graph._add_edge(prev_node, node)
                    prev_node = node
                
                else:
                    curr_node = graph.get_node_by_seq(new_kmer)
                    # twin = graph._get_twin(curr_node)
                    graph._add_edge(prev_node, curr_node)
                    prev_node = curr_node
                
                # slide window
                start = end - hash_length + 2
                consecutive_seq = consecutive_seq[1:]
            
            consecutive_seq += seq[end]
            
            # if we have reached the end of the read and there is no overlap create a node
            # for this rightmost run of uninterrupted kmers
            if end == len(seq) - 1 and len(kmer_table[can_kmer]) == 1:
                node = Node(consecutive_seq)
                # pre_nodes.append(node)
                graph._add_node(node)
                twin = graph._get_twin(node)

                graph._add_edge(prev_node, twin)
                prev_node = node

def concatenate_nodes(graph):
    # sliding window approach, get runs of consecutive nodes that can be merged
    for start in range(len(graph.nodes)-1):
        end = start
        while end < len(graph.nodes) - 1:
            if graph.can_merge(graph.nodes[end], graph.nodes[end + 1]):
                end += 1
            else:
                break
        
        if start != end:
            graph.merge_nodes(start, end)

        if end == len(graph.nodes) - 1: break