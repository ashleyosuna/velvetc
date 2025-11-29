from utils import canonical_form, reverse_complement
from typing import List, Tuple, Dict
from collections import defaultdict

# TODO: convert this into a graph class

class Node:
    descriptor: str
    in_edges: List[Node]
    out_edges: List[Node]

    def __init__(self, descriptor):
        self.descriptor: str = descriptor

        self.out_edges: List[Node] = []   # adjacency list
        self.in_edges:  List[Node] = []
    
    def __repr__(self):
        return self.descriptor
 
    def _has_single_out_edge(self):
        # TODO: make this more efficient
        # incoming_edges = 0

        # for (descriptor, _) in self._map_to_nodes.items():
        #     if descriptor == node.descriptor: continue

        #     for (_, end) in self.edges[descriptor]:
        #         if end.descriptor == node.descriptor: incoming_edges += 1
        #         if incoming_edges > 1:
        #             return False
        # return True
        return len(self.out_edges) == 1
    
    def _has_single_in_edge(self):
        return len(self.in_edges) == 1
    
class Graph:
    def __init__(self):
        self.source = Node('-1')
        # self.blocks : List[Tuple[Node, Node]] = []
        # self.edges: List[Tuple[Node, Node]] = []
        # self.edges: Dict[str, List[Tuple[Node, Node]]] = {}
        # self._map_to_nodes: Dict[str, Node] = {}

        self.nodes : Dict[str, Node] = {}

    def _add_node(self, node: Node):
        # self._map_to_nodes[node.descriptor] = node

        # rc = reverse_complement(node.descriptor)
        # rc_node = Node(rc)
        # self._map_to_nodes[rc] = rc_node

        # self.blocks.append((node, rc_node))

        self.nodes[node.descriptor] = node
        
        rc = reverse_complement(node.descriptor)
        rc_node = Node(rc)
        self.nodes[rc] = rc_node

    def _get_twin(self, node: Node):
        rc = reverse_complement(node.descriptor)
        return self.nodes[rc]
    
    def _add_edge(self, node_a, node_b):
        # self.edges.append((node_a, node_b))
        # self.edges.append((node_b, node_a))
        # if node_a.descriptor not in self.edges:
        #     self.edges[node_a.descriptor] = []
        # if node_b.descriptor not in self.edges: 
        #     self.edges[node_b.descriptor] = []

        node_a.out_edges.append(node_b)
        node_b.in_edges.append(node_a)

        if node_a.descriptor != "-1":
            node_a.in_edges.append(node_b)
            node_b.out_edges.append(node_a)
        # node_b.in_edges.append((node_a))
        
        # self.edges[node_a.descriptor].append((node_a, node_b))
        # self.edges[node_b.descriptor].append((node_b, node_a))
    
    def _has_single_edge(self):
        # TODO: make this more efficient
        # incoming_edges = 0

        # for (descriptor, _) in self._map_to_nodes.items():
        #     if descriptor == node.descriptor: continue

        #     for (_, end) in self.edges[descriptor]:
        #         if end.descriptor == node.descriptor: incoming_edges += 1
        #         if incoming_edges > 1:
        #             return False
        # return True
        return len(self.in_edges) == 1
    
    def _merge_nodes(self, node_a, node_b):
        twin = self._get_twin(node_b)

        new_descriptor = node_a.descriptor + twin.descriptor[-1]
        print('new descriptor is', new_descriptor)


        # node_a.out_edges = node_b.out_edges
        outgoing_edges: List[Node] = []
        for n in node_a.out_edges:
            if n.descriptor != node_b.descriptor: outgoing_edges.append(n)
        for n in node_b.out_edges:
            if n.descriptor != node_a.descriptor: outgoing_edges.append(n)
        # node_a.in_edges = [n for n in node_a.in_edges if n.descriptor != node_b.descriptor]
        node_a.out_edges = outgoing_edges

        incoming_edges: List[Node] = []
        for n in node_a.in_edges:
            if n.descriptor != node_b.descriptor: incoming_edges.append(n)
        for n in node_b.in_edges:
            if n.descriptor != node_a.descriptor: incoming_edges.append(n)
        
        node_a.in_edges = incoming_edges
        
        del self.nodes[node_b.descriptor]
        del self.nodes[node_a.descriptor]

        node_a.descriptor = new_descriptor
        self.nodes[new_descriptor] = node_a

        # set node a's descriptor to be the new descriptor
        # avoid having to delete all node x -> node a and node a -> node x edges

        # get all node_b -> node y edges

        # delete previous entries in edges, _map_to_nodes
        return

def create_pre_nodes(reads, kmer_table, hash_length, graph):
    for i in range(len(reads)):
        seq = reads[i]
        start = 0

        prev_node = graph.source
        
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

                    twin = graph._get_twin(node)

                    graph._add_edge(prev_node, twin)
                    prev_node = node
                    # pre_nodes.append(node)
                
                # create a node for the overlapping kmer
                if (i, dir, end - hash_length + 1) == first_occurrence:
                    node = Node(new_kmer)
                    # pre_nodes.append(node)
                    graph._add_node(node)
                    twin = graph._get_twin(node)

                    graph._add_edge(prev_node, twin)
                    prev_node = node
                
                else:
                    curr_node = graph.nodes[new_kmer]
                    twin = graph._get_twin(curr_node)
                    graph._add_edge(prev_node, twin)
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
    # go through nodes

    # if node a has only one outgoing edge to node b and
    # node b has only one incoming edge

    runs = [[]]

    for n in graph.nodes.values():
        print(n)
        # print(n, n.out_edges, n.in_edges)
        if n._has_single_out_edge():
            other_end = n.out_edges[0]

            if other_end._has_single_in_edge():
                # print('concatenating', n, other_end)
                # graph._merge_nodes(n, other_end)
                if len(runs[-1]) == 0: runs[-1].append(n)
                runs[-1].append(other_end)
            else:
                runs.append([])

    print(runs)

    return