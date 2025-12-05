from utils import reverse_complement, canonical_form
from typing import List
import numpy as np

class Node:
    seq: str
    id: int
    in_edges: dict[int, int]
    out_edges: dict[int, int]

    def __init__(self, seq, id):
        self.seq = seq
        self.id = id
        self.in_edges = {}
        self.out_edges = {}
    
    def __repr__(self):
        return f"{self.seq}"
    
    def add_incoming_edge(self, from_node, multiplicity):
        if self.in_edges.get(from_node, -1) == -1:
            self.in_edges[from_node] = multiplicity
        else:
            self.in_edges[from_node] += multiplicity
    
    def add_outgoing_edge(self, to_node, multiplicity):
        if self.out_edges.get(to_node, -1) == -1:
            self.out_edges[to_node] = multiplicity
        else:
            self.out_edges[to_node] += multiplicity

    def simple_edge(self):
        if len(self.out_edges) > 1: return False

        if len(self.out_edges) == 0:
            # print("shouldn't get here")
            return False
        
        other_end = list(self.out_edges.keys())[0]
        return abs(other_end) != abs(self.id)
    
    def length(self, k):
        # return number of kmers
        return len(self.seq) - k + 1
    
    def compare(self, other):
        """
        Simple comparison of nodes; if more than 50% of bases match then return True,
        otherwise nodes are not considered similar enough.
        """
        u = self.seq
        v = other.seq

        n, m = len(u) + 1, len(v) + 1
        network = np.full((n, m), 0, dtype=float)

        # fill first row
        for i in range(n): network[i][0] = 0 * i

        # fill first column
        for j in range(m): network[0][j] = 0 * j

        # fill matrix row by row
        for i in range(1, n):
            for j in range(1, m):
                network[i][j] = max(
                    network[i-1][j],
                    network[i][j-1],
                    network[i-1][j-1] + (1 if u[i-1] == v[j-1] else 0)
                )

        score = network[i-1][j-1]

        if score > max(len(u), len(v)) / 2: return True
        return False

class Graph:
    def __init__(self, k):
        self.k: int = k
        self.nodes: dict[int, Node] = {}
        self.map_to_nodes: dict[str, int] = {}
        self.starts: List[int] = []

    @property
    def node_count(self):
        return len(self.nodes) / 2
    
    def add_node(self, node):
        self.nodes[node.id] = node
        self.map_to_nodes[node.seq] = node.id

        rc = reverse_complement(node.seq)

        twin_node = Node(rc, -node.id)
        self.nodes[twin_node.id] = twin_node
        self.map_to_nodes[rc] = twin_node.id

    def create_init_nodes(self, reads, kmers):
        """
        Creates initial nodes and their twins as given by the canonical kmers.
        """
        # next_id = 1
        
        # for kmer in kmers:
        #     # creating canonical node
        #     self.nodes[next_id] = Node(kmer, next_id)
        #     self.map_to_nodes[kmer] = next_id

        #     # creating twin node
        #     rc = reverse_complement(kmer)
        #     self.nodes[-next_id] = Node(rc, -next_id)
        #     self.map_to_nodes[rc] = -next_id
            
        #     next_id += 1
        next_id = 1
        for i in range(len(reads)):
            seq = reads[i]
            start = 0

            prev_node = None
        
            # get initial kmer
            new_kmer = seq[start:self.k]

            # initialize sequence of uninterrupted kmers
            consecutive_seq = seq[start:self.k - 1]

            for end in range(self.k - 1, len(seq)):
                # if not in initial kmer, slide kmer window
                if end >= self.k: new_kmer = new_kmer[1:] + seq[end]
                
                can_kmer, dir = canonical_form(new_kmer)
                first_occurrence = kmers[can_kmer][0]

                # if newly added kmer to the window overlaps with other reads
                if len(kmers[can_kmer]) > 1:
                    # create a new node for the previously uninterrupted sequence of kmers
                    if (end - 1) - start + 1 >= self.k:
                        node = Node(consecutive_seq, next_id)
                        self.add_node(node)

                        # twin = graph._get_twin(node)


                        # graph._add_edge(prev_node, twin)
                        # prev_node = node
                        # pre_nodes.append(node)

                        self.add_edge(prev_node, node)
                        prev_node = node

                        next_id += 1
                    
                    # create a node for the overlapping kmer
                    if (i, dir, end - self.k + 1) == first_occurrence:
                        node = Node(new_kmer, next_id)
                        # pre_nodes.append(node)
                        self.add_node(node)
                        # twin = graph._get_twin(node)

                        self.add_edge(prev_node, node)
                        prev_node = node
                        next_id += 1
                    
                    else:
                        curr_node = self.nodes[self.map_to_nodes[new_kmer]]
                        # twin = graph._get_twin(curr_node)
                        self.add_edge(prev_node, curr_node)
                        prev_node = curr_node
                    
                    # slide window
                    start = end - self.k + 2
                    consecutive_seq = consecutive_seq[1:]
                
                consecutive_seq += seq[end]
                
                # if we have reached the end of the read and there is no overlap create a node
                # for this rightmost run of uninterrupted kmers
                if end == len(seq) - 1 and len(kmers[can_kmer]) == 1:
                    node = Node(consecutive_seq, next_id)
                    # pre_nodes.append(node)
                    self.add_node(node)
                    # twin = graph._get_twin(node)

                    self.add_edge(prev_node, node)
                    prev_node = node
                    next_id += 1

    def add_edge(self, node_a, node_b):
        # print(node_a, node_b)
        if node_a is None: 
            self.starts.append(node_b.id)
            return
        # add outgoing edge to b into node a
        node_a.add_outgoing_edge(node_b.id, 1)

        # add incoming edge from a into b
        node_b.add_incoming_edge(node_a.id, 1)
        
        twin_a = self.nodes[-node_a.id]
        twin_b = self.nodes[-node_b.id]

        # add outgoing edge to a's twin into b's twin
        twin_b.add_outgoing_edge(twin_a.id, 1)
        
        # add incoming edge from b's twin into a
        twin_a.add_incoming_edge(twin_b.id, 1)
    
    def map_through_reads(self, reads):
        """
        Goes through the reads and adds edges between corresponding nodes.
        """
        k = self.k
        
        for read in reads:
            if len(read) <= k: continue

            prev_kmer = read[:k]

            self.starts.append(self.map_to_nodes[prev_kmer])
            
            curr_kmer = read[1:k]

            for i in range(k, len(read)):
                if i > k: curr_kmer = curr_kmer[1:]
                curr_kmer += read[i]

                prev_node = self.nodes[self.map_to_nodes[prev_kmer]]
                curr_node = self.nodes[self.map_to_nodes[curr_kmer]]

                self.add_edge(prev_node, curr_node)

                prev_kmer = curr_kmer

    def remove_from_starts(self, node_id):
        filtered_starts = []
        for s in self.starts:
            if s != node_id: filtered_starts.append(s)
        self.starts = filtered_starts

    def concatenate_two_nodes(self, node_a, node_b):
        # node a only has one outgoing edge, so simply inherit all node b's edges
        node_a.out_edges = node_b.out_edges

        # all nodes with incoming edges from b should now have incoming edges from a
        for other_end, mult in node_b.out_edges.items():
            node = self.nodes[other_end]
            node.in_edges[node_a.id] = mult
            del node.in_edges[node_b.id]
        
        twin_a = self.nodes[-node_a.id]
        twin_b = self.nodes[-node_b.id]

        # twin a must have only one incoming edge, so simply inherit all twin b's incoming edges
        twin_a.in_edges = twin_b.in_edges

        # all nodes with outgoing edges into twin b should now point to twin a
        for other_end, mult in twin_b.in_edges.items():
            node = self.nodes[other_end]
            node.out_edges[twin_a.id] = mult
            del node.out_edges[twin_b.id]

        # delete old sequences from map
        del self.map_to_nodes[node_a.seq]
        del self.map_to_nodes[twin_a.seq]
        
        # modify sequences
        node_a.seq += node_b.seq[self.k-1:]
        twin_a.seq = twin_b.seq + twin_a.seq[self.k-1:]

        # add updated sequences to map
        self.map_to_nodes[node_a.seq] = node_a.id
        self.map_to_nodes[twin_a.seq] = twin_a.id

        # delete node b
        del self.map_to_nodes[node_b.seq]
        del self.nodes[node_b.id]

        # delete node b's twin
        del self.map_to_nodes[twin_b.seq]
        del self.nodes[twin_b.id]

        self.remove_from_starts(node_b.id)
        self.remove_from_starts(twin_b.id)

    def concatenate_nodes(self):
        modified = True

        while modified:
            modified = False

            for node in self.nodes.values():
                if len(node.out_edges) == 1:
                    other_end = self.nodes[list(node.out_edges.keys())[0]]

                    if abs(node.id) != abs(other_end.id) and len(other_end.in_edges) == 1:
                        
                        self.concatenate_two_nodes(node, other_end)
                        modified = True
                
                if modified: break
        return
    
    def find_start_points(self):
        start_points = []

        for node in self.nodes.values():
            if len(node.in_edges) == 0: start_points.append(node.id)
        
        return start_points
    
    def destroy_tip(self, tip):
        """
        Destroys all nodes in the tip, along with their twin nodes.
        Update edges for node connected to last node in the tip.
        """
        last_node = self.nodes[tip[-1]]
        last_twin = self.nodes[-tip[-1]]

        # unlink outgoing and incoming edges from last node in the tip
        for other_end in last_node.out_edges.keys():
            node = self.nodes[other_end]
            del node.in_edges[last_node.id]

        for other_end in last_node.in_edges.keys():
            node = self.nodes[other_end]
            del node.out_edges[last_node.id]

        # unlink outgoing and incoming edges from last twin in the tip
        for other_end in last_twin.out_edges.keys():
            node = self.nodes[other_end]
            del node.in_edges[last_twin.id]

        for other_end in last_twin.in_edges.keys():
            node = self.nodes[other_end]
            del node.out_edges[last_twin.id]

        # delete entries for all nodes in the tip
        for node in tip:
            node_obj = self.nodes[node]
            twin = self.nodes[-node]

            del self.map_to_nodes[node_obj.seq]
            del self.map_to_nodes[twin.seq]

            del self.nodes[node]
            del self.nodes[twin.id]
        
            self.remove_from_starts(node)
            self.remove_from_starts(-node)
        return
    
    def clip_tips(self):
        modified = True

        while modified:
            modified = False
            for s in self.starts:
                curr_node = self.nodes[s]

                tip = []
                length = 0

                while curr_node.simple_edge():
                    tip.append(curr_node.id)                 
                    length += curr_node.length(self.k)

                    other_end = list(curr_node.out_edges.keys())[0]
                    curr_node = self.nodes[other_end]

                if length > 0 and length < 2 * self.k:
                    self.destroy_tip(tip)
                    modified = True

                if modified: break
            
            # if modified: start_points = start_points[1:]        
        
        return
    
    def merge_nodes(self, node_a, node_b):
        """
        Merges nodes that were determined to be similar enough when merging paths.
        Node a should inherit node b's edges.
        """
        for other_end, mult in node_b.in_edges.items(): 
            node_a.add_incoming_edge(other_end, mult)

            # remove pointers to b and re-map to node a
            other_end_node = self.nodes[other_end]
            # other_end.replace_edge(node_b.id, node_a.id, False)
            del other_end_node.out_edges[node_b.id]
            other_end_node.add_outgoing_edge(node_a.id, mult)

        for other_end, mult in node_b.out_edges.items(): 
            node_a.add_outgoing_edge(other_end, mult)

            other_end_node = self.nodes[other_end]
            del other_end_node.in_edges[node_b.id]
            other_end_node.add_incoming_edge(node_a.id, mult)

        twin_a = self.nodes[-node_a.id]
        twin_b = self.nodes[-node_b.id]

        for other_end, mult in twin_b.in_edges.items(): 
            twin_a.add_incoming_edge(other_end, mult)

            other_end_node = self.nodes[other_end]
            del other_end_node.out_edges[twin_b.id]
            other_end_node.add_outgoing_edge(twin_a.id, mult)
        
        for other_end, mult in twin_b.out_edges.items(): 
            twin_a.add_outgoing_edge(other_end, mult)

            other_end_node = self.nodes[other_end]
            del other_end_node.in_edges[twin_b.id]
            other_end_node.add_incoming_edge(twin_a.id, mult)

        # delete all entries to node b and its twin
        del self.nodes[node_b.id]
        del self.nodes[-node_b.id]
        del self.map_to_nodes[node_b.seq]
        del self.map_to_nodes[twin_b.seq]
        
        self.remove_from_starts(node_b.id)
        self.remove_from_starts(twin_b.id)
        return

    def merge_paths(self, path1, path2):
        """
        Merges two paths forming a bubble that were considered similar enough.
        Path 1 should be the path that was reached first, i.e., the higher coverage path.
        """
        path2_ptr = 0
        for node_id in path1:
            path1_node = self.nodes[node_id]
            path2_node = self.nodes[path2[path2_ptr]]
            
            similar = path1_node.compare(path2_node)
            if similar:
                # print(f"node {node_id} will be merged with {path2_node.id}")
                self.merge_nodes(path1_node, path2_node)
                path2_ptr += 1
        return
    
    def get_contigs(self):
        """
        Outputs long nodes in the graph.
        """
        contigs = {}
        for id, node in self.nodes.items():
            if id < 0 or node.length(self.k) < self.k:
                continue
            # contigs.append(node.seq)
            contigs[node.id] = node.seq
        return contigs
