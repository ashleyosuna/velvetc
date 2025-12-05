from graph import Graph, Node
from tourbus import tourbus
from utils import n50

def test_initial_node_and_edges():
    kmers = ["ATC", "ATG", "GCA", "CAA", "AAA", "AAT"]

    reads = [
        "GATGCA", "GCAAA", "AAATT"
    ]

    g = Graph(3)
    g.create_init_nodes(kmers)

    g.map_through_reads(reads)

    g.concatenate_nodes()

    g.clip_tips()

    tourbus(g)

    g.concatenate_nodes()

    contigs = g.get_contigs()
    n50_res = n50(contigs, g.k)

    print(contigs, n50_res)

def test_tour_bus():
    g = Graph(3)

    g.nodes = {
        1: Node("ATC", 1),
        -1: Node("GAT", -1),
        2: Node("TCA", 2),
        -2: Node("TGA", -2),
        3: Node("CAA", 3),
        -3: Node("TTG", -3),
        4: Node("AAT", 4),
        -4: Node("ATT", -4),
        5: Node("TCC", 5),
        -5: Node("GGA", -5),
        6: Node("CCA", 6),
        -6: Node("TGG", -6),
        7: Node("CAA", 7),
        -7: Node("TTG", -7)
    }

    g.map_to_nodes = {
        "ATC": 1, "GAT": -1,
        "TCA": 2, "TGA": -2,
        "CAA": 3, "TTG": -3,
        "AAT": 4, "ATT": -4,
        "TCC": 5, "GGA": -5,
        "CCA": 6, "TGG": -6,
        "CAA": 7, "TTG": -7
    }

    g.nodes[1].out_edges = {2:1, 5:1}
    g.nodes[-1].in_edges = {-2:1, -5:1}

    g.nodes[2].out_edges = {3:1}
    g.nodes[2].in_edges = {1:1}
    g.nodes[-2].out_edges = {-1:1}
    g.nodes[-2].in_edges = {-3:1}

    g.nodes[3].out_edges = {4:1}
    g.nodes[3].in_edges = {2:1}
    g.nodes[-3].in_edges = {-4:1}
    g.nodes[-3].out_edges = {-2:1}


    g.nodes[4].in_edges = {3:1, 7:1}
    g.nodes[-4].out_edges = {-3:1, -7:1}

    g.nodes[5].in_edges = {1:1}
    g.nodes[5].out_edges = {6:1}
    g.nodes[-5].in_edges = {-6:1}
    g.nodes[-5].out_edges = {-1:1}

    g.nodes[6].in_edges = {5:1}
    g.nodes[6].out_edges = {7:1}
    g.nodes[-6].in_edges = {-7:1}
    g.nodes[-6].out_edges = {-5:1}

    g.nodes[7].in_edges = {6:1}
    g.nodes[7].out_edges = {4:1}
    g.nodes[-7].in_edges = {-4:1}
    g.nodes[-7].out_edges = {-6:1}

    g.starts = [1]

    # g.tour_bus()
    tourbus(g)

    for id, node in g.nodes.items():
        print(id, node, node.in_edges, node.out_edges)


test_initial_node_and_edges()
# test_tour_bus()