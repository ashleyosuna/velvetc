import graph
from tourbus import tourbus

def test_tour_bus():
    print()
    g = graph.Graph(hash_length=3)

    g.nodes = {
        1: graph.Node("ATC", 1),
        -1: graph.Node("GAT", -1),
        2: graph.Node("TCA", 2),
        -2: graph.Node("TGA", -2),
        3: graph.Node("CAA", 3),
        -3: graph.Node("TTG", -3),
        4: graph.Node("AAT", 4),
        -4: graph.Node("ATT", -4),
        5: graph.Node("TCC", 5),
        -5: graph.Node("GGA", -5),
        6: graph.Node("CCA", 6),
        -6: graph.Node("TGG", -6),
        7: graph.Node("CAA", 7),
        -7: graph.Node("TTG", -7)
    }

    g.map_to_id = {
        "ATC": 1, "GAT": -1,
        "TCA": 2, "TGA": -2,
        "CAA": 3, "TTG": -3,
        "AAT": 4, "ATT": -4,
        "TCC": 5, "GGA": -5,
        "CCA": 6, "TGG": -6,
        "CAA": 7, "TTG": -7
    }

    g.nodes[1].out_edges = [graph.Edge(2), graph.Edge(5)]
    g.nodes[-1].in_edges = [graph.Edge(-2), graph.Edge(-5)]

    g.nodes[2].out_edges = [graph.Edge(3)]
    g.nodes[2].in_edges = [graph.Edge(1)]

    g.nodes[-2].out_edges = [graph.Edge(1)]
    g.nodes[-2].in_edges = [graph.Edge(-3)]

    g.nodes[3].out_edges = [graph.Edge(4)]
    g.nodes[3].in_edges = [graph.Edge(2)]

    g.nodes[-3].in_edges = [graph.Edge(-4)]
    g.nodes[-3].out_edges = [graph.Edge(-2)]

    g.nodes[4].in_edges = [graph.Edge(3), graph.Edge(7)]
    g.nodes[-4].out_edges = [graph.Edge(-3), graph.Edge(-7)]

    g.nodes[5].in_edges = [graph.Edge(1)]
    g.nodes[5].out_edges = [graph.Edge(6)]
    g.nodes[-5].in_edges = [graph.Edge(-6)]
    g.nodes[-5].out_edges = [graph.Edge(-1)]

    g.nodes[6].in_edges = [graph.Edge(5)]
    g.nodes[6].out_edges = [graph.Edge(7)]
    g.nodes[-6].in_edges = [graph.Edge(-7)]
    g.nodes[-6].out_edges = [graph.Edge(-5)]

    g.nodes[7].in_edges = [graph.Edge(6)]
    g.nodes[7].out_edges = [graph.Edge(4)]
    g.nodes[-7].in_edges = [graph.Edge(-4)]
    g.nodes[-7].out_edges = [graph.Edge(-6)]

    g.starts = [1]

    # g.tour_bus()
    tourbus(g)

test_tour_bus()