"""Plotting the reaction enzymes."""
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from netgraph import Graph

rate_constant_dict = {"k1": 100, "k2": 2, "k11": 10, "k4": 2, "k3": 100, "k33": 10}

G = nx.DiGraph()
edge_list = [
    (r"$S_{1} + E_{1}$", r"$C_{1}$", {"w": "k1"}),
    (r"$C_{1}$", r"$S_{1} + E_{1}$", {"w": "k11"}),
    (r"$C_{1}$", r"$P_{1} + E_{2}$", {"w": "k2"}),
    (r"$P_{1} + E_{2}$", r"$C_{2}$", {"w": "k3"}),
    (r"$C_{2}$", r"$P_{1} + E_{2}$", {"w": "k33"}),
    (r"$C_{2}$", r"$P_{2}$", {"w": "k4"}),
]

G.add_edges_from(edge_list)
pos = {
    "$S_{1} + E_{1}$": np.array([0.0, 0.5]),
    r"$C_{1}$": np.array([0.3, 0.5]),
    r"$P_{1} + E_{2}$": np.array([0.5, 0.5]),
    r"$C_{2}$": np.array([0.75, 0.5]),
    r"$P_{2}$": np.array([1, 0.5]),
}
edge_labels = {
    (r"$S_{1} + E_{1}$", r"$C_{1}$"): r"$k_{1}: $" + str(rate_constant_dict["k1"]),
    (r"$C_{1}$", r"$S_{1} + E_{1}$"): r"$k_{-1}: $" + str(rate_constant_dict["k11"]),
    (r"$C_{1}$", r"$P_{1} + E_{2}$"): r"$k_{2}: $" + str(rate_constant_dict["k2"]),
    (r"$P_{1} + E_{2}$", r"$C_{2}$"): r"$k_{3}: $" + str(rate_constant_dict["k3"]),
    (r"$C_{2}$", r"$P_{1} + E_{2}$"): r"$k_{-3}: $" + str(rate_constant_dict["k33"]),
    (r"$C_{2}$", r"$P_{2}$"): r"$k_{4}: $" + str(rate_constant_dict["k4"]),
}

fig, ax = plt.subplots(figsize=(30, 10))
plot_instance = Graph(
    G,
    node_labels=True,
    edge_labels=edge_labels,
    edge_label_fontdict=dict(size=14),
    # edge_label_position=0.6,
    node_layout=pos,
    edge_layout="curved",
    node_size=6,
    edge_width=0.5,
    arrows=True,
    ax=ax,
)
plt.title(r"$S_{1}$ Input S2")
plt.savefig("Case1.png")
