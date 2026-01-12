import itertools
import random

import networkx as nx
import matplotlib.pyplot as plt
from pprint import pprint

#edge list plotting function
def plot_graph(edge_list):
    G = nx.MultiGraph()
    G.add_edges_from(edge_list)

    pos = nx.spring_layout(G)
    fig, ax = plt.subplots(figsize=(8, 6))

    nx.draw_networkx_nodes(G, pos, node_color="lightblue", node_size=500)
    nx.draw_networkx_labels(G, pos)

    for edge in G.edges(data=True, keys=True):
        u, v, key = edge[0], edge[1], edge[2]

        ax.annotate(
            "",
            xy=pos[u],
            xycoords="data",
            xytext=pos[v],
            textcoords="data",
            arrowprops=dict(
                arrowstyle="-",
                color="black",
                connectionstyle=f"arc3,rad={0.2 * key}",
                shrinkA=10,
                shrinkB=10,
            ),
        )
    plt.show()

#builders graph with element property attached
def build_graph(edges):
    G = nx.MultiGraph()
    G.add_edges_from(edges)

    for node in G:
        G.nodes[node]["element"] = node[0]
    return G

#checks if two graphs are isomorphs based on element property
def are_isomorphic(edges1, edges2):
    G1 = build_graph(edges1)
    G2 = build_graph(edges2)

    def match(n1, n2):
        return n1["element"] == n2["element"]

    return nx.is_isomorphic(G1, G2, node_match=match)

#creates a compound with valences and frequency
def build_compound(elements, valences):
    compound_valence = {}
    for el in elements.keys():
        number = elements[el]
        for i in range(1, number + 1):
            compound_valence[f"{el}{i}"] = valences[el]
    return compound_valence

#check if every valence is 0 
def balanced(valence):
    return all([val == 0 for val in valence.values()])


valences = {"C": 4, "O": 2, "H": 1, "P": 3}
compound = build_compound({"C": 2,"H":6,"O":1}, valences)


def unique_isomers(
    compound_valence, MAX_BONDS=3
):  # takes elements + valences, and max bonds
    combinations = list(
        itertools.combinations(compound_valence.keys(), 2)
    )  # make combinations of all possible bonds between elements
    result = []  # results with isomorphic graphs (same shape)
    u_result = []  # results without isomorphic graphs

    def backtracking(
        index, current_edges, current_valences
    ):  # index represents which combination function is on

        if index == len(
            combinations
        ):  # base case if you've finished all the combinations
            G = nx.MultiGraph()
            G.add_edges_from(current_edges)
            components = nx.number_connected_components(G)
            remaining_bonds = sum(current_valences.values()) // 2 

            if remaining_bonds < components - 1: # checks if molecule is a disconnected graph, number of connections n-1, number of nodes should be n
                print("disconnected")
                return

            if balanced(
                current_valences
            ):  # check if valences are = 0 
                result.append(current_edges[:])  # add valid isomer to results
            return  # end current function

        u, v = combinations[index]  # two atoms (nodes) of the current combination
        backtracking(
            index + 1, current_edges, current_valences
        )  # makes 0 bonds with the combination. This function will stack till you reach an empty []. Then you will backtrack.

        max_bonds = min(
            current_valences[u], current_valences[v], MAX_BONDS
        )  # find max bonds by getting min valence of element: e.g. (C,H) max bond is min(4,1) = 1

        for amount_b in range(
            1, max_bonds + 1
        ):  # make 1 ->  max_bonds with the combination. This function will stack and backtrack

            current_valences[u] -= amount_b  # updates valences
            current_valences[v] -= amount_b  # updates valences
            current_edges.extend(
                [(u, v)] * amount_b
            )  # adds combination n number of times based on iteration of loop

            backtracking(
                index + 1, current_edges, current_valences
            )  # checks next combination 
            for _ in range(amount_b):  # undo's last move by removing edges (BACKTRACK)
                current_edges.pop()
            current_valences[u] += amount_b  # restores valence
            current_valences[v] += amount_b  # restores valence

    backtracking(0, [], compound_valence)  # Invokes helper function

    for mol in result:  # creates new list with no graphical isomorphs
        is_new = True
        for umol in u_result:
            if are_isomorphic(umol, mol):
                is_new = False
        if is_new:
            u_result.append(mol)

    pprint(u_result)
    print(len(u_result))
    return u_result


for mol in unique_isomers(compound):
    plot_graph(mol)  # plots graph for each molecule
