import matplotlib.pyplot as plt
import math as maths
import copy
import numpy
from collections import deque
import networkx as nx
from networkx.algorithms import isomorphism
import time

VALENCES = {"C": 4, "O": 2, "N": 3, "H": 1, "P": 3}


class Molecule:
    def __init__(self, atoms):
        self.atoms = atoms
        self.num_atoms = len(self.atoms)
        self.max_valence = numpy.array([VALENCES[i] for i in self.atoms])
        self.adj_matrix = numpy.zeros((self.num_atoms, self.num_atoms))

    def _valid_cell(self, i, j):
        return i < j

    def set_bond(self, i, j, w=1):
        if self._valid_cell(i, j):
            self.adj_matrix[i, j] = w
            self.adj_matrix[j, i] = w

        else:
            raise IndexError

    def get_used_valence(self):
        used_val = self.adj_matrix.sum(axis=1)
        return used_val

    def get_free_valences(self):
        return self.max_valence - self.get_used_valence()

    def get_free_valence_at(self, row):
        return self.max_valence[row] - int(numpy.sum(self.adj_matrix[row, :]))

    def is_valid(self):
        return all(self.get_free_valences() >= 0)

    def is_complete(self):
        if self.is_valid():
            return sum(self.get_free_valences()) == 0

    def get_max_valence(self, r):
        return VALENCES[self.atoms[r]]

    def is_correct(self):
        return self.is_connected() and self.is_complete()

    def is_connected(self):
        visited = set()
        stack = deque()
        visited.add(0)
        stack.append(0)
        while len(stack) > 0:
            top = stack.pop()
            # print(top)
            for j in range(self.num_atoms):
                cell = self.adj_matrix[top, j]
                if j not in visited and cell > 0:
                    visited.add(j)
                    stack.append(j)
        return len(visited) == self.num_atoms


valid = []


def build_compound(atoms: dict):
    mol = []
    for k, v in atoms.items():
        for _ in range(v):
            mol.append(k)
    return mol


def triangleArea(n):
    return (n**2 - n) // 2


def generateIsomers(mol, max_bond):
    results = []

    def generateIsomerRec(mol: Molecule, i, max_bond):

        if not mol.is_valid():
            return

        if i == triangleArea(mol.num_atoms):
            if mol.is_correct():
                results.append(copy.deepcopy(mol))
            return

        row = None
        for j in range(1, triangleArea(mol.num_atoms) + 1):
            if i - (triangleArea(mol.num_atoms) - triangleArea(mol.num_atoms - j)) < 0:
                row = j - 1
                break
        total_to_row = triangleArea(mol.num_atoms) - triangleArea(
            mol.num_atoms - row - 1
        )
        col = (mol.num_atoms) - (total_to_row - i)

        max_bond_check = min(
            mol.get_free_valence_at(row), mol.get_free_valence_at(col), max_bond
        )

        for v in range(0, max_bond_check + 1):
            old = mol.adj_matrix[row, col]
            mol.set_bond(row, col, v)
            generateIsomerRec(mol, i + 1, max_bond)
            mol.set_bond(row, col, old)

    generateIsomerRec(mol, 0, max_bond)
    return results


def to_graph(mol: Molecule):
    G = nx.from_numpy_array(mol.adj_matrix)
    for i, atom in enumerate(mol.atoms):
        G.nodes[i]["element"] = atom
    return G


def is_mol_isomorphic(G1, G2):
    node_match = isomorphism.categorical_node_match("element", None)
    edge_match = isomorphism.categorical_edge_match("weight", 1)
    return nx.is_isomorphic(G1, G2, node_match=node_match, edge_match=edge_match)


def draw(G):
    MG = nx.MultiGraph()
    MG.add_nodes_from(G.nodes(data=True))

    for u, v, data in G.edges(data=True):
        for _ in range(int(data.get("weight", 1))):
            MG.add_edge(u, v)

    pos = nx.spring_layout(MG, seed=42)
    labels = nx.get_node_attributes(MG, "element")

    groups = {}
    for u, v, key in MG.edges(keys=True):
        pair = frozenset((u, v))
        groups.setdefault(pair, []).append((u, v, key))

    for edges in groups.values():
        n = len(edges)
        styles = [f"arc3,rad={(i - (n - 1) / 2) * 0.2}" for i in range(n)]

        nx.draw_networkx_edges(
            MG, pos, edgelist=edges, connectionstyle=styles, arrows=True, arrowstyle="-"
        )

    nx.draw_networkx_nodes(MG, pos, node_size=500)
    nx.draw_networkx_labels(MG, pos, labels=labels, font_size=10, font_color="white")

    plt.axis("off")
    plt.show()


def alkane(n):
    return Molecule(build_compound({"C": n, "H": 2 * n + 2}))


def alkene(n):
    return Molecule(build_compound({"C": n, "H": 2 * n}))


test_comp = build_compound({"P": 4})
mol = alkene(3)
start = time.time()
valid = generateIsomers(mol, 3)
end = time.time()
print(end - start)

unique_graphs = []

for m in valid:
    G = to_graph(m)
    iso_flag = False
    for g in unique_graphs:
        if is_mol_isomorphic(g, G):
            iso_flag = True
            break
    if not iso_flag:
        unique_graphs.append(G)

for g in unique_graphs:
    draw(g)
