#!/bin/env python3

import itertools
import collections
import glob
from pprint import pprint
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
import networkx as nx


def generate_fragments(mol, min_num_atoms=5, max_num_atoms=0):
    
    # create molecule graph
    graph = nx.Graph()
    graph.add_edges_from([
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
    ])

    if not nx.is_connected(graph):
        return
    
    visited = set()
    
    # find subgraphs
    for size in range(min_num_atoms, max_num_atoms or len(graph.nodes)):
        for nodes in itertools.combinations(graph, size):
            subgraph = graph.subgraph(nodes)
            if nx.is_connected(subgraph):
                # create fragment molecule from subgraph
                frag = Chem.RWMol(mol)
                
                to_be_removed = set(graph.nodes) - set(subgraph.nodes)
                for node in sorted(to_be_removed, reverse=True):
                    frag.RemoveAtom(node)
                
                smi = Chem.MolToSmiles(frag)
                if smi in visited:
                    continue
                visited.add(smi)
                frag = Chem.MolFromSmiles(smi)
                if not frag:
                    continue
                yield frag, smi


import sys
import subprocess


def print_fragments_from_molecules(filename):
    try:
        mol = rdmolfiles.MolFromPDBFile(filename)
        for _, smi in generate_fragments(mol):
            print(smi)
    except:
        print("Failed",  filename)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        subprocess.run("""
            ls murray_ligands/lig_*.pdb | head -n104 |
                xargs -I@ python fragmentator.py @ |
                sort | uniq -c > output.txt
        """, shell=True)
    
    elif len(sys.argv) == 2:
        print_fragments_from_molecules(sys.argv[1])
