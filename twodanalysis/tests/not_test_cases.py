import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pytest
import time


#from twodanalysis import Memb2D
#import imageio
import os
import time
import networkx as nx

from twodanalysis import MembProp
from twodanalysis import Cumulative2D
from twodanalysis import Voronoi2D
from twodanalysis import analysis


from twodanalysis.data.files import MEMBRANE_TPR, MEMBRANE_XTC


tpr = MEMBRANE_TPR
traj = MEMBRANE_XTC



u = mda.Universe(MEMBRANE_TPR, MEMBRANE_XTC, )
import re

def natural_key(node_name):
    """
    Key function to sort strings like C1, C2, C10 numerically.
    Extracts the prefix and numeric part.
    """
    match = re.match(r"([A-Za-z]+)(\d+)", str(node_name))
    if match:
        prefix, number = match.groups()
        return (prefix, int(number))
    else:
        return (node_name, 0)  # fallback


nbins = 150
start = 50
final = 100
layer = "top"  # or "bot", depending on the layer you want to analyze
membrane = Voronoi2D(u, nbins = nbins)
lipid_list = membrane.lipid_list.copy()
lipid_list.remove("CHL1")
print(lipid_list)
carbons, scd, error = analysis.OrderParameters.window_scd(u, selection="resname POPE", start=0, final=77, window = 10, chain ="sn2", step = 1)
plt.errorbar(carbons, scd, yerr=error, fmt='o', label='SCD Order Parameter')
plt.show()

print(scd)


#membrane._store_fist_lipids()
#lipids_ids = membrane.first_lipids
"""
for key in lipid_list:
    print(membrane.connection_chains[key])
    _, chain0 = membrane.cut_structure(
                    membrane.memb.select_atoms(f"resid {lipids_ids[key]}"),
                    membrane.connection_chains[key][0],
                    )
    _, chain1 =     membrane.cut_structure(
                    membrane.memb.select_atoms(f"resid {lipids_ids[key]}"),
                    membrane.connection_chains[key][0],
                    )


    c_in_node = [n for n in chain0.nodes if "C" in str(n)]
    c_in_node = sorted(c_in_node, key=natural_key)
    neighbors = [list(chain0.neighbors(n)) for n in c_in_node]
    neighbors = [[c_in_node[i]] + [item for item in neighbors[i] if "C" not in item] for i in range(len(c_in_node))]

"""


    #temp_chain0 = chain0.copy()
    #temp_chain0.remove_nodes_from(nodes_to_remove)

    #pos = nx.spring_layout(chain0)
    #nx.draw(temp_chain0, pos, with_labels=True, node_size=50, font_size=8)
    #nodes = temp_chain0.nodes
    #sorted_nodes = sorted(nodes, key=natural_key)

#    print(neighbors)
    #plt.show()


#scd, edges = membrane.voronoi_all_lip_order(
#                                            lipid_list=lipid_list,
#                                            layer = layer,
#                                            start = start,
#                                            final=final,
#                                            )

#plt.imshow(scd, cmap='viridis',extent =edges)
#plt.colorbar(label='SCD Order Parameter')
#plt.show()

