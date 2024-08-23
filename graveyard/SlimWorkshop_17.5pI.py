# Keywords: Python, tree-sequence recording, tree sequence recording
import os
# This is a Python recipe; note that it runs the SLiM model internally, below

import subprocess, tskit
import matplotlib.pyplot as plt
import numpy as np

slim_folder = "/home/tamsen/Data/SLiM_scripts/"
slim_folder = "/home/tamsen/Data/SLiM_scripts/"
my_SLiM_script = os.path.join("../SLiM_scripts", "SlimWorkshop_17.5p1.slim")

# Run the SLiM model and load the resulting .trees file
#subprocess.check_output(["slim", "-m", "-s", "0", "./admix.slim"])
subprocess.check_output(["slim", "-m", "-s", "0", my_SLiM_script])
trees_file = os.path.join(slim_folder, "Tex1_TS.trees")
ts = tskit.load("./admix.trees")

# Load the .trees file and assess true local ancestry
breaks = np.zeros(ts.num_trees + 1)
ancestry = np.zeros(ts.num_trees + 1)
for tree in ts.trees():
    subpop_sum, subpop_weights = 0, 0
    for root in tree.roots:
        leaves_count = tree.num_samples(root) - 1  # subtract one for the root, which is a sample
        subpop_sum += tree.population(root) * leaves_count
        subpop_weights += leaves_count
    breaks[tree.index] = tree.interval[0]
    ancestry[tree.index] = subpop_sum / subpop_weights
breaks[-1] = ts.sequence_length
ancestry[-1] = ancestry[-2]

# Make a simple plot
plt.plot(breaks, ancestry)
plt.show()
