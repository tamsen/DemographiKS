import os
import shutil
import sys
from datetime import datetime#
#import subprocess, msprime, pyslim
import subprocess, pyslim

import msprime
import numpy as np
import tskit


def run_sim():

    #conf = setup(sys.argv)
    #if not conf:
    #    return

    print("hello!!")
    # Run the SLiM model

    slim_folder="/home/tamsen/Data/SLiM_scripts/"
    my_SLiM_script= os.path.join("SLiM_scripts", "TreeSeqRecording.slim")
    trees_file=os.path.join(slim_folder,"Tex1_TS.trees")
    sequence_file=os.path.join(slim_folder,"ex1_TS_overlaid.trees")
    subprocess.check_output(["slim", "-m", "-s", "0", my_SLiM_script])
    # Overlay neutral mutations
    ts = tskit.load(trees_file)

    #https://tskit.dev/msprime/docs/stable/mutations.html
    mts = msprime.sim_mutations(ts, rate=1e-11, random_seed=42, keep=True)
    mts.dump(sequence_file)

    i=0

    samp_ids = mts.samples()
    print("  ID of diploid individual: ", " ".join([f"{ts.node(s).individual:3}" for s in samp_ids]))
    print("       ID of (sample) node: ", " ".join([f"{s:3}" for s in samp_ids]))
    for var in mts.variants():
        #print(str(i))
        site=var.site
        alleles = np.array(var.alleles)
        print(f"Site {site.id} (ancestral state '{site.ancestral_state}')", alleles[var.genotypes])
        i=i+1
        #if i>5:
        #    return
    return

    print("Sequence file written here: " + sequence_file)


if __name__ == '__main__':
    run_sim()
