import os
import unittest

import numpy as np
from matplotlib import pyplot as plt


class TestTCoalAggregation(unittest.TestCase):

    def test_TCoal_with_varying_RC(self):

        #make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx"

        RC_run_list=["RC09_m12d18y2024_h14m30s54","RC10_m12d18y2024_h13m36s12",
                     "RC11_m12d18y2024_h14m35s15","RC12_m12d18y2024_h14m35s17",
                     "RC13_m12d18y2024_h14m35s19","RC14_m12d18y2024_h14m39s19"]
        recombination_rates=[8e-9,8e-10,8e-11,8e-12,8e-13,8e-14]#8e-8,
        theory_output_file=os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                        RC_run_list[0],"theoretical_ancestral_gene_mrcas.csv")
        num_runs=len(RC_run_list)
        bin_size = 10000
        xmax = 160_000
        png_out=os.path.join(aggragate_output_folder,"mrcsa.png")
        fig, ax = plt.subplots(1, num_runs,figsize=(20,5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
            "burnin_time = 5000000 generations")

        for i in range(0,num_runs):
            run_name=RC_run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)
            slim_csv_file=os.path.join(local_output_folder,"simulated_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene=read_data_csv(slim_csv_file)
            loci, theory_mrcas_by_gene=read_data_csv(theory_output_file)

            plot_mrca(ax[i],slim_mrcas_by_gene, [],theory_mrcas_by_gene,
                recombination_rates[i],bin_size,xmax)
        
        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()
        
        self.assertEqual(True, True)  # add assertion here


def read_data_csv(csv_file):

    loci=[]
    mrcas=[]

    with open(csv_file, "r") as f:

        while True:
            line = f.readline()
            if "Start" in line:
                continue
            if len(line)==0:
                break

            data = line.strip().split(",")
            loci.append(int(data[0]))
            mrcas.append(float(data[1]))

    return loci,mrcas


def plot_mrca(this_ax,slim_mrcas_by_gene, slim_mrcas_by_tree, theoretical_mrcas_by_gene,
              recombination_rate, bin_size,xmax):

    #fig = plt.figure(figsize=(10, 10), dpi=350
    #Co.T=(1/2N)*e^-((t-1)/2N))

    max_mrca = max(slim_mrcas_by_gene)
    num_genes=len(slim_mrcas_by_gene)
    num_segments=len(slim_mrcas_by_tree)
    bins = np.arange(0, xmax , bin_size)

    if len(slim_mrcas_by_gene) > 0:
        this_ax.hist(slim_mrcas_by_gene, bins=bins, facecolor='b', alpha=0.25,
                                label='SLiM Tcoal by gene',
                                density=False)
    #label = 'SLiM Tcoal by gene (total: ' + str(num_genes) + ')',

    if len(slim_mrcas_by_tree) > 0:
        this_ax.hist(slim_mrcas_by_tree, bins=bins, facecolor='r', alpha=0.25,
                 label='SLiM Tcoal by segment (total: '+str(num_segments)+')',
                 density=False)
    
    if len(theoretical_mrcas_by_gene) > 0:
        this_ax.hist(theoretical_mrcas_by_gene, bins=bins, facecolor='c', alpha=0.25,
                 label='Theoretical Tcoal by gene',
                 density=False)

    this_ax.set(ylim=[0, 2_000])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(title="Recom. rate: " + str(recombination_rate))
    # plt.xlim([0, max_mrca])
    # plt.ylim([0, 300])
    this_ax.legend()
    #plt.xlabel("MRCA time")
    #plt.ylabel("# genes in bin")
    #plt.savefig(png_out)

if __name__ == '__main__':
    unittest.main()
