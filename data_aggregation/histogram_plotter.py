import os
import unittest

import numpy as np
from matplotlib import pyplot as plt


class Test_Plot_Histogram(unittest.TestCase):
    def test_plot_histogram(self):
        specks_full_path='foo'
        specks_ks_results = read_Ks_csv(specks_full_path)
        out_png1 = os.path.join("hist_comparison_", "specks_DemographiKS_out.png")

        specks_hist_data = make_simple_histogram(specks_ks_results, species_run_name, bin_size, color, WGD_ks,
                                                 max_Ks, density, out_png1)

        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()

def make_simple_histogram(Ks_results, species_name, bin_size, color,WGD_ks, max_Ks, density, out_png):

    # MBE says: 600 - 1200 dpi for line drawings
    # and 350 dpi for color and half-tone artwork)
    fig = plt.figure(figsize=(10, 10), dpi=350)
    x = Ks_results
    # print(PAML_hist_out_file)
    label="hist for " + os.path.basename(out_png).replace("_out.png","")
    if max_Ks:
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25,
                                    label=label, density=density)
        plt.xlim([0, max_Ks * (1.1)])


    plt.axvline(x=WGD_ks, color='b', linestyle='-', label="WGD paralog start")
    num_pairs=sum(n)
    num_after_wgd=sum([n[i] for i in range(0,len(n)) if bins[i] > WGD_ks])
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for {0}.\n{1} pairs of genes. ~{2} retained from WGD.".format(
        species_name,num_pairs,num_after_wgd))
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return [n,bins]

def read_Ks_csv(csv_file):

    ks_results = []
    with open(csv_file, "r") as f:

        reading_header=True
        while True:
            line = f.readline()
            if "ersion" in line:
                continue
            if "Git" in line:
                continue
            if "leaf names" in line:
                continue
            if not line:
                break
            if len(line)==0:
                break
            if reading_header:
                reading_header=False
                continue
            data = line.split(",")
            #print(data)
            ks_value=float(data[2])
            ks_results.append(ks_value)

    return ks_results
