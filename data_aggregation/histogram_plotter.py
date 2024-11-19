import os
import unittest

import numpy as np
from matplotlib import pyplot as plt


class Test_Plot_Histogram(unittest.TestCase):
    def test_plot_histogram(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        run_name='DGKS_1000_1000_m08d29y2024_h16m08s09'
        #DGKS_1000_1000_m08d29y2024_h16m08s09
        #DGKS_100_100_m08d29y2024_h16m09s38
        #DGKS_1000_100_m08d29y2024_h15m27s59
        #DGKS_10_10_m08d29y2024_h16m13s27
        #DGKS_1000_10_m08d29y2024_h16m12s08
        run_name_splat=run_name.split("_")
        nickname="_".join(run_name_splat[0:3])
        run_path=os.path.join(demographiKS_out_path,run_name)
        csv_file_name='allotetraploid_bottleneck.csv'
        demographiKS_ks_results = read_Ks_csv(os.path.join(run_path,csv_file_name))
        out_png1 = os.path.join(run_path,"DemographiKS_out.png")
        WGD_time_in_Ks=25*0.01*10**-6
        DIV_time_in_Ks=75*0.01*10**-6
        max_Ks = 0.1
        bin_size = 0.001
        make_simple_histogram(demographiKS_ks_results,nickname, bin_size,
        'k', WGD_time_in_Ks,DIV_time_in_Ks,
                                                 max_Ks, 0.5, out_png1)
        self.assertEqual(os.path.exists(out_png1), True)  # add assertion here


if __name__ == '__main__':
    unittest.main()

def make_simple_histogram(Ks_results, title, bin_size, color,WGD_ks,
                          DIV_ks, max_Ks, density, out_png):

    # MBE says: 600 - 1200 dpi for line drawings
    # and 350 dpi for color and half-tone artwork)
    fig = plt.figure(figsize=(10, 10), dpi=350)
    x = Ks_results
    # print(PAML_hist_out_file)
    label="hist for " + os.path.basename(out_png).replace("_out.png","")

    plt.axvline(x=WGD_ks, color='b', linestyle='-', label="WGD (" +str(WGD_ks)+ ")")
    plt.axvline(x=DIV_ks, color='r', linestyle='-', label="DIV (" +str(DIV_ks) + ")")

    if max_Ks:
        bins = np.arange(0, max_Ks + 0.1, bin_size)
        n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25,
                                    label=label, density=density)
        plt.xlim([0, max_Ks * (1.1)])

    #n, bins, patches = plt.hist(x, bins=100,facecolor=color, alpha=0.25,
    #                            label=label, density=density)
    plt.title(title)
    plt.xlim([0, max_Ks])
    plt.ylim([0, 700])
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")

    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return

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
