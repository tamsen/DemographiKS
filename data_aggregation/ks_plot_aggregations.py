import os
import unittest

import numpy as np
from matplotlib import pyplot as plt

from data_aggregation.coalescent_plot_aggregation import get_run_time_in_minutes
from data_aggregation.histogram_plotter import read_Ks_csv,make_simple_histogram


class TestKsPlotAgg(unittest.TestCase):

    def test_Ks_for_TE_sims(self):



        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'

        TE_run_list = ['TE01_m12d20y2024_h14m16s21', 'TE02_m12d20y2024_h14m22s23',
                       'TE03_m12d20y2024_h14m26s56'] 

        num_runs=len(TE_run_list)
        png_out = os.path.join(demographiKS_out_path,"ks_hist_by_TE_test.png")

        fig, ax = plt.subplots(1, num_runs, figsize=(20, 5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
            "Recombination rate = 8e-9, Ne varies, BI varies")
        max_Ks = 0.001  # max(demographiKS_ks_results)
        ymax=100
        bin_size = 0.00001  # max_Ks/10.0
        Ne=[10, 100, 1000]
        burnin_times_in_generations=[5e6, 5e7, 5e7]
        for i in range(0,num_runs):
            run_name=TE_run_list[i]
            run_name_splat=run_name.split("_")
            nickname="_".join(run_name_splat[0:3])
            run_path=os.path.join(demographiKS_out_path,run_name)
            csv_file_name='allotetraploid_bottleneck.csv'
            demographiKS_ks_results = read_Ks_csv(os.path.join(run_path,csv_file_name))
            print("run path: " + run_path)
            #WGD_time_in_Ks=25*0.01*10**-6
            #DIV_time_in_Ks=75*0.01*10**-6
            WGD_time_in_Ks=False#20*0.01
            DIV_time_in_Ks=False#20*0.01
            plot_title="burnin time=" + str(burnin_times_in_generations[i]) + " gen\n"\
                        + "Ne=" + str(Ne[i])
            run_duration_in_m = get_run_time_in_minutes(run_path)
            plot_ks(ax[i], demographiKS_ks_results,
                    run_duration_in_m, plot_title, bin_size, max_Ks, ymax)

        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here


def plot_ks(this_ax, slim_ks_by_gene,
            run_duration_in_m, title, bin_size, xmax, ymax):

    num_genes = len(slim_ks_by_gene)

    if not xmax:
        xmax = max(slim_ks_by_gene)
    bins = np.arange(0, xmax, bin_size)

    if len(slim_ks_by_gene) > 0:
        this_ax.hist(slim_ks_by_gene, bins=bins, facecolor='b', alpha=1.0,
                     label='SLiM Ks by gene\n'
                           + "(" + str(num_genes) + " genes in genome)",
                     density=False)

    x_axis_label = "Ks \n" + "run time: " + str(round(run_duration_in_m, 2)) + " min"

    if ymax:
        this_ax.set(ylim=[0, ymax])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel=x_axis_label)
    this_ax.set(title=title)
    # this_ax.set(title="Recom. rate: " + str(recombination_rate))
    this_ax.legend()
    # plt.ylabel("# genes in bin")
    # plt.savefig(png_out)



if __name__ == '__main__':
    unittest.main()
