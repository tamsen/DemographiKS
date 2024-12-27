import os
import unittest

import numpy as np
from matplotlib import pyplot as plt

from data_aggregation.coalescent_plot_aggregation import get_run_time_in_minutes, read_data_csv, plot_mrca
from data_aggregation.histogram_plotter import read_Ks_csv,make_simple_histogram


class TestKsPlotAgg(unittest.TestCase):

    def test_Ks_for_TE_sims(self):



        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'

        #TE1_run_list = ['TE01_m12d20y2024_h14m16s21', 'TE02_m12d20y2024_h14m22s23',
        #               'TE03_m12d20y2024_h14m26s56']

        TE5_run_list=[
            'TE03_m12d20y2024_h14m26s56','TE05_m12d23y2024_h08m52s16','TE07_m12d23y2024_h09m18s26',
            'TE08_m12d24y2024_h09m31s26','TE09_m12d26y2024_h09m10s55']
        #'TE06_m12d23y2024_h09m10s28',
        #             'TE07_m12d23y2024_h09m18s26']

        Ne=[1000, 1000, 1000, 1000, 1000]
        burnin_times_in_generations=[5e7, 5e7, 5e7, 5e7, 5e7]
        time_since_DIV=[1000, 10000,100000,500000,1000000]
        bin_sizes_Tc = [200,200, 200, 200,200]
        #bin_sizes_Ks = [0.00005,0.0001,0.0004,0.0002,0.0001]
        bin_sizes_Ks = [0.0002, 0.0002, 0.0002, 0.0002, 0.0002]
        run_list_num="_5to9"
        num_runs=len(TE5_run_list)
        png_out = os.path.join(demographiKS_out_path,"ks_hist_by_TE{0}_test.png".format(run_list_num))

        fig, ax = plt.subplots(2, num_runs, figsize=(20, 10))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
                     "Recombination rate = 8e-9, Ne and BI constant")
        #    "Recombination rate = 8e-9, Ne varies, BI varies")
        xmax_Ks = 0.025#0.001  # max(demographiKS_ks_results)
        xmax_Tc = False
        #ymax=100
        ymax=False



        for i in range(0,num_runs):
            run_name=TE5_run_list[i]
            run_name_splat=run_name.split("_")
            nickname="_".join(run_name_splat[0:3])
            run_path=os.path.join(demographiKS_out_path,run_name)
            csv_file_name='allotetraploid_bottleneck.csv'
            demographiKS_ks_results = read_Ks_csv(os.path.join(run_path,csv_file_name))
            plot_title="burnin time=" + str(burnin_times_in_generations[i]) + " gen,\n"\
                        + "Ne=" + str(Ne[i]) + ", Tdiv=" + str(time_since_DIV[i])
            run_duration_in_m = get_run_time_in_minutes(run_path)
            plot_ks(ax[0,i], demographiKS_ks_results,time_since_DIV[i],
                    run_duration_in_m, plot_title, bin_sizes_Ks[i], xmax_Ks, ymax)

            slim_csv_file = os.path.join(run_path, "simulated_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)

            theory_output_file = os.path.join(run_path, "theoretical_ancestral_gene_mrcas.csv")
            loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
            plot_title = "burnin time=" + str(burnin_times_in_generations[i]) + " gen\n" \
                         + "Ne=" + str(Ne[i])
            plot_mrca(ax[1,i], slim_mrcas_by_gene, [], theory_mrcas_by_gene,
                      run_duration_in_m, plot_title, bin_sizes_Tc[i], xmax_Tc, ymax)

        ax[0,0].set(ylabel="# genes in bin")
        ax[1,0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here


def plot_ks(this_ax, slim_ks_by_gene, t_div,
            run_duration_in_m, title, bin_size, xmax, ymax):

    num_genes = len(slim_ks_by_gene)

    if not xmax:
        xmax = max(slim_ks_by_gene)
    bins = np.arange(0, xmax, bin_size)

    if len(slim_ks_by_gene) > 0:
        this_ax.hist(slim_ks_by_gene, bins=bins, facecolor='b', alpha=0.25,
                     label='SLiM Ks by gene\n'
                           + "(" + str(num_genes) + " genes in genome)",
                     density=False)

    if t_div:
        t_div_as_ks= t_div * 0.01 * 10**-6
        this_ax.axvline(x=t_div_as_ks, color='b', linestyle='--', label="input Tdiv as Ks")

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
