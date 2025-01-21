import unittest
import glob
import math
import os
import unittest
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import config
from data_aggregation import curve_fitting
from data_aggregation.coalescent_plot_aggregation import plot_mrca, read_data_csv, get_run_time_in_minutes
from data_aggregation.curve_fitting import fit_curve_to_xs_and_ys, travelling_exp
from data_aggregation.histogram_plotter import read_Ks_csv
from data_aggregation.ks_plot_aggregations import plot_ks, plot_expository_images, predict_Ks

class MyTestCase2(unittest.TestCase):
    def test_modeling(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Modeling2'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Modeling2'

        demographics_TE_run_list = ["DGKS_10_10_m5_RC7_m01d14y2025_h09m14s47",
                                    "DGKS_100_100_m5_RC7_m01d14y2025_h09m14s18",
                                    "DGKS_1000_1000_m5_BI40_RC7_m01d13y2025_h15m36s22",
                                    "DGKS_5000_5000_m5_BI_40K_RC7_m01d14y2025_h09m13s52"]

        specks_TE_run_list = ['specks_TE10_m01d13y2025_h13m18s28',
                              'specks_TE100_m01d13y2025_h13m17s56',
                              'specks_TE1000_m01d13y2025_h13m17s53',
                              'specks_TE5000_m01d13y2025_h13m18s40']

        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate
        Ks_per_YR = 0.833 * 10 ** -5
        bin_sizes_Tc = [10, 40, 400, 800]  # looks good
        xmax_Ks = [0.04, 0.04, 0.1, 0.4]  # [0.01,0.01,0.01,0.1,0.2]#False#0.08  # for mut rate e-5
        bin_sizes_Ks = [0.001, 0.001, 0.002, 0.008]
        xmax_Tc = [400, 2000, 20000, 40000]
        run_list_num = "_modeling2"
        ymax_Ks = [False, False, False, False, False]

        suptitle = "SLiM vs SpecKS, Tcoal and Ks\n" + \
                   "Recombination rate = 1.26e-7, mut rate 1.0e-5"

        #kingman, Ks_model_exponential, ks_model_smoothed_exponential, ks_model_as_gaussian
        #ks_prediction_Tnow, ks_prediction_Tdiv, ks_model_with_dispersion
        show_KS_predictions = [False, True, False,True]
        total_num_genes = 333
        make_Tc_Ks_model_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                     demographiKS_out_path, demographics_TE_run_list, run_list_num,
                                     specks_TE_run_list, specks_out_path, Ks_per_YR,
                                     xmax_Ks, xmax_Tc, ymax_Ks, suptitle, show_KS_predictions, total_num_genes)

        self.assertEqual(True, True)  # add assertion here

def make_Tc_Ks_model_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                 demographiKS_out_path, demographics_TE9_run_list, run_list_name,
                                 specks_TE9_run_list, specks_out_path, Ks_per_YR,
                                 xmax_Ks, xmax_Tc, ymax_Ks, suptitle, show_KS_predictions, total_num_genes):
    ymax_Tc = False
    num_runs = len(demographics_TE9_run_list)
    png_out = os.path.join(demographiKS_out_path, "ks_hist_by_TE{0}_test.png".format(run_list_name))
    fig, ax = plt.subplots(2, num_runs, figsize=(20, 10))
    fig.suptitle(suptitle)
    for i in range(0, num_runs):
        dgx_run_name = demographics_TE9_run_list[i]

        if dgx_run_name:

            dgx_run_path = os.path.join(demographiKS_out_path, dgx_run_name)
            print("dgx_run_path: " +dgx_run_path )
            glob_results=glob.glob(dgx_run_path + '/*.used.xml')
            input_xml_file = glob_results[0]
            config_used = config.DemographiKS_config(input_xml_file)
            csv_file_name = 'allotetraploid_bottleneck.csv'
            ks_file = os.path.join(dgx_run_path, csv_file_name)
            print("reading " + ks_file)
            demographiKS_ks_results = read_Ks_csv(ks_file, False)
            dgx_run_duration_in_m = get_run_time_in_minutes(dgx_run_path)
            plot_title = "burnin time=" + str(config_used.burnin_time) + " gen,\n" \
                     + "Ne=" + str(config_used.ancestral_Ne) +\
                         ", Tdiv=" + str(config_used.DIV_time_Ge)
        else:
            config_used = False
            demographiKS_ks_results = []
            dgx_run_duration_in_m = 0
            plot_title = "foo - didnt load a config"


        spx_run_name = specks_TE9_run_list[i]
        if spx_run_name:
            spx_run_nickname = spx_run_name.split('_')[1]
            spx_run_path = os.path.join(specks_out_path, spx_run_name)
            csv_file_name = 'Allo_' + spx_run_nickname + '_ML_rep0_Ks_by_GeneTree.csv'
            spx_ks_results = read_Ks_csv(os.path.join(spx_run_path,csv_file_name), True)
            spx_run_duration_in_m = get_run_time_in_minutes(spx_run_path)
            specks_csv_file = os.path.join(spx_run_path, "variations_in_div_time.txt")
            loci, specks_mrcas_by_gene = read_data_csv(specks_csv_file)

        else:
            spx_ks_results = []
            spx_run_duration_in_m = 0
            specks_mrcas_by_gene = False

        plot_ks_for_modelling(ax[0, i], config_used, [], spx_ks_results, config_used.DIV_time_Ge,
                config_used.ancestral_Ne, Ks_per_YR,
                dgx_run_duration_in_m, spx_run_duration_in_m,
                plot_title, bin_sizes_Ks[i], xmax_Ks[i], ymax_Ks[i], show_KS_predictions)

        plot_ks_for_modelling(ax[1, i], config_used, demographiKS_ks_results, [], config_used.DIV_time_Ge,
                config_used.ancestral_Ne, Ks_per_YR,
                dgx_run_duration_in_m, spx_run_duration_in_m,
                plot_title, bin_sizes_Ks[i], xmax_Ks[i], ymax_Ks[i], show_KS_predictions)

        slim_csv_file = os.path.join(dgx_run_path, "simulated_ancestral_gene_mrcas.csv")
        loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)

        theory_output_file = os.path.join(dgx_run_path, "theoretical_ancestral_gene_mrcas.csv")
        loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
        burnin_times_in_generations=config_used.burnin_time
        plot_title = "Tcoal at Tdiv\nburnin time=" + str(burnin_times_in_generations) + " gen, " \
                     + "Ne=" + str(config_used.ancestral_Ne)

        #theory_mrcas_by_gene=False
        #plot_mrca(ax[2, i], slim_mrcas_by_gene, specks_mrcas_by_gene, theory_mrcas_by_gene,
        #          dgx_run_duration_in_m, plot_title, config_used.ancestral_Ne,
        #          Ks_per_YR, bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc, total_num_genes)

    ax[0, 1].set(ylabel="# paralog pairs in bin")
    ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=550)
    plt.clf()
    plt.close()

def plot_ks_for_modelling(this_ax, config_used, slim_ks_by_gene, spx_ks_by_gene, t_div,Ne, Ks_per_YR,
            slim_run_duration_in_m, specks_run_duration_in_m, title, bin_size, xmax, ymax,
            show_KS_predictions):

    num_slim_genes = len(slim_ks_by_gene)
    num_specks_genes = len(spx_ks_by_gene)
    mean_ks_from_Tc = 2.0 * config_used.ancestral_Ne * Ks_per_YR

    if not xmax:
        xmax = max(slim_ks_by_gene)
    bins = np.arange(0, xmax, bin_size)

    if len(slim_ks_by_gene) > 0:
        this_ax.hist(slim_ks_by_gene, bins=bins, facecolor='b', alpha=0.25,
                     label='SLiM Ks by gene\n'
                           + "(" + str(num_slim_genes) + " paralogs in genome)",
                     density=False)

    if len(spx_ks_by_gene) > 0:
        this_ax.hist(spx_ks_by_gene, bins=bins, facecolor='c', alpha=0.25,
                     label='SpecKS Ks by gene\n'
                           + "(" + str(num_specks_genes) + " paralogs in genome)",
                     density=False)

    if t_div:
        t_div_as_ks= config_used.DIV_time_Ge * Ks_per_YR
        this_ax.axvline(x=t_div_as_ks, color='b', linestyle='--', label="input Tdiv as Ks")
        total_ks_shift=mean_ks_from_Tc+t_div_as_ks
        this_ax.axvline(x=total_ks_shift, color='k', linestyle='--', label="Expected Ks mean")

    expected_num_genes=int(config_used.total_num_bases / config_used.gene_length)

    ks_prediction_Tdiv, ks_model_expontial_Tnow, ks_model_smoothed_Tnow, ks_model_gaussian = (
        predict_Ks(Ne, mean_ks_from_Tc, t_div_as_ks, bin_size, bins, config_used, expected_num_genes))

    if show_KS_predictions[0]:
        this_ax.plot(bins,ks_prediction_Tdiv,c='red', label='Ks_exp at Tdiv (Kingman assumption)',alpha=1)
    if show_KS_predictions[1]:
        this_ax.plot(bins,ks_model_expontial_Tnow[0:len(bins)],
                 c='gray', label='Ks_exp at Tnow (Kingman)',alpha=1,linestyle='--',)
    if show_KS_predictions[2]:
        this_ax.plot(bins, ks_model_smoothed_Tnow[0:len(bins)], c='r',
                 label='Ks_exp at Tnow (Smoothed Kingman)',alpha=1,linestyle=':',)
    if show_KS_predictions[3]:
        this_ax.plot(bins, ks_model_gaussian[0:len(bins)], c='gray',
                 label='Ks_exp at Tnow (Gaussian)',alpha=1,linestyle=':',)

    x_axis_label = "Ks \n" + "SLiM run time: " + str(round(slim_run_duration_in_m, 2)) + " min\n" + \
                   "SpecKS run time: " + str(round(specks_run_duration_in_m,2)) + " min"
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
