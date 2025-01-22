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


class TestMathematicalModelsOld(unittest.TestCase):

    def test_mathematical_models(self):

        #read and plot specKS data
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Modeling'
        demographics_run_list=['DGKS_100_100_m5_RC7_m01d14y2025_h09m14s18',
                               'DGKS_1000_1000_m5_BI40_RC7_m01d13y2025_h15m36s22']
        #demographics_run_list=['TE07_fix__m01d08y2025_h15m09s22','TE08_fix_m01d14y2025_h09m17s20']
        #  ,'TE09_fix_m01d13y2025_h14m13s12']
        #'TE05fix__m01d06y2025_h11m17s15'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'
        specks_run_list = ["specks_TE100_m01d13y2025_h13m17s56","specks_TE1000_m01d13y2025_h13m17s53"]
        #specks_run_list = [ "specks_TE07_m01d14y2025_h09m50s16",
        #                       "specks_TE08_m01d14y2025_h09m50s16"]
        #                       "specks_TE09_m01d14y2025_h09m50s16"]
        #"specks_TE05_m01d14y2025_h09m50s16",


        bin_sizes_Tc = [200,200, 200, 200,200]
        #bin_sizes_Ks = [0.0002, 0.0002, 0.0002, 0.0002, 0.0002]
        num_bins_Ks=[200 for f in demographics_run_list]
        xmax_Ks = [False for f in demographics_run_list]#[0.025,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [False,False,False,False,False]
        run_list_name="Ks_Modeling_for_varying_varying_Tdiv"
        #<mutation_rate>0.000000012</mutation_rate>
        #<mutation_rate>0.000000012</mutation_rate>
        Ks_per_YR = 0.00000001
        #0.01 * 10**-6 #Mut rate is 1.2 to Ks rate of 1.0 in SpecKS
        ymax_KS = [False for f in demographics_run_list]

        #ks_prediction_Tnow, ks_prediction_Tdiv, ks_model_with_dispersion
        show_KS_predictions=[False,True,False]

        total_num_genes=3333
        suptitle = "SLiM and SpecKS Ks histograms\n" + \
                                  "Recombination rate = 8e-9, Ne and BI constant"
        make_Tc_Ks_MathModel_fig_with_subplots(num_bins_Ks, bin_sizes_Tc,
                                          demographiKS_out_path, demographics_run_list, run_list_name,
                                          specks_run_list, specks_out_path,Ks_per_YR,
                                     xmax_Ks, xmax_Tc, ymax_KS, suptitle,
                                     show_KS_predictions,total_num_genes)

        #test model
        self.assertEqual(True, True)  # add assertion here


def make_Tc_Ks_MathModel_fig_with_subplots(num_bins_Ks, bin_sizes_Tc,
                                 demographiKS_out_path, demographics_TE9_run_list, run_list_name,
                                 specks_TE9_run_list, specks_out_path, Ks_per_YR,
                                 xmax_Ks, xmax_Tc, ymax_Ks, suptitle, show_KS_predictions, total_num_genes):
    ymax_Tc = False
    num_runs = len(demographics_TE9_run_list)
    png_out = os.path.join(demographiKS_out_path, "ks_hist_by_TE{0}_test.png".format(run_list_name))
    fig, ax = plt.subplots(3, num_runs, figsize=(20, 15))
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

        plot_ks_model(ax[0, i], config_used, demographiKS_ks_results, [],
                Ks_per_YR,
                dgx_run_duration_in_m, spx_run_duration_in_m,
                plot_title, num_bins_Ks[i], xmax_Ks[i], ymax_Ks[i], show_KS_predictions)

        plot_ks_model(ax[1, i], config_used, [], spx_ks_results,
                Ks_per_YR,
                dgx_run_duration_in_m, spx_run_duration_in_m,
                plot_title, num_bins_Ks[i], xmax_Ks[i], ymax_Ks[i], show_KS_predictions)

        slim_csv_file = os.path.join(dgx_run_path, "simulated_ancestral_gene_mrcas.csv")
        loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)

        theory_output_file = os.path.join(dgx_run_path, "theoretical_ancestral_gene_mrcas.csv")
        loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
        burnin_times_in_generations=config_used.burnin_time
        plot_title = "Tcoal at Tdiv\nburnin time=" + str(burnin_times_in_generations) + " gen, " \
                     + "Ne=" + str(config_used.ancestral_Ne)

        #theory_mrcas_by_gene=False
        plot_mrca(ax[2, i], slim_mrcas_by_gene, specks_mrcas_by_gene, theory_mrcas_by_gene,
                  dgx_run_duration_in_m, plot_title, config_used.ancestral_Ne,
                  Ks_per_YR, bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc, total_num_genes)

    ax[0, 1].set(ylabel="# paralog pairs in bin")
    ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=550)
    plt.clf()
    plt.close()


def plot_ks_model(this_ax, config_used, slim_ks_by_gene, spx_ks_by_gene, Ks_per_YR,
            slim_run_duration_in_m, specks_run_duration_in_m, title, num_bins_Ks, xmax, ymax,
            show_KS_predictions):

    num_slim_genes = len(slim_ks_by_gene)
    num_specks_genes = len(spx_ks_by_gene)
    mean_ks_from_Tc = 2.0 * config_used.ancestral_Ne * Ks_per_YR
    t_div_as_ks= config_used.DIV_time_Ge * Ks_per_YR
    total_ks_shift=mean_ks_from_Tc+t_div_as_ks

    #if not xmax:
    #    if len(slim_ks_by_gene)>0:
    #        xmax = max(slim_ks_by_gene)
    #    else:
    #        xmax = max(spx_ks_by_gene)
    #bins = np.arange(0, xmax, bin_size)

    if len(slim_ks_by_gene) > 0:
        n, bins, patches = this_ax.hist(slim_ks_by_gene, bins=num_bins_Ks, facecolor='b', alpha=0.25,
                     density=False, label="DemographiKS")
        bin_midpoints=[0.5*(bins[i]+bins[i+1]) for i in range(0,len(bins)-1)]
        #gaussian_fit_curve_ys1, xs, popt = fit_curve_to_xs_and_ys(bin_midpoints, n, curve_fitting.wgd_gaussian)
        #if gaussian_fit_curve_ys1:
        #    this_ax.plot(xs, gaussian_fit_curve_ys1, c='k', label='DemographiKS curve fit\npopt: ' + str(popt))

    if len(spx_ks_by_gene) > 0:


        n, bins, patches = this_ax.hist(spx_ks_by_gene, bins=num_bins_Ks, facecolor='c', alpha=0.5,
                     density=False, label="SpecKS")
        bin_midpoints=[0.5*(bins[i]+bins[i+1]) for i in range(0,len(bins)-1)]
        #Amp, K, loc, scale):
        #p0 = [num_specks_genes, 2*config_used.bottleneck_Ne,
        #      total_ks_shift, config_used.DIV_time_Ge*0.1]
        #gaussian_fit_curve_ys1, xs, popt = fit_curve_to_xs_and_ys(bin_midpoints,n,
        #                                        curve_fitting.gaussian_modified_exponential,
        #                                                          p0=p0)
        #if gaussian_fit_curve_ys1:
        #    this_ax.plot(xs,gaussian_fit_curve_ys1,c='k',  label='SpecKS curve fit\npopt: ')#+ str(popt))

        two_Ne= 2.0*config_used.bottleneck_Ne
        bin_size=bins[1]-bins[0]
        Amp = bin_size * num_specks_genes / two_Ne
        K = two_Ne
        O = config_used.DIV_time_Ge * Ks_per_YR
        #ys=  [travelling_exp(x, Amp,K,O) for x in xs]
        #bin_size_in_time = bin_size / config_used.mutation_rate
        #kingman = [min(num_specks_genes,
        #                   (bin_size_in_time  * num_specks_genes / two_Ne) * math.e ** (
        #                               (-1 * (i / config_used.mutation_rate)) / two_Ne))
        #               for i in bins]
        #this_ax.plot(bins, kingman, c='r', label='Travelling Ks')  # + str(popt))

    #expected_Ks_peak_shift = config_used.DIV_time_Ge * config_used.mutation_rate
    #t_div_as_ks = config_used.DIV_time_Ge * Ks_per_YR
    this_ax.axvline(x=t_div_as_ks, color='b', linestyle='--', label="input Tdiv as Ks: " + str(t_div_as_ks) )
    this_ax.axvline(x=total_ks_shift, color='k', linestyle=':', label="Expected Ks mean: "
                                                                      + str(total_ks_shift))
    bin_size = bins[1] - bins[0]
    ks_prediction_Tnow, ks_prediction_Tdiv, ks_model_with_dispersion = (
        predict_Ks(config_used.bottleneck_Ne, mean_ks_from_Tc, Ks_per_YR, bin_size, bins, config_used, 3333))

    if show_KS_predictions[0]:
        this_ax.plot(bins,ks_prediction_Tdiv,c='red', label='Ks_exp at Tdiv (Kingman assumption)',alpha=1)
    if show_KS_predictions[1]:
        this_ax.plot(bins,ks_prediction_Tnow[0:len(bins)],
                 c='gray', label='Ks_exp at Tnow (raw)',alpha=1)
    if show_KS_predictions[2]:
        this_ax.plot(bins, ks_model_with_dispersion[0:len(bins)], c='k',
                 label='Ks_exp at Tnow (with dispersion)',alpha=1)

    x_axis_label = "Ks \n" + "SLiM run time: " + str(round(slim_run_duration_in_m, 2)) + " min\n" + \
                   "SpecKS run time: " + str(round(specks_run_duration_in_m,2)) + " min"
    if ymax:
        this_ax.set(ylim=[0, ymax])
    #this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel=x_axis_label)
    this_ax.set(title=title)
    # this_ax.set(title="Recom. rate: " + str(recombination_rate))
    this_ax.legend()
    # plt.ylabel("# genes in bin")
    # plt.savefig(png_out)


if __name__ == '__main__':
    unittest.main()
