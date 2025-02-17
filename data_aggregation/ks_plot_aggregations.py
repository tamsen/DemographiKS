import glob
import math
import os
import unittest
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from scipy.optimize import curve_fit

import config
import ks_modeling
from data_aggregation import curve_fitting
from data_aggregation.coalescent_plot_aggregation import get_run_time_in_minutes, read_data_csv, plot_mrca, \
    add_mrca_annotations
from data_aggregation.curve_fitting import gaussian_modified_exponential
from data_aggregation.histogram_plotter import read_Ks_csv,make_simple_histogram


class TestKsPlotAgg(unittest.TestCase):

    def test_Ks_for_varying_Ne_early_runs(self):

        print('foo')

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'

        #<mutation_rate>1.2e-8</mutation_rate>,<WGD_time_Ge>1000</WGD_time_Ge>
        demographics_TE_run_list=[False, "DGKS_10_10_v2_m01d06y2025_h15m35s38",
                        "DGKS_100_100_v2_m01d06y2025_h15m35s43",
                            "DGKS_1000_1000_v2_m01d06y2025_h15m35s46"]


        #<mutation_rate>1.0e-5</mutation_rate>, <DIV_time_Ge>75</DIV_time_Ge>
        demographics_TE_run_list=[False, "DGKS_10_10_v2_m01d06y2025_h13m09s31",
                        "DGKS_100_100_v2_m01d06y2025_h13m09s35",
                            "DGKS_1000_1000_v2_m01d06y2025_h13m05s18"]

        specks_TE_run_list=[False,False,False,False]

        #Ks_per_YR = 0.01 * 10**-6
        Ks_per_YR = 10 ** -5
        Ne = [10,10, 100, 1000]
        #Ne=[500, 500, 1000]
        burnin_times_in_generations=[2e4,2e4, 2e4,2e4, 2e4]
        #time_since_DIV=[1000,1000,1000,1000]
        time_since_DIV = [75, 75, 75, 75]

        bin_sizes_Tc = [80, 80, 80, 80, 80]#looks good

        xmax_Ks = 0.05 #for mut rate e-5
        bin_sizes_Ks = [0.001, 0.001,0.001, 0.001, 0.001]


        #xmax_Ks = False # 0.00001  #for mut rate 1.2e-8
        #bin_sizes_Ks = [0.000001, 0.000001, 0.000001, 0.000001, 0.000001]
        xmax_Tc = 5000
        run_list_num = "_early_DGKS_by_Ne"
        ymax = False

        suptitle = "SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
                                  "Recombination rate = 1.26e-6, Ne and BI constant"

        make_Tc_Ks_fig_with_subplots(Ne, bin_sizes_Ks, bin_sizes_Tc, burnin_times_in_generations,
                                          demographiKS_out_path, demographics_TE_run_list, run_list_num,
                                          specks_TE_run_list, specks_out_path, time_since_DIV,Ks_per_YR,
                                          xmax_Ks, xmax_Tc, ymax,suptitle)

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_for_varying_Ne(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'


        demographics_TE_run_list=['DGKS_10_10_v2_m01d06y2025_h13m09s31','DGKS_10_10_v2_m01d06y2025_h13m09s31',
                                   'DGKS_100_100_v2_m01d06y2025_h13m09s35',
                                  'DGKS_1000_1000_v2_m01d06y2025_h13m05s18']

        specks_TE_run_list=['specks_TE07_m12d30y2024_h12m10s15','specks_TE07_m12d30y2024_h12m10s15',
                             'specks_TE07_m12d30y2024_h12m10s15']


        specks_TE_run_list=[False,False,False,False]


        Ne = [10,10, 100, 1000]
        #Ne=[500, 500, 1000]
        burnin_times_in_generations=[5e7, 5e7, 5e7, 5e7, 5e7]
        time_since_DIV=[25,100000, 100000,100000]

        bin_sizes_Tc = [40, 40, 40, 40, 40]#looks good
        bin_sizes_Ks = [0.005, 0.005,0.005, 0.005, 0.005]
        xmax_Ks = 0.1#0.001  # max(demographiKS_ks_results)
        xmax_Tc = 10000
        run_list_num = "_9to11_by_Ne"
        # end
        ymax = False

        make_Tc_Ks_fig_with_subplots(Ne, bin_sizes_Ks, bin_sizes_Tc, burnin_times_in_generations,
                                          demographiKS_out_path, demographics_TE_run_list, run_list_num,
                                          specks_TE_run_list, specks_out_path, time_since_DIV, xmax_Ks, xmax_Tc, ymax)

        self.assertEqual(True, True)  # add assertion here


    def test_show_Ks_for_varying_Tdiv_times(self):



        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'

        #TE 5 to 9 parameters - varies across time since DIV
        #start

        demographics_TE5_run_list=['TE15_m01d03y2025_h10m23s32',
                                   'TE15_m01d03y2025_h10m23s32','TE07_fix__m01d05y2025_h09m07s41']
        demographics_TE5_run_list=['TE03_m12d20y2024_h14m26s56','TE05fix__m01d03y2025_h11m36s57','TE07_fix__m01d05y2025_h09m07s41',
           'TE08_m12d24y2024_h09m31s26','TE09_m12d26y2024_h09m10s55']

        #specks_TE5_run_list=['specks_TE07_m12d30y2024_h12m10s15','specks_TE07_m12d30y2024_h12m10s15',
        #                     'specks_TE07_m12d30y2024_h12m10s15']
        specks_TE5_run_list=['specks_TE05_m12d30y2024_h11m50s03','specks_TE05_m12d30y2024_h11m50s03',
                             'specks_TE07_m12d30y2024_h12m10s15',
                             'specks_TE08_m12d30y2024_h12m10s13','specks_TE09_m12d30y2024_h12m10s11']


        #specks_TE5_run_list=['specks_TE05_m12d31y2024_h09m10s39','specks_TE05_m12d31y2024_h09m10s39',
        #                     'specks_TE07_m12d31y2024_h09m10s28',
        #                    'specks_TE08_m12d31y2024_h09m10s32',
        #                    'specks_TE09_m12d31y2024_h09m10s34']

        Ne = [1000, 1000, 1000, 1000, 1000]
        burnin_times_in_generations = [5e7, 5e7, 5e7, 5e7, 5e7]
        time_since_DIV = [1000, 10000, 100000, 500000, 1000000]

        bin_sizes_Tc = [200,200, 200, 200,200]
        bin_sizes_Ks = [0.0002, 0.0002, 0.0002, 0.0002, 0.0002]
        xmax_Ks = 0.025#0.001  # max(demographiKS_ks_results)
        xmax_Tc = False
        run_list_num="_5to9_vary_Tdiv_fix2"
        Ks_per_YR = 0.01 * 10**-6
        #end


        ymax = False

        suptitle = "SLiM Tcoal and Ks\n" + \
                                  "Recombination rate = 1.26e-6, Ne and BI constant"
        make_Tc_Ks_fig_with_subplots(Ne, bin_sizes_Ks, bin_sizes_Tc, burnin_times_in_generations,
                                          demographiKS_out_path, demographics_TE5_run_list, run_list_num,
                                          specks_TE5_run_list, specks_out_path, time_since_DIV,
                                            Ks_per_YR, xmax_Ks, xmax_Tc, ymax, suptitle)



        self.assertEqual(True, True)  # add assertion here


def make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                 demographiKS_out_path, demographics_TE9_run_list, run_list_name,
                                 specks_TE9_run_list, specks_out_path,
                                 xmax_Ks, xmax_Tc, ymax_Ks, ymax_Tc,
                                 suptitle, show_KS_predictions):

    num_runs = len(demographics_TE9_run_list)
    png_out = os.path.join(demographiKS_out_path, run_list_name)
    par_dir = Path(__file__).parent.parent
    image_folder = os.path.join(par_dir, "images")
    png_Tnow = os.path.join(image_folder, 'Ks_now_time_slice.jpg')
    png_Tdiv = os.path.join(image_folder, 'Tdiv_TimeSlice.jpg')
    png_with_migration_Tnow = os.path.join(image_folder, 'Migration_now.png')
    png_with_migration_Tdiv = os.path.join(image_folder, 'Migration_Tc.png')

    fig, ax = plt.subplots(2, num_runs, figsize=(40, 20))
    fig.suptitle(suptitle)

    for i in range(1, num_runs):
        dgx_run_name = demographics_TE9_run_list[i]

        if dgx_run_name:

            dgx_run_path = os.path.join(demographiKS_out_path, dgx_run_name)
            print("dgx_run_path: " +dgx_run_path )
            glob_results=glob.glob(dgx_run_path + '/*.used.xml')
            input_xml_file = glob_results[0]
            config_used = config.DemographiKS_config(input_xml_file)
            csv_file_name = config_used.sim_name + '.csv'
            ks_file = os.path.join(dgx_run_path, csv_file_name)
            print("reading " + ks_file)
            demographiKS_ks_results = read_Ks_csv(ks_file, False)
            dgx_run_duration_in_m = get_run_time_in_minutes(dgx_run_path)
            plot_title = "Ks at Tnow\n" + \
                         "burnin time=" + str(config_used.burnin_time) + " gen, " \
                     + "Na=" + str(config_used.ancestral_Ne)  + ", Nb=" + str(config_used.bottleneck_Ne) +\
                         ",\nTdiv=" + str(config_used.DIV_time_Ge) + ", RC=" + str(config_used.recombination_rate)
            if config_used.mig_rate:
                plot_title = plot_title + ", MigStart=" + str(config_used.mig_start) + \
                              ", MigEnd=" + str(config_used.mig_stop) +  ", MigRate=" + str(config_used.mig_rate)
            else:
                plot_title = plot_title + ", MigRate=0"

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


        plot_ks(ax[0, i], config_used, demographiKS_ks_results, spx_ks_results,
                plot_title, bin_sizes_Ks[i], xmax_Ks[i], ymax_Ks[i], show_KS_predictions)

        slim_csv_file = os.path.join(dgx_run_path, "simulated_ancestral_gene_mrcas.csv")
        loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)

        theory_output_file = os.path.join(dgx_run_path, "theoretical_ancestral_gene_mrcas.csv")
        loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
        plot_title = "Tcoal at Tdiv\nburnin time=" + str(config_used.burnin_time) + " gen, " \
                     + "Na=" + str(config_used.ancestral_Ne)

        theory_mrcas_by_gene=False
        avg_slim_Tc = plot_mrca(ax[1, i], slim_mrcas_by_gene, specks_mrcas_by_gene, theory_mrcas_by_gene,
                  plot_title, config_used.ancestral_Ne,
                  bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc[i], config_used.num_genes)

        add_mrca_annotations(ax[1, i], config_used, avg_slim_Tc, dgx_run_duration_in_m, spx_run_duration_in_m)

    if config_used.mig_rate > 0:
        plot_expository_images(ax, png_with_migration_Tdiv , png_with_migration_Tnow)
    else:
        plot_expository_images(ax, png_Tdiv, png_Tnow)

    ax[0, 1].set(ylabel="# paralog pairs in bin")
    ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=550)
    plt.clf()
    plt.close()


def plot_expository_images(ax, png_Tdiv, png_Tnow):

        im = mpimg.imread(png_Tnow)
        ax[0, 0].imshow(im)
        ax[0, 0].get_xaxis().set_visible(False)
        ax[0, 0].get_yaxis().set_visible(False)
        # Selecting the axis-X making the bottom and top axes False.
        ax[0, 0].tick_params(axis='x', which='both', bottom=False,
                             top=False, labelbottom=False)
        ax[0, 0].tick_params(axis='y', which='both', right=False,
                             left=False, labelleft=False)
        for pos in ['right', 'top', 'bottom', 'left']:
            ax[0, 0].spines[pos].set_visible(False)
        ax[0, 0].set(title="polyploid Ks at T_now")

        #img = Image.open(png_Tdiv)
        #im = plt.imread(get_sample_data(png_Tdiv))
        #img.close()
        im = mpimg.imread(png_Tdiv)
        ax[1, 0].imshow(im)
        ax[1, 0].get_xaxis().set_visible(False)
        ax[1, 0].get_yaxis().set_visible(False)
        ax[1, 0].set(title="ancestral Tc at T_div")
        for pos in ['right', 'top', 'bottom', 'left']:
            ax[1, 0].spines[pos].set_visible(False)


def plot_ks(this_ax, config_used, slim_ks_by_gene, spx_ks_by_gene,
            title, bin_size, xmax, ymax, show_predictions):

    num_slim_genes = len(slim_ks_by_gene)
    num_specks_genes = len(spx_ks_by_gene)

    if not xmax:
        xmax = max(slim_ks_by_gene)
    bins = np.arange(0, xmax, bin_size)

    if len(slim_ks_by_gene) > 0:
        dgx_hist_ys, bins, patches = this_ax.hist(slim_ks_by_gene, bins=bins, facecolor='b', alpha=0.25,
                     label='SLiM Ks by gene\n'
                           + "(" + str(num_slim_genes) + " paralogs in genome)",
                     density=False)

    if len(spx_ks_by_gene) > 0:
        this_ax.hist(spx_ks_by_gene, bins=bins, facecolor='c', alpha=0.25,
                     label='SpecKS Ks by gene\n'
                           + "(" + str(num_specks_genes) + " paralogs in genome)",
                     density=False)

    this_ax.axvline(x=config_used.t_div_as_ks, color='b', linestyle='-.', label="input Tdiv as Ks")
    theoretical_ks_mean_now=config_used.mean_Ks_from_Tc+config_used.t_div_as_ks
    this_ax.axvline(x=theoretical_ks_mean_now, color='r', linestyle='--', label="Expected Ks mean")
    #mean_Ks_from_Nb_string=  "({:.2E})".format(config_used.mean_Ks_from_Nb)
    #this_ax.axvline(x=config_used.mean_Ks_from_Nb,
    #                color='g', linestyle='--', label="Tc due to Nb " + mean_Ks_from_Nb_string )
    if config_used.mig_rate > 0:
        mig_start_as_Ks=config_used.t_div_as_ks- (config_used.Ks_per_YR*config_used.mig_start)
        mig_stop_as_Ks=config_used.t_div_as_ks- (config_used.Ks_per_YR*config_used.mig_stop)
        this_ax.axvline(x=mig_start_as_Ks, color='g', linestyle=':', label="Mig start")
        this_ax.axvline(x=mig_stop_as_Ks, color='g', linestyle='--', label="Mig stop")
        this_ax.axvspan(mig_stop_as_Ks, mig_start_as_Ks, alpha=0.5, color='g')
    add_Ks_annotations(this_ax, config_used, slim_ks_by_gene,dgx_hist_ys, bins, show_predictions)

    if ymax:
        this_ax.set(ylim=[0, ymax])
    this_ax.set(xlim=[0, xmax])
    this_ax.set(xlabel="Ks")
    this_ax.set(title=title)
    this_ax.legend()


def add_Ks_annotations(this_ax, config_used,dgks_Ks_results,dgks_hist_results,
                       bins,show_predictions):

    ks_predictions = ks_modeling.Ks_modeling_predictions(config_used, bins)
    ks_fits=ks_modeling.Ks_modeling_fits(dgks_Ks_results,dgks_hist_results,bins)

    #compare simulated and predicted Ks means
    theoretical_ks_mean_now_as_string="Expected Ks mean ({:.2E})".format(ks_predictions.theoretical_ks_mean_now)
    simulated_ks_mean_now_as_string="Simulated Ks mean ({:.2E})".format(ks_fits.mean_ks_now_from_slim)

    #compare simulated and predicted Ks sigma
    theoretical_sigma_from_kingman_as_string="Kingman Ks sigma ({:.2E})".format(ks_predictions.theoretical_Kingman_sigma_now )
    theoretical_sigma_from_subsampling_genes_as_string = "Theoretical Ks sigma from subsampling ({:.2E})".format(
        ks_predictions.theoretical_sigma_due_to_kingman_from_subsampling_genes)
    simulated_ks_sigma_now_as_string="Simulated Ks sigma ({:.2E})".format(ks_fits.variance_from_slim)

    #plot predictions
    if show_predictions[0]:
        this_ax.plot(bins, ks_predictions.travelling_kingman_ys, label='expectation due to Kingman',
                 linestyle='solid', color='b', alpha=1)
    if show_predictions[1]:
        this_ax.plot(bins,ks_predictions.travelling_gaussian_ys,label='expectation due to CLT',
                 linestyle='solid', color='r',alpha=1)

    #add annotations
    annotation_txt = "\n".join([theoretical_ks_mean_now_as_string,
                                simulated_ks_mean_now_as_string,
                                theoretical_sigma_from_kingman_as_string,
                                theoretical_sigma_from_subsampling_genes_as_string,
                                simulated_ks_sigma_now_as_string])

    this_ax.annotate(annotation_txt, (0, 0), (0, -30), xycoords='axes fraction', textcoords='offset points', va='top')


if __name__ == '__main__':
    unittest.main()
