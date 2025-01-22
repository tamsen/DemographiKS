import glob
import os
import unittest
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import config
import ks_modeling
from data_aggregation.coalescent_plot_aggregation import get_run_time_in_minutes, read_data_csv, plot_mrca
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

        #demographics_TE_run_list=['TE15_m01d03y2025_h10m23s32',
        #                           'TE15_m01d03y2025_h10m23s32','TE07_fix__m01d05y2025_h09m07s41']
        #demographics_TE9_run_list = ['TE10_m12d26y2024_h10m35s41','TE10_m12d26y2024_h10m35s41',
        #    'TE09_m12d26y2024_h09m10s55', 'TE11_m12d26y2024_h10m35s44']

        #demographics_TE9_run_list = ['TE11_m12d26y2024_h10m35s44','TE11_m12d26y2024_h10m35s44',
        #    'TE09_m12d26y2024_h09m10s55']

        #specks_TE9_run_list=['specks_TE07_m12d30y2024_h12m10s15','specks_TE07_m12d30y2024_h12m10s15',
        #                     'specks_TE07_m12d30y2024_h12m10s15']
        #specks_TE9_run_list = ['specks_TE10_m12d31y2024_h09m30s26', 'specks_TE10_m12d31y2024_h09m30s26',
        #                       'specks_TE09_m12d31y2024_h09m10s34',
        #                       'specks_TE11_m12d31y2024_h09m30s22']

        #specks_TE9_run_list = ['specks_TE11_m12d31y2024_h09m30s22', 'specks_TE11_m12d31y2024_h09m30s22',
        #                       'specks_TE09_m12d31y2024_h09m10s34']

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
                                 specks_TE9_run_list, specks_out_path, Ks_per_YR,
                                 xmax_Ks, xmax_Tc, ymax_Ks, suptitle, show_KS_predictions, total_num_genes):
    ymax_Tc = False
    num_runs = len(demographics_TE9_run_list)
    png_out = os.path.join(demographiKS_out_path, "ks_hist_by_TE{0}_test.png".format(run_list_name))
    par_dir = Path(__file__).parent.parent
    image_folder = os.path.join(par_dir, "images")
    png_Tnow = os.path.join(image_folder, 'Ks_now_time_slice.jpg')
    png_Tdiv = os.path.join(image_folder, 'Tdiv_TimeSlice.jpg')
    fig, ax = plt.subplots(2, num_runs, figsize=(20, 10))
    fig.suptitle(suptitle)
    plot_expository_images(ax, png_Tdiv, png_Tnow)
    for i in range(1, num_runs):
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

        plot_ks(ax[0, i], config_used, demographiKS_ks_results, spx_ks_results, config_used.DIV_time_Ge,
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

        theory_mrcas_by_gene=False
        plot_mrca(ax[1, i], slim_mrcas_by_gene, specks_mrcas_by_gene, theory_mrcas_by_gene,
                  dgx_run_duration_in_m, plot_title, config_used.ancestral_Ne,
                  Ks_per_YR, bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc, total_num_genes)

    ax[0, 1].set(ylabel="# paralog pairs in bin")
    ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=550)
    plt.clf()
    plt.close()


def plot_expository_images(ax, png_Tdiv, png_Tnow):

        #img = Image.open(png_Tnow)
        #im = plt.imread(get_sample_data(png_Tnow))
        #img.close()

        #file = open(png_Tnow, "rb")  # Open file in binary read mode
        #im = file.read()
        #file.close()  # Close the file
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


def plot_ks(this_ax, config_used, slim_ks_by_gene, spx_ks_by_gene, t_div,Ne, Ks_per_YR,
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

    modeling_result= ks_modeling.Ks_modeling_result(config_used, bins)
    #Ks_modeling.Ks_modeling_result (Ne, mean_ks_from_Tc, t_div_as_ks, bin_size, bins, config_used, num_slim_genes))

    if show_KS_predictions[0]:
        this_ax.plot(bins,modeling_result.initial_kingman_as_ks,c='red', label='Ks_exp at Tdiv (Kingman assumption)',alpha=0.25)
    if show_KS_predictions[1]:
        this_ax.plot(bins,modeling_result.ks_model_exponential[0:len(bins)],
                 c='gray', label='Ks_exp at Tnow (raw)',alpha=0.25)
    if show_KS_predictions[2]:
        this_ax.plot(bins, modeling_result.ks_model_as_gaussian[0:len(bins)], c='k',
                 label='Ks_exp at Tnow (gaussian)',alpha=0.25)

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
