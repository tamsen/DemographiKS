import os
import unittest
import glob
from matplotlib import pyplot as plt
from pathlib import Path
import config
from data_aggregation.coalescent_plot_aggregation import plot_mrca, read_data_csv, get_run_time_in_minutes
from data_aggregation.ks_plot_aggregations import plot_expository_images


class MyTestCase(unittest.TestCase):

    def test_TCoal_with_varying_BI(self):

        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/DGKS_Tc_vs_BI"

        BI_run_list = ["BI2_m12d19y2024_h10m57s18", "BI3_m12d19y2024_h11m02s23",
                       "BI4_m12d19y2024_h11m02s24", "BI5_m12d19y2024_h11m02s27"]
        # "BI7_m12d19y2024_h11m02s31"] <-crashed

        burnin_times_in_generations = [5e2, 5e3, 5e4, 5e5, 5e6]  # , 5e7]#8e-8, <- longer BI times still running..
        #theory_output_file = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
        #                                  BI_run_list[0], "theoretical_ancestral_gene_mrcas.csv")

        bin_sizes = [10000,10000,10000,10000,10000]
        xmaxs = [160_000,160_000,160_000,160_000,160_000]
        suptitle="SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
                     "Recombination rate = 8e-10, Ne=1e4"
        total_num_genes=3333
        ymax_Ks=False

        # Mut rate is 1.2 to Ks rate of 1.0 in SpecKS
        #<mutation_rate>0.000000012</mutation_rate>
        Ks_per_YR=0.00000001
        show_KS_predictions=False
        make_Tc_fig_with_subplots(bin_sizes,
                                  aggragate_output_folder,
                                  BI_run_list,"mrcsa_by_burnin_time_BI",
                                  False, False, Ks_per_YR,
                                  xmaxs , xmaxs , ymax_Ks, suptitle,
                                  show_KS_predictions, total_num_genes)

        self.assertEqual(True, True)  # add assertion here


def make_Tc_fig_with_subplots(bin_sizes_Tc,
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
    fig, ax = plt.subplots(1, num_runs, figsize=(20, 10))
    fig.suptitle(suptitle)
    for i in range(0, num_runs):
        dgx_run_name = demographics_TE9_run_list[i]

        if dgx_run_name:

            dgx_run_path = os.path.join(demographiKS_out_path, dgx_run_name)
            print("dgx_run_path: " +dgx_run_path )
            glob_results=glob.glob(dgx_run_path + '/*.used.xml')
            input_xml_file = glob_results[0]
            config_used = config.DemographiKS_config(input_xml_file)

            dgx_run_duration_in_m = get_run_time_in_minutes(dgx_run_path)

        else:
            config_used = False
            demographiKS_ks_results = []
            dgx_run_duration_in_m = 0
            plot_title = "foo - didnt load a config"



        slim_csv_file = os.path.join(dgx_run_path, "simulated_ancestral_gene_mrcas.csv")
        loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)

        theory_output_file = os.path.join(dgx_run_path, "theoretical_ancestral_gene_mrcas.csv")
        loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
        burnin_times_in_generations=config_used.burnin_time
        plot_title = "Tcoal at Tdiv\nburnin time=" + str(burnin_times_in_generations) + " gen, " \
                     + "Ne=" + str(config_used.ancestral_Ne)

        theory_mrcas_by_gene=False
        plot_mrca(ax[i], slim_mrcas_by_gene, False, theory_mrcas_by_gene,
                  dgx_run_duration_in_m, plot_title, config_used.ancestral_Ne,
                  Ks_per_YR, bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc, total_num_genes)

    ax[0].set(ylabel="# paralog pairs in bin")
    ax[1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=550)
    plt.clf()
    plt.close()

if __name__ == '__main__':
    unittest.main()
