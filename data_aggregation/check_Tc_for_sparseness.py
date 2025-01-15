import os
import unittest
import glob
from matplotlib import pyplot as plt
from pathlib import Path
import config
from data_aggregation.coalescent_plot_aggregation import plot_mrca, read_data_csv, get_run_time_in_minutes
from data_aggregation.ks_plot_aggregations import plot_expository_images


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)  # add assertion here



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

            dgx_run_duration_in_m = get_run_time_in_minutes(dgx_run_path)

        else:
            config_used = False
            demographiKS_ks_results = []
            dgx_run_duration_in_m = 0
            plot_title = "foo - didnt load a config"


        spx_run_name = specks_TE9_run_list[i]



        slim_csv_file = os.path.join(dgx_run_path, "simulated_ancestral_gene_mrcas.csv")
        loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)

        theory_output_file = os.path.join(dgx_run_path, "theoretical_ancestral_gene_mrcas.csv")
        loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
        burnin_times_in_generations=config_used.burnin_time
        plot_title = "Tcoal at Tdiv\nburnin time=" + str(burnin_times_in_generations) + " gen, " \
                     + "Ne=" + str(config_used.ancestral_Ne)

        theory_mrcas_by_gene=False
        plot_mrca(ax[1, i], slim_mrcas_by_gene, False, theory_mrcas_by_gene,
                  dgx_run_duration_in_m, plot_title, config_used.ancestral_Ne,
                  Ks_per_YR, bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc, total_num_genes)

    ax[0, 1].set(ylabel="# paralog pairs in bin")
    ax[1, 1].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=550)
    plt.clf()
    plt.close()

if __name__ == '__main__':
    unittest.main()
