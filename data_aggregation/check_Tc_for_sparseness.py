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


        bin_sizes = [10000,10000,10000,10000,10000]
        xmaxs = [160_000,160_000,160_000,160_000,160_000]
        suptitle="SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
                     "Recombination rate = 8e-10, Ne=1e4"
        total_num_genes=3333

        # Mut rate is 1.2 to Ks rate of 1.0 in SpecKS
        #<mutation_rate>0.000000012</mutation_rate>
        Ks_per_YR=0.00000001
        make_Tc_fig_with_subplots(bin_sizes,
                                  aggragate_output_folder,
                                  BI_run_list,"mrcsa_by_burnin_time_BI",
                                  Ks_per_YR,
                                  xmaxs , suptitle,
                                  total_num_genes)

        self.assertEqual(True, True)  # add assertion here




    def test_TCoal_with_varying_RC(self):

        # make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/DGKS_Tc_vs_RC"

        RC_run_list = ["RC08_m12d18y2024_h14m30s58",
            "RC09_m12d18y2024_h14m30s54", "RC10_m12d18y2024_h13m36s12",
                       "RC11_m12d18y2024_h14m35s15", "RC12_m12d18y2024_h14m35s17",
                       "RC13_m12d18y2024_h14m35s19", "RC14_m12d18y2024_h14m39s19"]

        bin_sizes = [10000,10000,10000,10000,10000,10000,10000,10000,10000,10000]
        xmaxs = [160_000,160_000,160_000,160_000,160_000,160_000,160_000,160_000,160_000,160_000]
        suptitle= "SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
                     "burnin_time = 5000000 generations"
        total_num_genes=3333
        Ks_per_YR=0.00000001

        make_Tc_fig_with_subplots(bin_sizes,
                                  aggragate_output_folder,
                                  RC_run_list,"mrcsa_by_recombination_rate_RC",
                                  Ks_per_YR,
                                  xmaxs , suptitle,
                                  total_num_genes)


        self.assertEqual(True, True)  # add assertion here


    def test_TCoal_with_varying_Ne(self):

        # make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx"

        Ne_run_list = ["Ne1_m12d18y2024_h16m08s35",
                     "Ne2_m12d18y2024_h16m08s38",
                     "Ne3_m12d18y2024_h16m08s40",
                     "Ne4_m12d18y2024_h16m10s52","Ne5_m12d18y2024_h16m10s54",]
        effective_population_size = [1e1,1e2, 1e3, 1e4,1e5]
        total_num_genes=3333
        Ks_per_YR=0.00000001
        num_runs = len(Ne_run_list )

        xmax = False#160_000
        ymax = 2_000
        png_out = os.path.join(aggragate_output_folder, "mrcsa_by_ancestral_pop_size_Ne.png")
        fig, ax = plt.subplots(1, num_runs, figsize=(20, 5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv, by Ne\n" + \
                     "burnin_time = 5000000 generations, recombination rate = 8e-9")

        for i in range(0, num_runs):
            bin_size = effective_population_size[i]/10.0
            run_name = Ne_run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)
            run_duration_in_m = get_run_time_in_minutes(local_output_folder)
            slim_csv_file = os.path.join(local_output_folder, "simulated_ancestral_gene_mrcas.csv")
            theory_output_file = os.path.join(local_output_folder, "theoretical_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)
            loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
            plot_title="pop size (Ne)=" + str(effective_population_size[i])
            plot_mrca(ax[i], slim_mrcas_by_gene, [], theory_mrcas_by_gene,
                      run_duration_in_m,plot_title, bin_size, xmax, ymax)

        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here

    def test_TCoal_with_varying_Ge(self):

        # make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx"

        GE_run_list =["GE4_m12d19y2024_h11m47s58","GE5_m12d19y2024_h11m47s58",
                      "GE6_m12d19y2024_h11m48s02","GE7_m12d19y2024_h13m30s32",
                      "GE8_m12d19y2024_h13m34s13"]
        genome_size = [1e4,1e5, 1e6, 1e7,1e8]
        #bin_sizes= [10,10,100,1000,10000]
        bin_sizes = [1000, 1000, 1000, 1000, 1000]
        num_runs = len(GE_run_list )
        #xmaxes = [20,160,160_0,160_00,160_000]# False#160_000
        xmaxes = [75_000, 75_000, 75_000, 75_000, 75_000]
        ymax = False#2_000
        png_out = os.path.join(aggragate_output_folder, "mrcsa_by_genome_size_GE.png")
        fig, ax = plt.subplots(1, num_runs, figsize=(20, 5))
        fig.suptitle("SLiM Tcoal by gene in ancestral species at Tdiv, by genome size\n" + \
                     "burnin_time = 5000000 generations, recombination rate = 8e-9")

        for i in range(0, num_runs):

            bin_size = bin_sizes[i]
            #ymax = 10**i
            run_name = GE_run_list[i]
            local_output_folder = os.path.join("/home/tamsen/Data/DemographiKS_output_from_mesx",
                                               run_name)
            run_duration_in_m = get_run_time_in_minutes(local_output_folder)
            slim_csv_file = os.path.join(local_output_folder, "simulated_ancestral_gene_mrcas.csv")
            theory_output_file = os.path.join(local_output_folder, "theoretical_ancestral_gene_mrcas.csv")
            loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)
            loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
            plot_title="genome size = " + str(genome_size[i])
            plot_mrca(ax[i], slim_mrcas_by_gene, [], theory_mrcas_by_gene,
                      run_duration_in_m,plot_title, bin_size, xmaxes[i], ymax)

        ax[0].set(ylabel="# genes in bin")
        plt.tight_layout()
        plt.savefig(png_out, dpi=550)
        plt.clf()
        plt.close()

        self.assertEqual(True, True)  # add assertion here




def make_Tc_fig_with_subplots(bin_sizes_Tc,
                                 demographiKS_out_path, demographics_TE9_run_list, run_list_name,
                                 Ks_per_YR,
                                 xmax_Tc, suptitle, total_num_genes):
    ymax_Tc = False
    num_runs = len(demographics_TE9_run_list)
    png_out = os.path.join(demographiKS_out_path, "ks_hist_for_{0}.png".format(run_list_name))
    par_dir = Path(__file__).parent.parent
    fig, ax = plt.subplots(1, num_runs, figsize=(20, 5))
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
            dgx_run_duration_in_m = 0



        slim_csv_file = os.path.join(dgx_run_path, "simulated_ancestral_gene_mrcas.csv")
        loci, slim_mrcas_by_gene = read_data_csv(slim_csv_file)

        theory_output_file = os.path.join(dgx_run_path, "theoretical_ancestral_gene_mrcas.csv")
        loci, theory_mrcas_by_gene = read_data_csv(theory_output_file)
        burnin_times_in_generations=config_used.burnin_time
        plot_title = "Tcoal at Tdiv\nburnin time=" + str(burnin_times_in_generations) + " gen, " \
                     + "Ne=" + str(config_used.ancestral_Ne)

        plot_mrca(ax[i], slim_mrcas_by_gene, False, theory_mrcas_by_gene,
                  dgx_run_duration_in_m, plot_title, config_used.ancestral_Ne,
                  Ks_per_YR, bin_sizes_Tc[i], xmax_Tc[i], ymax_Tc, total_num_genes)

    ax[0].set(ylabel="# genes in bin")
    plt.tight_layout()
    plt.savefig(png_out, dpi=550)
    plt.clf()
    plt.close()

if __name__ == '__main__':
    unittest.main()
