import os
import unittest
import glob
from matplotlib import pyplot as plt
import config
from data_aggregation.coalescent_plot_aggregation import plot_mrca, read_data_csv, get_run_time_in_minutes, \
    make_Tc_fig_with_subplots


class TestTCoalForSparseness(unittest.TestCase):

    def test_TCoal_with_varying_BI(self):

        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/DGKS_Tc_vs_BI"

        BI_run_list = ["BI2_m12d19y2024_h10m57s18", "BI3_m12d19y2024_h11m02s23",
                       "BI4_m12d19y2024_h11m02s24", "BI5_m12d19y2024_h11m02s27"]
        # "BI7_m12d19y2024_h11m02s31"] <-crashed


        bin_sizes = [10000,10000,10000,10000,10000]
        xmaxs = [160_000,160_000,160_000,160_000,160_000]
        suptitle="SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
                     "Recombination rate = 8e-10, Ne=1e4"
        total_num_genes=[3333 for r in BI_run_list]

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

        RC_run_list_1 = ["RC08_m12d18y2024_h14m30s58",
            "RC09_m12d18y2024_h14m30s54", "RC10_m12d18y2024_h13m36s12",
                       "RC11_m12d18y2024_h14m35s15"]
        RC_run_list_2 = [ "RC12_m12d18y2024_h14m35s17",
                       "RC13_m12d18y2024_h14m35s19", "RC14_m12d18y2024_h14m39s19"]

        bin_sizes = [10000,10000,10000,10000,10000,10000,10000,10000,10000,10000]
        xmaxs = [160_000,160_000,160_000,160_000,160_000,160_000,160_000,160_000,160_000,160_000]
        suptitle= "SLiM Tcoal by gene in ancestral species at Tdiv\n" + \
                     "burnin_time = 5000000 generations"
        total_num_genes=[3333 for r in RC_run_list_1]
        Ks_per_YR=0.00000001

        make_Tc_fig_with_subplots(bin_sizes,
                                  aggragate_output_folder,
                                  RC_run_list_1,"mrcsa_by_recombination_rate_RC_1",
                                  Ks_per_YR,
                                  xmaxs , suptitle,
                                  total_num_genes)

        make_Tc_fig_with_subplots(bin_sizes,
                                  aggragate_output_folder,
                                  RC_run_list_2,"mrcsa_by_recombination_rate_RC_2",
                                  Ks_per_YR,
                                  xmaxs , suptitle,
                                  total_num_genes)

        self.assertEqual(True, True)  # add assertion here


    def test_TCoal_with_varying_Ne(self):

        # make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/DGKS_Tc_vs_Ne"

        #skip this one, to simplify plot "Ne1_m12d18y2024_h16m08s35",
        Ne_run_list = [
                     "Ne2_m12d18y2024_h16m08s38",
                     "Ne3_m12d18y2024_h16m08s40",
                     "Ne4_m12d18y2024_h16m10s52","Ne5_m12d18y2024_h16m10s54",]

        total_num_genes=[10000 for r in Ne_run_list ]
        Ks_per_YR=0.00000001
        num_runs = len(Ne_run_list )

        xmaxs = [False for i in range(0,num_runs)] #160_000
        bin_sizes = [40,400,4000,40000] #effective_population_size[i]/10.0
        suptitle="SLiM Tcoal by gene in ancestral species at Tdiv, by Ne\n" + \
                     "burnin_time = 5000000 generations, recombination rate = 8e-9"

        make_Tc_fig_with_subplots(bin_sizes,
                                  aggragate_output_folder,
                                  Ne_run_list, "mrcsa_by_ancestral_pop_size_Ne",
                                  Ks_per_YR,
                                  xmaxs, suptitle,
                                  total_num_genes)

        self.assertEqual(True, True)  # add assertion here

    def test_TCoal_with_varying_Ge(self):

        # make 1 histogram.
        aggragate_output_folder = "/home/tamsen/Data/DemographiKS_output_from_mesx/DGKS_Tc_vs_GE"

        #"GE4_m12d19y2024_h11m47s58" save if needed...
        GE_run_list =["GE4_m12d19y2024_h11m47s58","GE5_m12d19y2024_h11m47s58",
                      "GE6_m12d19y2024_h11m48s02","GE7_m12d19y2024_h13m30s32"]
         #             "GE8_m12d19y2024_h13m34s13"]
        genome_sizes = [1e4,1e5, 1e6, 1e7,1e8]
        genome_sizes_string = [str(int(g)) for g in genome_sizes]
        bin_sizes = [1000, 1000, 1000, 1000, 1000]
        xmaxs = [75_000, 75_000, 75_000, 75_000, 75_000]
        total_num_genes=[3,33,333,3333]
        Ks_per_YR=0.00000001
        suptitle="SLiM Tcoal by gene in ancestral species at Tdiv,\n"+\
                 "by genome size: " + ", ".join(genome_sizes_string) +\
                     "\nburnin_time = 5000000 generations, recombination rate = 8e-9"

        make_Tc_fig_with_subplots(bin_sizes,
                                  aggragate_output_folder,
                                  GE_run_list, "mrcsa_by_genome_size_GE",
                                  Ks_per_YR,
                                  xmaxs, suptitle,
                                  total_num_genes)

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
