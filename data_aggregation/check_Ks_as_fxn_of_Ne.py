import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestKsByNe(unittest.TestCase):

    def test_Ks_for_varying_Ne_Tdiv_1000(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'

        # <mutation_rate>1.0e-5</mutation_rate>, <DIV_time_Ge>1000</DIV_time_Ge>


        #no bottleneck
        demographics_TE_run_list = [False, "DGKS_10_10_v2_m01d06y2025_h14m09s31",
                                    "DGKS_100_100_v2_m01d06y2025_h14m09s27",
                                    "DGKS_1000_1000_v2_m01d06y2025_h14m09s23",
                                    "DGKS_5000_5000_m5_BI_40K_m01d13y2025_h10m53s49"]

        demographics_TE_run_list = [False, "DGKS_10_10_v2_m01d06y2025_h14m09s31",
                                    "DGKS_100_100_v2_m01d06y2025_h14m09s27",
                                    "DGKS_1000_1000_v2_m01d06y2025_h14m09s23",
                                    "DGKS_5000_10_m5_BI_40K_m01d13y2025_h13m24s34"]
        #                            "DGKS_5000_2_m5_BI_40K_m01d13y2025_h13m53s40"]

        specks_TE_run_list = [False,'specks_TE10_m01d13y2025_h13m18s28',
                              'specks_TE100_m01d13y2025_h13m17s56',
                              'specks_TE1000_m01d13y2025_h13m17s53',
                              'specks_TE5000_m01d13y2025_h13m18s40']

        # since mutation rate is 1.0e-5
        # we multiply by 1/1.2 since thats syn / total mut rate
        Ks_per_YR = 0.833*10**-5
        burnin_times_in_generations = [2e4, 2e4, 2e4, 2e4, 2e4, 2e4]
        #time_since_DIV = [2000, 2000, 2000, 2000, 2000, 2000]

        bin_sizes_Tc = [80, 80, 80, 400, 800]  # looks good

        xmax_Ks = [0.04,0.04,0.04,0.1,0.2] #[0.01,0.01,0.01,0.1,0.2]#False#0.08  # for mut rate e-5
        bin_sizes_Ks = [0.001, 0.001, 0.001, 0.002, 0.003]

        # xmax_Ks = False # 0.00001  #for mut rate 1.2e-8
        # bin_sizes_Ks = [0.000001, 0.000001, 0.000001, 0.000001, 0.000001]
        xmax_Tc = [2000,2000,2000,20000,40000]
        run_list_num = "_early_DGKS_1000_gen_by_Ne_fast_mut_rate"
        ymax = False

        suptitle = "SLiM Tcoal and Ks\n" + \
                   "Recombination rate = 1.26e-6, mut rate 1.0e-5"
        show_KS_predictions=[False,False,False]
        total_num_genes = 333
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc, burnin_times_in_generations,
                                     demographiKS_out_path, demographics_TE_run_list, run_list_num,
                                     specks_TE_run_list, specks_out_path,Ks_per_YR,
                                     xmax_Ks, xmax_Tc, ymax, suptitle, show_KS_predictions, total_num_genes)

        self.assertEqual(True, True)  # add assertion here

    def test_Ks_for_varying_Ne_Tdiv_75(self):
        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'
        # <mutation_rate>1.0e-5</mutation_rate>, <DIV_time_Ge>75</DIV_time_Ge>
        demographics_TE_run_list = [False, "DGKS_10_10_v2_m01d06y2025_h13m09s31",
                                    "DGKS_100_100_v2_m01d06y2025_h13m09s35",
                                    "DGKS_1000_1000_v2_m01d06y2025_h13m05s18"]
        specks_TE_run_list = [False, False, False, False]
        Ks_per_YR = 10 ** -5
        burnin_times_in_generations = [2e4, 2e4, 2e4, 2e4, 2e4]
        time_since_DIV = [75, 75, 75, 75]
        bin_sizes_Tc = [80, 80, 80, 80, 80]  # looks good
        xmax_Ks = [0.05,0.05,0.05,0.05]  # for mut rate e-5
        bin_sizes_Ks = [0.001, 0.001, 0.001, 0.001, 0.001]
        xmax_Tc = [5000,5000,5000,5000]
        run_list_num = "_early_DGKS_75_gen_by_Ne"
        ymax = False
        suptitle = "SLiM Tcoal and Ks\n" + \
                   "Recombination rate = 1.26e-6, Ne and BI constant"
        show_KS_predictions = [True, False, True]
        total_num_genes = 333
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc, burnin_times_in_generations,
                                     demographiKS_out_path, demographics_TE_run_list, run_list_num,
                                     specks_TE_run_list, specks_out_path, Ks_per_YR,
                                     xmax_Ks, xmax_Tc, ymax, suptitle, show_KS_predictions,
                                     total_num_genes)
        self.assertEqual(True, True)  # a

if __name__ == '__main__':
    unittest.main()
