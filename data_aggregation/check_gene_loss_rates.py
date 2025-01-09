import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestGeneLossRates(unittest.TestCase):

    def test_gene_loss_for_varying_Tdiv_times(self):



        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'
        demographics_TE5_run_list=['TE03_m12d20y2024_h14m26s56',
          'TE05fix__m01d06y2025_h11m17s15','TE07_fix__m01d08y2025_h15m09s22',
           'TE08_m12d24y2024_h09m31s26','TE09_m12d26y2024_h09m10s55']

        specks_TE5_run_list=['specks_TE05_m12d30y2024_h11m50s03','specks_TE05_m12d30y2024_h11m50s03',
                             'specks_TE07_m12d30y2024_h12m10s15',
                             'specks_TE08_m12d30y2024_h12m10s13','specks_TE09_m12d30y2024_h12m10s11']

        burnin_times_in_generations = [5e7, 5e7, 5e7, 5e7, 5e7]
        time_since_DIV = [1000, 10000, 100000, 500000, 1000000]

        bin_sizes_Tc = [200,200, 200, 200,200]
        bin_sizes_Ks = [0.0002, 0.0002, 0.0002, 0.0002, 0.0002]
        xmax_Ks = 0.025#0.001  # max(demographiKS_ks_results)
        xmax_Tc = False
        run_list_num="_5to9_vary_Tdiv_gene_shedding_fix"
        Ks_per_YR = 0.01 * 10**-6 #Mut rate is 1.2 to Ks rate of 1.0 in SpecKS
        ymax = False
        show_KS_predictions=[False,False,False]
        suptitle = "SLiM and SpecKS Ks histograms\n" + \
                                  "Recombination rate = 1.26e-6, Ne and BI constant"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc, burnin_times_in_generations,
                                          demographiKS_out_path, demographics_TE5_run_list, run_list_num,
                                          specks_TE5_run_list, specks_out_path, time_since_DIV,
                                            Ks_per_YR, xmax_Ks, xmax_Tc, ymax, suptitle,
                                     show_KS_predictions)



        self.assertEqual(True, True)  # add assertion here

if __name__ == '__main__':
    unittest.main()

