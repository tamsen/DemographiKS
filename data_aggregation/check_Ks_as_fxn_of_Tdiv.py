import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestGeneLossRates(unittest.TestCase):

    def test_Ks_for_varying_varying_Tdiv_times(self):



        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'
        demographics_TE5_run_list=['TE03_m12d20y2024_h14m26s56',
          'TE05fix__m01d06y2025_h11m17s15','TE07_fix__m01d08y2025_h15m09s22',
           'TE08_m12d24y2024_h09m31s26','TE09_fix_m01d13y2025_h14m13s12']

        specks_TE5_run_list=['specks_TE05_m12d30y2024_h11m50s03',
                             "specks_TE05_m01d14y2025_h09m50s16" ,
        "specks_TE07_m01d14y2025_h09m50s16",
        "specks_TE08_m01d14y2025_h09m50s16" ,
        "specks_TE09_m01d14y2025_h09m50s16" ]


        bin_sizes_Tc = [200,200, 200, 200,200]
        bin_sizes_Ks = [0.0002, 0.0002, 0.0002, 0.0002, 0.0002]
        xmax_Ks = [0.025,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [False,False,False,False,False]
        run_list_name="Ks_for_varying_varying_Tdiv"
        Ks_per_YR = 0.01 * 10**-6 #Mut rate is 1.2 to Ks rate of 1.0 in SpecKS
        ymax_KS = [800,800,800,800,800]
        show_KS_predictions=[False,False,False]
        total_num_genes=3333
        suptitle = "SLiM and SpecKS Ks histograms\n" + \
                                  "Recombination rate = 8e-9, Ne and BI constant"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                          demographiKS_out_path, demographics_TE5_run_list, run_list_name,
                                          specks_TE5_run_list, specks_out_path,Ks_per_YR,
                                     xmax_Ks, xmax_Tc, ymax_KS, suptitle,
                                     show_KS_predictions,total_num_genes)

        self.assertEqual(True, True)  # add assertion here

if __name__ == '__main__':
    unittest.main()

