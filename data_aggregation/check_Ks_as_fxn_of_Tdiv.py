import unittest

from data_aggregation.ks_plot_aggregations import make_Tc_Ks_fig_with_subplots


class TestGeneLossRates(unittest.TestCase):

    def test_Ks_for_varying_varying_Tdiv_times(self):

        demographiKS_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'
        #specks_out_path = '/home/tamsen/Data/Specks_output_from_mesx'
        specks_out_path = '/home/tamsen/Data/DemographiKS_output_from_mesx/SPKS_vs_DGKS_Tdiv'
        demographics_TE5_run_list=[False,
          'TE05fix__m01d06y2025_h11m17s15','TE07_fix__m01d08y2025_h15m09s22',
           'TE08_fix_m01d14y2025_h09m17s20','TE09_fix_m01d13y2025_h14m13s12']

        specks_TE5_run_list=[False,
                             "specks_TE05_m01d14y2025_h09m50s16" ,
        "specks_TE07_m01d14y2025_h09m50s16",
        "specks_TE08_m01d14y2025_h09m50s16" ,
        "specks_TE09_m01d14y2025_h09m50s16" ]


        bin_sizes_Tc = [200,200, 200, 200,200]
        bin_sizes_Ks = [0.0002, 0.0002, 0.0002, 0.0002, 0.0002]
        xmax_Ks = [0.025,0.025,0.025,0.025,0.025] #0.001  # max(demographiKS_ks_results)
        xmax_Tc = [False,False,False,False,False]
        run_list_name="Ks_for_varying_varying_Tdiv"
        ymax_KS = [800,800,800,800,800]
        show_KS_predictions=[False,False,False]
        suptitle = "SLiM and SpecKS Ks histograms\n" + \
                                  "Recombination rate = 8e-9, Ne and BI constant"
        make_Tc_Ks_fig_with_subplots(bin_sizes_Ks, bin_sizes_Tc,
                                          demographiKS_out_path, demographics_TE5_run_list, run_list_name,
                                          specks_TE5_run_list, specks_out_path,
                                     xmax_Ks, xmax_Tc, ymax_KS, suptitle,
                                     show_KS_predictions)

        self.assertEqual(True, True)  # add assertion here

if __name__ == '__main__':
    unittest.main()

